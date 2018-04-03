#' Calculate protein or peptide labeling rate
#'
#' Calculate protein or peptide labeling rate by non-linear least squares fitting to exponential equation \eqn{f_{labeled} = 1 - e^{-k\cdot t}} where \eqn{k} is the protein/peptide labeling rate (in units of reciprocal time, depending on what the time column is recorded in). Forces fit through the origin (i.e. assumes no labeling at time point 0) and thus explicitly ignores 0 time points.
#' @param data the data with frac_lab/frac_ulab calculated (\link{tor_calculate_labeled_fraction})
#' @param time_col the name of the column that holds the time, assumes "hours" by default.
#' @param min_num_timepoints what is the minimum number of time points to try to fit a curve? 0 time point do not count towards making this minimum
#' @param combine_peptides whether to combine all peptides for a fit (i.e. fit the whole protein), or to calculate a fit for each peptide. Default is to combine across peptides.
#' @return data frame with one line for each peptide or protein (depending on \code{combine_peptides}) and the following added columns:
#' \itemize{
#' \item{\code{nested_data}: }{all the data for the protein or peptide}
#' \item{\code{num_peptides}: }{the number of peptides combined for each protein (always 1 if \code{combine_peptides = FALSE})}
#' \item{\code{num_timepoints}: }{the number of unique time points that went into each fit - does not include any 0 time points}
#' \item{\code{num_datapoints}: }{the number of individual data points that went into each fit}
#' \item{\code{enough_data}: }{whether there was enough data to fit a curve at all - i.e. whether there were at least 2 time points}
#' \item{\code{fit_error}: }{whether the non-linear squares fit succeeded or failed}
#' \item{\code{label_rate}: }{the estimated label rate, in units of reciprocal time depending on the \code{time_col}}
#' \item{\code{label_rate_se}: }{the standard error of the estimated label rate based on the least squares fit}
#' \item{\code{fit_rse}: }{the residual standard error of the fit}
#' \item{\code{fit}: }{the actual regression model fit}
#' }
#' @export
tor_calculate_label_rate <- function(data, time_col = "hours", min_num_timepoints = 2, combine_peptides = TRUE, quiet = FALSE) {

  # safety checks for data
  if (missing(data)) stop("need to supply a data set", call. =FALSE)
  if (!is.data.frame(data)) {
    glue("wrong data type supplied: {class(data)[1]}") %>% stop(call. = FALSE)
  }


  # safety checks for specific variables
  columns <- c("frac_lab", "prot_id", "peptide_seq", time_col)
  missing <- setdiff(columns, names(data))
  if (length(missing) > 0) {
    glue("missing column(s) '{collapse(missing, sep = ', ', last = ' and ')}' do(es) not exist in the dataset") %>% stop(call. = FALSE)
  }

  # specify grouping columns: group by protein only or by individual peptides
  group_columns <- "prot_id"
  if ("uniprot_id" %in% names(data))
    group_columns <- c(group_columns, "uniprot_id")
  if (!combine_peptides)
    group_columns <- c(group_columns, "peptide_seq")

  # nest and calculate number of peptides
  data <-
    data %>% nest(!!!map(group_columns, ~quo(-!!sym(.x))), .key = "nested_data") %>%
    mutate(
      num_peptides = map_int(
        nested_data,
        ~if(combine_peptides) { length(unique(.x$peptide_seq)) } else { 1L }
      )
    )

  # information for user
  if (!quiet) {
    type <- if(combine_peptides) 'proteins' else 'peptides'
    glue("Info: processing data for {nrow(data)} ",
         "{type}, this may take a few seconds... ") %>%
      message(appendLF = FALSE)
  }

  # make sure to catch non-convergent NLS
  safe_nls <- safely(nls)

  #perform curve generation and curve fitting
  data <- data %>% mutate(
    nested_data = map(nested_data, ~mutate(.x, ..time.. = !!sym(time_col))),
    filtered_data = map(nested_data, ~filter(.x, ..time.. > 0)),
    num_timepoints = map_int(filtered_data, ~length(unique(.x$..time..))),
    num_datapoints = map_int(filtered_data, ~length(.x$..time..)),

    enough_data = num_timepoints >= min_num_timepoints,
    lm_fit = map2(filtered_data, enough_data, ~if(.y){ lm(-log(1 - frac_lab) ~ ..time.. - 1, data = .x) } else {(NULL)}),
    lm_coefficients = map(lm_fit, tidy),
    lm_summary = map(lm_fit, glance),
    # # use linear fit as a starting estimate for non-linear fit
    nls_safe_fit = map2(filtered_data, lm_coefficients, ~safe_nls(frac_lab ~ 1 - exp(-k_synth * ..time..), start = list(k_synth = .y$estimate), data = .x)),
    nls_error = map_lgl(nls_safe_fit, ~!is.null(.x$error)),
    nls_fit = map2(nls_safe_fit, enough_data, ~if(.y){.x$result} else {NULL}),
    nls_coefficients = map(nls_fit, tidy),
    nls_summary = map(nls_fit, glance),
    nested_data = map(nested_data, ~select(.x, -..time..))
  ) %>%
    # don't need nls safe fit
    select(-nls_safe_fit, -filtered_data)

  if (!quiet) {
    good_fit_n  <- nrow(filter(data, !nls_error & enough_data))
    bad_fit_n_error <- nrow(filter(data, nls_error))
    bad_fit_n_data <- nrow(filter(data, !enough_data))
    if (bad_fit_n_data == 0 && bad_fit_n_error == 0) {
      glue("all of the {type} could be fit to a labeling curve.") %>%
        message()
    } else {
      glue("{good_fit_n} of the {type} could be fit to a labeling curve",
           "{if (bad_fit_n_data > 0) paste0(', ', bad_fit_n_data, ' did not have enough time points') else ''}",
           "{if (bad_fit_n_error > 0) paste0(', ', bad_fit_n_error, ' could not be fit succesfully') else ''}") %>%
        message()
    }
  }

  # basic columns
  core_columns <- c(group_columns, "nested_data", "num_peptides", "num_timepoints", "num_datapoints", "enough_data", fit_error = "nls_error")

  # data without fit
  bad_data <- filter(data, enough_data == FALSE | nls_error == TRUE) %>%
    select(!!!map(core_columns, sym))

  # data with fit
  good_data <- filter(data, enough_data == TRUE, nls_error == FALSE)

  if(nrow(good_data) > 0){
    good_data <- good_data %>%
      unnest(nls_summary, .drop = FALSE) %>%
      unnest(nls_coefficients, .drop = FALSE) %>%
      select(!!!map(core_columns, sym), label_rate = estimate,
             label_rate_se = std.error, fit_rse = sigma,
             fit = nls_fit, fit_error = nls_error)
    final_data <-bind_rows(good_data, bad_data)
  }else{
    final_data <- bad_data
  }

  return(final_data)
}

#' Calculate growth parameters
#'
#' Calculate the generation time and growth rate in units of time and reciprocal time, respectively (same units as the denominator of the flow rate) as well as . The flow rate and volume must have the same mass units (e.g. g, kg, mL, or L). Calculates the standard errors as well if error estimates for flow rate and volume are provided.
#'
#' @param flow_rate the flow rate (mass/time)
#' @param flow_rate_se the standard error of the flow rate (should be same units as flow_rate)
#' @param volume the reactor volume (mass)
#' @param volume_se the standard error of the volume (same units as volume)
#' @return a data frame with growth rate and its error as well as doubling time and its error
#' @export
tor_calculate_growth_params <- function(flow_rate, volume, flow_rate_se = NA_real_, volume_se = NA_real_) {

  if (missing(flow_rate)) stop("no flow rate supplied", call. = FALSE)
  if (missing(volume)) stop("no volume supplied", call. = FALSE)

  data_frame(
    growth_rate = flow_rate / volume,
    growth_rate_se = abs(growth_rate) * sqrt( (flow_rate_se/flow_rate)^2 + (volume_se/volume)^2),
    gen_time = log(2) / growth_rate,
    gen_time_se = abs(gen_time * growth_rate_se/growth_rate)
  )
}


#' Calculate degradation rate
#'
#' Calculate degradataion and dissipation, including propagated error. Protein degradation is calculated as the labeling rate minus the growth rate (deg_rate = label_rate - growth_rate). Protein dissipation is calculated as the degradation rate divided by the labeling rate, and is displayed as a percent (dissipation = (deg_rate / label_rate)*100).
#'
#' @param data data with label rate calculated
#' @param growth_rate growth rate from experiment
#' @param growth_rate_se standard error of growth rate from experiment
#' @export
tor_calculate_degradation_dissipation <- function(data, growth_rate, growth_rate_se = 0, quiet = FALSE) {

  # safety checks for data
  if (missing(data)) stop("need to supply a data set", call. =FALSE)
  if (!is.data.frame(data)) {
    glue("wrong data type supplied: {class(data)[1]}") %>% stop(call. = FALSE)
  }
  if (missing(growth_rate)) stop("need to supply experiment specific growth rate", call. =FALSE)


  #safety check for required variables (need label rate, label reate error columns)
  if (!"label_rate" %in% names(data)){
    stop ("Label rate column does not exist in the dataset. Run function tor_calculate_label_rate first.") %>% stop(call. = FALSE)
  }
  if (!"label_rate_se" %in% names(data)){
    stop ("Label rate error column does not exist in the dataset. Run function tor_calculate_label_rate first.") %>% stop(call. = FALSE)
  }

  data<- data %>%
    mutate(deg_rate = label_rate - growth_rate,
           deg_rate_se = sqrt(label_rate_se^2 + growth_rate_se^2),
           dissipation = (deg_rate / label_rate)*100 , #label_rate or growth rate
           dissipation_se = abs(dissipation * sqrt((deg_rate_se/deg_rate)^2 + (label_rate_se/label_rate)^2)),
           growth_rate = growth_rate,
           growth_rate_se = growth_rate_se
    )

  if (!quiet) {
    glue("Info: calculated degradation rate and dissipation for {nrow(filter(data, !is.na(label_rate)))} records") %>%
      message()
  }

  return(data)

}
