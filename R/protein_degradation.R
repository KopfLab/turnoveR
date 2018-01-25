#' Calculate protein or peptide labeling rate
#' @description calculate protein or peptide labeling rate from svm_data (fit each peptide to exponential equation y = Ae^dx)
#' @param data the svm_data with fraclab/fraculab calculated \link{calculate_fraculab}
calculate_label_rate <- function(data, combine_peptides = TRUE, quiet = FALSE) {

  # safety checks for data, specific variables
  if (missing(data)) stop("need to supply a data set", call. =FALSE)
  if (!is.data.frame(data)) {
    glue("wrong data type supplied: {class(data)[1]}") %>% stop(call. = FALSE)
  }

  columns <- c("frac_lab", "protein", "isopep", "hours")
  missing <- setdiff(columns, names(data))
  if (length(missing) > 0) {
    glue("columns '{collapse(missing, sep = ', ', last = ' and ')}' do not exist in the dataset") %>% stop(call. = FALSE)
  }

  # make sure to catch non-convergent NLS
  safe_nls <- safely(nls)

  if(combine_peptides){
    data <-
      data %>%
      nest(-protein, .key = "nested_data") %>%
      mutate(num_isopep = map_int(nested_data, ~length(unique(.x$isopep))))
    group_columns = "protein"
  }else{
    data<-
      data %>%
      nest(-protein, -isopep, .key = "nested_data") %>%
      mutate(num_isopep = 1)
    group_columns = c("protein", "isopep")
  }

  # information for user
  if (!quiet) {
    number <- nrow(data)
    glue("Info: Processing data for {number} {if(combine_peptides)'proteins' else 'peptides'}...this may take a few seconds") %>%
      message()
  }

  #perform curve generation and curve fitting
  data <- data %>% mutate(
    filtered_data = map(nested_data, ~filter(.x, hours > 0)),
    num_timepoints = map_int(filtered_data, ~length(unique(.x$hours))),
    num_datapoints = map_int(filtered_data, ~length(.x$hours)),

    enough_data = num_timepoints > 1,
    lm_fit = map2(filtered_data, enough_data, ~if(.y){ lm(-log(1 - frac_lab) ~ hours - 1, data = .x) } else {(NULL)}),
    lm_coefficients = map(lm_fit, tidy),
    lm_summary = map(lm_fit, glance),
    # # use linear fit as a starting estimate for non-linear fit
    nls_safe_fit = map2(filtered_data, lm_coefficients, ~safe_nls(frac_lab ~ 1 - exp(-k_synth * hours), start = list(k_synth = .y$estimate), data = .x)),
    nls_error = map_lgl(nls_safe_fit, ~!is.null(.x$error)),
    nls_fit = map2(nls_safe_fit, enough_data, ~if(.y){.x$result} else {NULL}),
    nls_coefficients = map(nls_fit, tidy),
    nls_summary = map(nls_fit, glance)
  ) %>%
    # don't need nls safe fit
    select(-nls_safe_fit, -filtered_data)

  if (!quiet) {
    glue("{nrow(filter(data, !nls_error))} of the {if(combine_peptides)'proteins' else 'peptides'} could be fit to a labeling curve.") %>%
      message()
  }

  bad_data <-filter(data, enough_data == FALSE | nls_error == TRUE) %>%
    select(!!!group_columns, num_isopep, nested_data, num_timepoints, num_datapoints, enough_data, nls_error)

  good_data <- filter(data, enough_data == TRUE, nls_error == FALSE)
  if(nrow(good_data) > 0){
    good_data <- good_data %>%
      unnest(nls_summary, .drop = FALSE) %>%
      unnest(nls_coefficients, .drop = FALSE) %>%
      select(!!!group_columns, num_isopep, nested_data, label_rate = estimate, label_rate_se = std.error, fit_rse = sigma,
             num_timepoints, num_datapoints, enough_data, nls_error)
    final_data <-bind_rows(good_data, bad_data)
  }else{
    final_data <- bad_data
  }








  #keep/return value of parameter("d"), stderror, rse of overall fit, number timepoints, number datapoints
  #also keep protein/isopep nested data for plotting
  #keep enough_data and nls_error for failed proteins/isopeps

  return(final_data)
}


#' Calculate degradation rate, chisquared value
#' @description calculate kdeg and chisquared (Step 12a-c in workflow, previously performed in Igor)
#' @param fraculab_data or fraculab_data_clean
calc_pep_degrate <- function(fraculab_data_clean) {

  # fit each peptide to exponential equation y = Ae^dx (timepoint, frac_lab)
  #extract d value for each peptide, store in new column
  #calculate chisquared value, store in new column

  # output example: return(deg_rate_data)
}
