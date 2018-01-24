#' Calculate protein or peptide labeling rate
#' @description calculate protein or peptide labeling rate from svm_data
#' @param data the svm_data with fraclab/fraculab calculated \link{calculate_fraculab}
calculate_label_rate <- function(data, combine_peptides = TRUE, quiet = FALSE) {

  # fit each peptide to exponential equation y = Ae^dx (timepoint, frac_lab)
  #extract d value for each peptide, store in new column

  #add checks for data (see safety checks dataformat, colums, etc)
  #add check for hours, proteins, isopep... fac_lab


  # make sure to catch non-convergent NLS
  safe_nls <- safely(nls)

  if(combine_peptides)
    data <-
      data %>%
      nest(-protein, .key = "nested_data")
  else
    data<-
      data %>%
      nest(-protein, -isopep, .key = "nested_data")

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
    nls_safe_fit = map2(nested_data, lm_coefficients, ~safe_nls(frac_lab ~ 1 - exp(-k_synth * hours), start = list(k_synth = .y$estimate), data = .x)),
    nls_error = map_lgl(nls_safe_fit, ~!is.null(.x$error)),
    nls_fit = map2(nls_safe_fit, enough_data, ~if(.y){.x$result} else {NULL}),
     nls_coefficients = map(nls_fit, tidy),
     nls_summary = map(nls_fit, glance)
  ) %>%
    # don't need nls safe fit
    select(-nls_safe_fit)

  if (!quiet) {
    glue("{nrow(filter(data, !nls_error))} of the {if(combine_peptides)'proteins' else 'peptides'} could be fit to a labeling curve.") %>%
      message()
  }
  return(data)
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
