#' Plot label rate histogram
#'
#' Plots the label rate histogram together with the growth rate as a reference point. Supports plotting of multiple datasets by coloring according to growth rate.
#'
#' @param data data with label rate calculated, including growth rate
#' @export
tor_plot_label_rate_hist <- function(data) {

  # safety checks for data
  if (missing(data)) stop("need to supply a data set", call. =FALSE)
  if (!is.data.frame(data)) {
    glue("wrong data type supplied: {class(data)[1]}") %>% stop(call. = FALSE)
  }

  # safety checks for specific variables
  columns <- c("growth_rate", "label_rate")
  missing <- setdiff(columns, names(data))
  if (length(missing) > 0) {
    glue("column '{collapse(missing, sep = ', ', last = ' and ')}' does not exist in the dataset") %>% stop(call. = FALSE)
  }

  data %>%
    filter(!is.na(label_rate)) %>%
    ggplot() +
    aes(x = label_rate, fill = factor(signif(growth_rate, 2))) +
    geom_histogram(position = "dodge", bins = 80) +
    geom_vline(data = function(x) unique(x %>% select(growth_rate)),
               mapping = aes(xintercept = growth_rate)) +
    scale_y_continuous(expand = c(0, 0)) +
    labs(x = "label rate", fill = "growth rate") +
    theme_bw()
}

#' Plot label rate error
#'
#' Visualize the label rate and residual standard error. Displays the residual standard error (i.e. the normalized error of the fit) for each calculated label rate to easily evaluate what quality cutoffs might be useful. Plots the growth rate as a reference and the individual standard errors of the parameters from the fit. Supports plotting of multiple datasets by coloring according to growth rate.
#'
#' @param data data with label rate calculated, including growth rate standard errors
#' @export
tor_plot_label_rate_error <- function(data) {

  # safety checks for data
  if (missing(data)) stop("need to supply a data set", call. =FALSE)
  if (!is.data.frame(data)) {
    glue("wrong data type supplied: {class(data)[1]}") %>% stop(call. = FALSE)
  }

  # safety checks for specific variables
  columns <- c("prot_id", "growth_rate", "label_rate", "fit_rse", "label_rate_se")
  missing <- setdiff(columns, names(data))
  if (length(missing) > 0) {
    glue("column '{collapse(missing, sep = ', ', last = ' and ')}' does not exist in the dataset") %>% stop(call. = FALSE)
  }

  data %>%
    filter(!is.na(label_rate)) %>%
    mutate(`growth rate` = factor(signif(growth_rate, 2))) %>%
    ggplot() +
    aes(text = prot_id, x = label_rate, y = fit_rse,
        color = `growth rate`) +
    geom_errorbarh(map = aes(xmin = label_rate - label_rate_se, xmax = label_rate + label_rate_se)) +
    geom_vline(data = function(x) unique(x %>% select(growth_rate)),
               mapping = aes(xintercept = growth_rate)) +
    geom_point() +
    scale_y_continuous(labels = function(x) str_c(100*x)) +
    expand_limits(y = 0) +
    #scale_color_gradient2(low = "blue", high = "red", mid = "green") +
    #theme(legend.position = "none") +
    labs(y = "residual standard error of fit [% labelled]", x = "label rate [1/time]") +
    theme_bw()
}

#' Plot labeling curves
#'
#' Visualizes the labeling fits for peptides or proteins. Selects a random subset from the provided data and plots a maximum of \code{plot_number} of plots.
#'
#' @param data data with label rate (and error), growth rate (and error)
#' @param grouping by protein or peptide for visualization ??
#' @param time_col the name of the column that holds the time, assumes "hours" by default.
#' @param plot_number designate the maximum number of curves you want to plot. Use \code{plot_number = NULL} to plot all - be careful with this option, the plot might take a long time to render with a larger number of proteins or peptides.
#' @param random_seed specify this to get reproducible selection of "random" plots
#' @export
tor_plot_labeling_curves <- function(data, grouping = protein, time_col = "hours", plot_number = 6, random_seed = random()) {

  # safety checks for data
  if (missing(data)) stop("need to supply a data set", call. =FALSE)
  if (!is.data.frame(data)) {
    glue("wrong data type supplied: {class(data)[1]}") %>% stop(call. = FALSE)
  }
  # safety checks for specific variables
  columns <- c("prot_id", "nested_data", "label_rate", "growth_rate", "fit")
  missing <- setdiff(columns, names(data))
  if (length(missing) > 0) {
    glue("column '{collapse(missing, sep = ', ', last = ' and ')}' does not exist in the dataset") %>% stop(call. = FALSE)
  }

  # safety checks for specific variables in nested data
  unnested_checks <- data$nested_data[[1]]
  columns <- c(time_col, "frac_lab")
  missing <- setdiff(columns, names(unnested_checks))
  if (length(missing) > 0) {
    glue("columns '{collapse(missing, sep = ', ', last = ' and ')}' do not exist in the dataset") %>% stop(call. = FALSE)
  }

  # check grouping mode
  by_peptide <- "peptide_seq" %in% names(data)

  # random seed
  random <- function() sample(.Random.seed, 1)
  set.seed(random_seed)

  # generate plot data df
  plot_data <- data %>%
    {
      if (is.null(plot_number))
        .
      else
        sample_n(., size = min(nrow(data), plot_number))
    } %>%
    # introduce a unique row number to group lines together easily (no matter if it's proteins or peptides)
    mutate(`#` = row_number()) %>%
    # calculate curves
    mutate(curves =
             pmap(list(
               data = nested_data,
               label_rate = label_rate,
               growth_rate = growth_rate
             ),
             function(data, label_rate, growth_rate) {
               data_frame(
                 !!time_col := seq(0, max(data[[time_col]]), length.out = 20),
                 fit =  1 - exp(-label_rate * hours),
                 expected = 1 - exp(-growth_rate * hours)
               ) %>%
                 gather(curve, frac_lab, !!quo(-!!sym(time_col)))
             })
    )

  p <- plot_data %>%
    ggplot() +
    aes_q(x = sym(time_col), y = sym("frac_lab"), group = sym("#")) +
    suppressWarnings(
      # warnings need to be surpressed because of the text aesthetic for plotly
      geom_point(data = function(x) unnest(x, nested_data) %>% filter(!is.na(frac_lab)),
                 mapping = aes(text = peptide_seq))
    ) +
    geom_line(data = function(x) unnest(x, curves) %>% filter(!is.na(frac_lab)),
              mapping = aes(color = NULL, group = curve, linetype = curve)) +
    theme_bw() +
    scale_y_continuous("fraction labeled", labels = scales::percent)

  # coloring and facetting
  if (by_peptide)
    p <- p + aes(color = peptide_seq) + facet_wrap(~peptide_seq)
  else
    p <- p + aes(color = prot_id) + facet_wrap(~prot_id)

  return(p)
}
