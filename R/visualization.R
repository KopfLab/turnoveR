#' Plot label rate histogram
#' @description creates histogram of labeing rate
#' @param data data with label rate calculated, including growth rate
plot_label_rate_hist <- function(data) {

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
    ggplot() +
    aes(x = label_rate) +
    geom_vline(data = function(x) unique(x %>% select(growth_rate)),
               mapping = aes(xintercept = growth_rate)) +
    geom_histogram(position = "dodge", bins = 80)

  ggplot2::last_plot() %+% aes(fill = factor(growth_rate))
}




#' plot label rate error
#' @description visualize the label rate and residual standard error
#' @param data data with label rate calculated, including growth rate standard errors
plot_label_rate_error <- function(data) {

  # safety checks for data
  if (missing(data)) stop("need to supply a data set", call. =FALSE)
  if (!is.data.frame(data)) {
    glue("wrong data type supplied: {class(data)[1]}") %>% stop(call. = FALSE)
  }

  # safety checks for specific variables
  columns <- c("growth_rate", "label_rate", "fit_rse", "label_rate_se")
  missing <- setdiff(columns, names(data))
  if (length(missing) > 0) {
    glue("column '{collapse(missing, sep = ', ', last = ' and ')}' does not exist in the dataset") %>% stop(call. = FALSE)
  }

  data %>%
    ggplot() +
    aes(x = label_rate, y = fit_rse, label = protein) + #, color = factor(num_datapoints)) +
    geom_errorbarh(map = aes(xmin = label_rate - label_rate_se, xmax = label_rate + label_rate_se)) +
    geom_vline(data = function(x) unique(x %>% select(growth_rate)),
               mapping = aes(xintercept = growth_rate)) +
    geom_point() +
    scale_y_continuous(labels = function(x) str_c(100*x)) +
    #scale_color_gradient2(low = "blue", high = "red", mid = "green") +
    #theme(legend.position = "none") +
    labs(y = "residual standard error of fit [% labelled]", x = "label rate [1/time]")

}



#' plot disippation curves
#' @description visualuze the dissipation curve of peptides or proteins
#' @param data data with label rate (and error), growth rate (and error)
#' @param grouping by protein or peptide for visualization ??
#' @param plot_number designate number of curves you want to see
plot_labeling_curves <- function(data, grouping = protein, plot_number = 6) {

  # safety checks for data
  if (missing(data)) stop("need to supply a data set", call. =FALSE)
  if (!is.data.frame(data)) {
    glue("wrong data type supplied: {class(data)[1]}") %>% stop(call. = FALSE)
  }
  # safety checks for specific variables
  columns <- c("nested_data", "label_rate", "growth_rate", "fit")
  missing <- setdiff(columns, names(data))
  if (length(missing) > 0) {
    glue("column '{collapse(missing, sep = ', ', last = ' and ')}' does not exist in the dataset") %>% stop(call. = FALSE)
  }

  # safety checks for specific variables in nested data
  unnested_checks <- data$nested_data[[1]]
  columns <- c("hours", "frac_lab" )
  missing <- setdiff(columns, names(unnested_checks))
  if (length(missing) > 0) {
    glue("columns '{collapse(missing, sep = ', ', last = ' and ')}' do not exist in the dataset") %>% stop(call. = FALSE)
  }

  plot_data <- data %>%
    sample_n(size = min(nrow(data), plot_number)) %>%
    # introduce a unique row number to group lines together easily (no matter if it's proteins or peptides)
    mutate(`#` = row_number()) %>%
    # calculate curves
    mutate(
      curves = pmap(list(data = nested_data, label_rate = label_rate, growth_rate = growth_rate),
                    function(data, label_rate, growth_rate) {
                      data_frame(
                        hours = seq(0, max(data$hours), length.out = 20),
                        fit =  1 - exp(-label_rate * hours),
                        expected = 1 - exp(-growth_rate * hours)
                      ) %>%
                        gather(curve, frac_lab, -hours)
                    })
    )

  plot_data %>%
    ggplot() +
    aes(x = hours, y = frac_lab, group = `#`) +
    geom_point(data = function(x) unnest(x, nested_data)) +
    geom_line(data = function(x) unnest(x, curves), mapping = aes(color = NULL, group = curve, linetype = curve)) +
    aes(color = protein)+
    facet_wrap(~protein)
}
