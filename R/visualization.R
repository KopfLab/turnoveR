#' Plot label rate histogram
#' @description calculate kdeg (previously performed in Igor)
#' @param data data with label rate calculated
#' @param growth_rate growth rate from experiment
#' @param growth_rate_se standard error of growth rate from experiment
plot_label_rate_hist <- function(data) {

#check for data frame
#check for growth_rate, label_rate

data %>%
  ggplot() +
  aes(x = label_rate) +
  geom_vline(data = function(x) unique(x %>% select(growth_rate)),
             mapping = aes(xintercept = growth_rate)) +
  geom_histogram(position = "dodge", bins = 80)

ggplot2::last_plot() %+% aes(fill = factor(growth_rate))
}
