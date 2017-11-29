

# convenience function for guessing metadata time frame from original sample names
# not exported because it might be a very special user case
guess_metadata <- function(smv_data, save_file = "metadata.xlsx") {
  regexp_pattern <- "^(\\d{1,2})(\\d{2})([ap])-?(\\d)?$"
  metadata <- svm_data %>% select(sample) %>% unique() %>%
    mutate(
      matches_pattern = str_detect(sample, regexp_pattern),
      hour = str_match(sample, regexp_pattern) %>% { .[,2] },
      minutes = str_match(sample, regexp_pattern) %>% { .[,3] },
      ampm = str_match(sample, regexp_pattern) %>% { .[,4] },
      # not sure that's what it is but probably
      replicate = str_match(sample, regexp_pattern) %>% { .[,5] } %>% { ifelse(is.na(.), "0", .) },
      time_str = glue("2017-11-17 {hour}:{minutes}{ampm}m"),
      time = as.POSIXct(time_str, format = "%Y-%m-%d %I:%M%p"),
      daily_hours = as.numeric(time - time[1], "hours"),
      day_switch = cumsum(c(0, diff(daily_hours)) < 0),
      hours = daily_hours + 24*day_switch
    ) %>%
    select(sample, replicate, hours)

  # save in the provided excel file
  metadata %>% openxlsx::write.xlsx(save_file)

  return(metadata)
}
