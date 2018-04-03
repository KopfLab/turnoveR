# general utilities ====

#' Filter label rate fits
#'
#' Filter protein or peptide label rate fit records for quality, visualization or focus on problematic records. This function is typically called some time after \link{tor_calculate_label_rate} or \link{tor_calculate_degradation_dissipation} before processing the label rate fits further. Use the \code{condition} parameter to set filtering conditions. Use the \code{select} parameter to keep only the most informative columns.
#'
#' @param dt data label rates calculated
#' @param condition filtering condition. Includes all records by default.
#' @param select which columns to select - use \code{c(...)} to select multiple, supports all \link[dplyr]{select} syntax including renaming columns. Includes all columns by default.
#' @export
tor_filter_label_rate_fits <- function(dt, condition = TRUE, select = everything(), quiet = FALSE) {

  # safety checks
  if (missing(dt)) stop("no data table supplied", call. = FALSE)

  # filtering and column selection
  filter_quo <- enquo(condition)
  dt_out <- dt %>%
    filter(!!filter_quo) %>%
    dplyr::select(!!enquo(select))


  # info
  if(!quiet) {
    glue("Info: fetching {nrow(dt_out)} of {nrow(dt)} ",
         "entries based on filter condition '{quo_text(filter_quo)}'") %>%
      message(appendLF = FALSE)
    if (ncol(dt_out) != ncol(dt))
      glue(", keeping {ncol(dt_out)} of {ncol(dt)} columns.") %>% message()
    else
      message()
  }

  return(dt_out)

}

#' Remove complex data
#'
#' Convenience function to remove nested and list columns in a data table (e.g. in preparation for printing to console or RMarkdown).
#' @export
tor_remove_list_columns <- function(dt) {

  if (missing(dt)) stop("no data table supplied", call. = FALSE)

  list_cols <- dt %>% map_lgl(is_list)
  dt[!list_cols]
}

# metadata =====

#' Add metadata
#'
#' @param data the mass spec data
#' @param metadata a data frame of metadata
#' @param join_by which column to join by
#' @export
tor_add_metadata <- function(data, metadata, join_by, quiet = FALSE) {

  if (missing(data)) stop("no data frame supplied", call. = FALSE)
  if (missing(metadata)) stop("no metadata supplied", call. = FALSE)
  if (missing(join_by)) stop("no join_by column specified", call. = FALSE)

  if (!quiet)
    glue("Info: adding metadata to mass spec data, joining by '{join_by}'...") %>%
      message(appendLF = FALSE)

  # combine
  data_combined <-
    left_join(data, mutate(metadata, ..rid.. = row_number()), by = join_by) %>%
    mutate(missing_metadata = is.na(..rid..)) %>%
    select(-..rid..)

  # info
  if (!quiet)
    glue("{nrow(metadata)} metadata entries successfully added to ",
         "{nrow(filter(data_combined, !missing_metadata))} data recors, ",
         "{nrow(filter(data_combined, missing_metadata))} could not be matched to metadata") %>%
    message()

  return(data_combined)
}



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
