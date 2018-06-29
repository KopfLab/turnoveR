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

# internal functions ======

# function to reliably turn a list of lists into a data frame
# @param lists_df a data frame with a column that has lists in each row
# @param column the name of the column that has the list
# @param unnest_single_values whether to unnest single values (values that have a none or only a single entry for all retrieve records)
unpack_lists_data_frame <- function(lists_df, column = lists, unnest_single_values = TRUE) {

  # convert lists into data frame format
  lists_df <-
    lists_df %>%
    rename(lists = !!enquo(column)) %>%
    mutate(
      nr = row_number(),
      name = map(lists, names)
    ) %>% unnest(name, .drop = FALSE) %>%
    mutate(
      class = map2_chr(lists, name, ~class(.x[[.y]])[1]),
      length = map2_int(lists, name, ~length(.x[[.y]])),
      value = map2(lists, name, ~.x[[.y]]),
      name = str_to_lower(name)
    )

  # data classes
  data_classes <-
    lists_df %>%
    group_by(name) %>%
    summarize(
      data_class = unique(class)[1],
      value_max_n = as.integer(max(length)))

  # lists wide
  lists_df_wide <- lists_df %>%
    select(nr, name, value) %>%
    spread(name, value) %>%
    select(-nr)

  # fill NULL values with NA to not loose records during unnesting (except for lists)
  for (i in 1:nrow(data_classes)) {
    lists_df_wide <-
      with(data_classes[i,], {
        # make sure the function exists
        if (exists(data_class)) {
          if (data_class %in% c("character", "integer", "numeric"))
            default_value <- do.call(str_c("as.", data_class), args = list(NA))
          else
            default_value <- do.call(data_class, args=list())
          # note, could also do this with a right_join back in (but perhaps slower?)
          mutate(lists_df_wide,
                 !!name := map(!!sym(name), ~if (is.null(.x)) { default_value } else { .x }))
        } else {
          # don't do anything if it's not a standard class
          lists_df_wide
        }
      })
  }

  # unnest all the ones that have only single value
  if (unnest_single_values) {
    unnest_cols <- data_classes %>%
      filter(value_max_n == 1, data_class %in% c("character", "integer", "numeric")) %>%
      {.$name}
    lists_df_wide <- unnest(lists_df_wide, !!!syms(unnest_cols), .drop = FALSE) %>%
      select(!!!syms(unnest_cols), everything())
  }

  return(lists_df_wide)
}
