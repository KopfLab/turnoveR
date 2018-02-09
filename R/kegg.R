#' Fetch KEGG pathway details
#'
#' @inheritParams fetch_kegg_details
fetch_kegg_pathway_details <- function(kegg_db_entry) {

  # safety checks
  if (missing(kegg_db_entry)) stop("no pathway kegg db entries provided", call. = FALSE)
  if (!all(str_detect(kegg_db_entry, "^path:")))
    stop("all pathway kegg db IDs must start with 'path:'", call. = FALSE)

  info <- fetch_kegg_details(kegg_db_entry)

  # discard unneded (or confusing) columns if they exist
  discard_cols <- c("organism")
  discard_cols <- intersect(discard_cols, names(pathway_details))
  if ( length(discard_cols) > 0)  {
    info[discard_cols] <- NULL
  }

  return(info)
}


#' Fetch KEGG details
#'
#' Helper function to convert output of \link[KEGGREST]{keggGet} to a nested data frame that is easy to work with in the tidyverse.
#' To unnest individual nested columns and preserve others use the \code{.drop=FALSE} parameter in \link[tidyr]{unnest}
#'
#' @param kegg_db_entry one or more (up to a maximum of 10) KEGG identifiers (see \link[KEGGREST]{keggGet})
#' @param unnest_single_values whether to unnest single values (values that have a none or only a single entry for all retrieve records)
fetch_kegg_details <- function(kegg_db_entry, unnest_single_values = TRUE, ...) {
  # safety checks
  if (missing(kegg_db_entry)) stop("no kegg db entries provided", call. = FALSE)

  # user info
  glue("Info: querying KEGG db for {length(kegg_db_entry)} entries ",
       "('{collapse(kegg_db_entry, sep=\"', '\")}')...") %>%
    message()

  # query kegg db
  tryCatch(
    info <- keggGet(kegg_db_entry),
    error = function(e) {
      stop("KEGG API could not process the request: ", e$message, call. = FALSE)
    })

  # convert info into data frame format
  info_df <-
    data_frame(kegg_info = info) %>%
    mutate(
      nr = row_number(),
      name = map(kegg_info, names)
    ) %>% unnest(name, .drop = FALSE) %>%
    mutate(
      class = map2_chr(kegg_info, name, ~class(.x[[.y]])[1]),
      length = map2_int(kegg_info, name, ~length(.x[[.y]])),
      value = map2(kegg_info, name, ~.x[[.y]]),
      name = str_to_lower(name)
    )

  # info classes
  info_classes <-
    info_df %>%
    group_by(name) %>%
    summarize(
      info_class = unique(class)[1],
      value_max_n = as.integer(max(length)))

  # wide info
  info_df_wide <- info_df %>%
    select(nr, name, value) %>%
    spread(name, value) %>%
    select(-nr)

  # fill NULL values with NA to not loose records during unnesting (except for lists)
  for (i in 1:nrow(info_classes)) {
    info_df_wide <-
      with(info_classes[i,], {
        # make sure the function exists
        if (exists(info_class)) {
          if (info_class %in% c("character", "integer", "numeric"))
            default_value <- do.call(str_c("as.", info_class), args = list(NA))
          else
            default_value <- do.call(info_class, args=list())
          mutate(info_df_wide,
                 !!name := map(!!sym(name), ~if (is.null(.x)) { default_value } else { .x }))
        } else {
          # don't do anything if it's not a standard class
          info_df_wide
        }
      })
  }

  # unnest all the ones that have only single value
  if (unnest_single_values) {
    unnest_cols <- info_classes %>%
      filter(value_max_n == 1, info_class %in% c("character", "integer", "numeric")) %>%
      {.$name}
    info_df_wide <- unnest(info_df_wide, !!!syms(unnest_cols), .drop = FALSE) %>%
      select(!!!syms(unnest_cols), everything())
  }


  # return data frame
  return(info_df_wide)
}
