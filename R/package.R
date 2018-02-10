#' @keywords internal
"_PACKAGE"

#' @importFrom rlang !! !!! sym syms
#' @importFrom readr read_csv cols_only col_character col_double
#' @importFrom readxl read_excel
#' @importFrom tibble enframe
#' @importFrom dplyr as_data_frame tbl_df %>% filter mutate data_frame left_join select rename bind_rows arrange group_by ungroup do everything summarize
#' @importFrom glue glue collapse
#' @importFrom stringr str_c str_detect str_replace str_to_lower
#' @importFrom purrr safely map map2 map_lgl map_dbl map_int map2_int map_chr map2_chr
#' @importFrom tidyr nest unnest spread gather
#' @importFrom broom tidy glance
#' @importFrom KEGGREST keggGet keggConv keggLink keggList
NULL
