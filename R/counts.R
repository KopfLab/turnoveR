#@ note: needs unit testing

#' Read protein counts data
#'
#' Reads the protein counts data from a `psms.csv` file. Note that the columns from the file get slightly renamed for compatibility with other functions.
#'
#' @param psms_file_path the path to the protein counts csv file - expected columns are \code{protein} (the protein ID) and \code{T0_counts} - the protein counts at the T0 condition.
#' @export
tor_read_protein_counts_data <- function(psms_file_path, quiet = FALSE) {

  if (missing(psms_file_path)) stop("no file path to the protein abundance sums provided", call. = FALSE)
  if (!file.exists(psms_file_path))
    glue("file path '{psms_file_path}' does not point to an existing file") %>% stop(call. = FALSE)

  if (!quiet) {
    glue("Info: reading protein counts data from '{basename(psms_file_path)}'... ") %>%
      message(appendLF = FALSE)
  }

  protein_counts <-
    read_csv(
      psms_file_path,
      col_types = cols(
        protein = col_character(),
        T0_counts = col_integer()
      )
    ) %>%
    select(prot_id = protein, prot_counts = T0_counts)

  if (!quiet) {
    glue("read {nrow(protein_counts)} records") %>% message()
  }

  return(protein_counts)
}

#' Add protein counts info
#'
#' Add protein counts to the labeling rate data. This function adds the protein counts and only keeps proteins that have counts > 0 (with a warning). It also calculates relative abundances by counts and mass and thus requires the \link{tor_add_uniprot_info} function to have provided molecular masses of the proteins beforehand.
#'
#' @param data the data with the calculated labeling rates
#' @param protein_counts the protein counts data frame
#' @export
tor_add_protein_counts_info <- function(data, protein_counts, quiet = FALSE) {

  if (missing(data)) stop("no data frame provided", call. = FALSE)
  if (missing(protein_counts)) stop("no protein counts provided", call. = FALSE)

  # safety checks
  columns <- c("prot_id")
  missing <- setdiff(columns, names(data))
  if (length(missing) > 0) {
    glue("column '{collapse(missing, sep = ', ', last = ' and ')}' does not exist in the dataset") %>% stop(call. = FALSE)
  }

  columns <- c("prot_id", "prot_counts")
  missing <- setdiff(columns, names(protein_counts))
  if (length(missing) > 0) {
    glue("column '{collapse(missing, sep = ', ', last = ' and ')}' does not exist in the protein counts data frame") %>% stop(call. = FALSE)
  }

  # combine
  data_combined <-
    full_join(
      mutate(protein_counts, ..pcid.. = row_number()),
      mutate(data, ..did.. = row_number()),
      by = "prot_id"
    )

  # user info
  if (!quiet) {
    n_zero_counts <- filter(data_combined, !is.na(..did..), prot_counts == 0) %>% nrow()
    n_no_counts <- filter(data_combined, !is.na(..did..), is.na(..pcid..)) %>% nrow()
    n_added <- filter(data_combined, !is.na(..did..), !is.na(..pcid..)) %>% nrow()
    glue("Info: protein counts added and weighted rates calculated ",
         "for {n_added}/{nrow(data)} datasets") %>% message(appendLF = FALSE)
    if (n_no_counts > 0)
      glue(", {n_no_counts} record(s) have missing counts") %>% message(appendLF = FALSE)
    if (n_zero_counts > 0)
      glue(", {n_zero_counts} record(s) discarded because of zero counts") %>% message()
    else
      message()
  }

  # focus on just what has been matched
  data_out <- data_combined %>%
    filter(!is.na(..did..)) %>%
    select(-..pcid.., -..did..)

  # calculate relative weighted values
  data_out <- data_out %>%
    mutate(
      # calculate relative counts and relative mass
      prot_rel_counts = prot_counts / sum(prot_counts, na.rm = TRUE),
      prot_rel_mass = (prot_counts * prot_mw) / sum(prot_counts * prot_mw, na.rm = TRUE),
      #calculate a weighted degradation rate
      deg_rate_weighted = deg_rate * prot_rel_mass,
      #calculate a weighted dissipation rate
      dissipation_weighted = dissipation * prot_rel_mass
    )
}
