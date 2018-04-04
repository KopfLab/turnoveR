# finished --

#' Read SVM data file
#'
#' This functions reads SVM preprocessed data files and selects the most important information
#'
#' @param filepath the path to the svm data file
#' @return data_frame with the read data
#' @export
tor_read_svm_data_file <- function(filepath, quiet = FALSE) {

  if(!file.exists(filepath)) {
    glue("This svm file does not seem to exist '{filepath}' in the current working directory") %>%
    stop(call. = FALSE)
  }

  data <-
    read_csv(filepath, col_types =
               # we're reading all general columns that are important for downstream data processing, skipping some that are not
               cols_only(
                 uniShort = col_character(),
                 prot = col_character(),
                 gene = col_character(),
                 seqz = col_character(),
                 mz = col_double(),
                 ppmerr = col_double(),
                 svmPred = col_double(),
                 ampU = col_double(),
                 ampL = col_double(),
                 sample = col_character()
               )
    ) %>%
    as_data_frame()

  # information for user
  if (!quiet) {
    glue("Info: successfully read {nrow(data)} records from SVM file '{basename(filepath)}'") %>%
      message()
  }

  # check for missing data
  columns <- c(
    prot_id = "prot", uniprot_id = "uniShort", peptide_seq = "seqz",
    peptide_mz = "mz", svm_pred = "svmPred", amp_ulab = "ampU",
    amp_lab = "ampL", sample = "sample")
  missing <- setdiff(unname(columns), names(data))
  if (length(missing) > 0) {
    glue("missing column(s) '{collapse(missing, sep = ', ', last = ' and ')}' do(es) not exist in the SVM file") %>% stop(call. = FALSE)
  }

  # select & rename the key field
  data %>% select(!!!map(columns, sym))
}


# finished but will have to work with filtering criteria other than pred_cutoff (once we determine which filter criteria are actually most useful)

#' Filter peptides
#'
#' Process and filter data (Step 10 in current workflow) by spectral quality evaluation.
#'
#' @param data the data set
#' @param condition the filtering condition
#' @export
tor_filter_peptides_by_spectral_fit_quality <- function(svm_data, condition, quiet = FALSE) {

  if (missing(svm_data)) stop("need to supply the svm data frame", call. =FALSE)
  if (!is.data.frame(svm_data)) {
    glue("wrong data type supplied: {class(svm_data)[1]}") %>% stop(call. = FALSE)
  }
  if (missing(condition)) stop("no filtering condition supplied", call. = FALSE)

  condition_quo <- enquo(condition)
  filtered_data <- svm_data %>%
    # filter out peptides that do have a low probability score
    filter(!!condition_quo)


  # information for user
  if (!quiet) {
    n_original <- nrow(svm_data)
    n_kept <- nrow(filtered_data)
    glue("Info: kept {n_kept} of {n_original} ({round(n_kept/n_original*100, 1)}%) peptide measurements during spectral fit quality filtering (condition '{quo_text(condition_quo)}')") %>%
      message()
  }

  return(filtered_data)
}


# -- finished

#' Rename proteins based on protein mapping file
#'
#' This function may be obsolete once uniprot integration is completed.
#'
#' @param data the data set
#' @param renaming_protein_map_file the filepath to the xlsx mapping file
#' @param prot_col the name of the column that has the protein ID/name
#' @export
tor_recode_protein_ids <- function(data, renaming_protein_map_file, prot_col = "prot_id", prot_new_col = "new_prot_id", quiet = FALSE) {

  # safety checks for data
  if (missing(data)) stop("need to supply a data set", call. =FALSE)
  if (!is.data.frame(data)) {
    glue("wrong data type supplied: {class(data)[1]}") %>% stop(call. = FALSE)
  }
  if (!prot_col %in% names(data)){
    glue("protein ID column '{prot_col}' does not exist in the dataset") %>% stop(call. = FALSE)
  }

  # safety checks for renaming
  if (missing(renaming_protein_map_file)) {
    stop("parameter 'renaming_protein_map_file' (the path to the mapping text file) is required", call. =FALSE)
  }
  if (!file.exists(renaming_protein_map_file))  {
    glue("Mapping text file '{renaming_protein_map_file}' does not exist. Please check that the path is correct.") %>% stop(call. = FALSE)
  }

  mapping_file <- read_excel(renaming_protein_map_file)
  if (!prot_col %in% names(mapping_file)){
    glue("protein ID column '{prot_col}' does not exist in the mapping file") %>% stop(call. = FALSE)
  }
  if (!prot_new_col %in% names(mapping_file)){
    glue("protein new ID column '{prot_new_col}' does not exist in the mapping file") %>% stop(call. = FALSE)
  }

  # do the mapping
  combined_data <-
    left_join(data, mapping_file, by = prot_col)

  # information for user
  if (!quiet) {
    #how_many_renamed <- sum(!is.na(combined_data[[prot_new_col]]))
    how_many_renamed <- combined_data[[prot_new_col]] %>% { !is.na(.) } %>% sum()
    how_many_unique <- combined_data[[prot_new_col]] %>% unique() %>% { !is.na(.) } %>% sum()
    glue("Info: renamed {how_many_renamed} protein entries for {how_many_unique} different proteins.") %>%
      message()
  }

  # because we have variable names for the columns, need to use quoted expressions (_qs)
  prot_col_q <- as.name(prot_col)
  prot_new_col_q <- as.name(prot_new_col)
  combined_data %>%
    mutate(!!prot_col_q := ifelse(!is.na(!!prot_new_col_q), !!prot_new_col_q, !!prot_col_q)) %>%
    select(-!!prot_new_col_q) %>%
    return()

}



# add option for columns being differently names (ampU, amplL)
#' Calculate labeled fraction
#'
#' Cacluated the labeled and unlabaled fraction based on labeld and unlabeld amplitudes.
#'
#' @param data the dataset
#' @param unlabeled_col the unlabeled singal column
#' @param labeled_col the labeled signal column
#' @export
tor_calculate_labeled_fraction <- function(data, unlabeled_col = "amp_ulab", labeled_col = "amp_lab", quiet=FALSE) {

  #checking whether data file was supplied, and in correct format
  if (missing(data)) stop("need to supply the data frame", call. =FALSE)
  if (!is.data.frame(data)) {
    glue("wrong data type supplied: {class(data)[1]}") %>% stop(call. = FALSE)
  }

  if (!unlabeled_col %in% names(data))
    glue("column '{unlabeled_col}' does not exist in this dataset") %>% stop(call. = FALSE)
  if (!labeled_col %in% names(data))
    glue("column '{labeled_col}' does not exist in this dataset") %>% stop(call. = FALSE)

  # adding columns for unlabeled and labeled fraction.
  data <- data %>%
    mutate(
      #create new columns frac_ulab and frac_lab using ampU and ampL values from filtered dataset
      frac_ulab = !!sym(unlabeled_col) / (!!sym(unlabeled_col) + !!sym(labeled_col)),
      frac_lab = 1 - frac_ulab)

  if (!quiet)
    glue("Info: calculated labeled/unlabeled fraction for {nrow(data)} peptides") %>%
    message()

  return(data)
}


