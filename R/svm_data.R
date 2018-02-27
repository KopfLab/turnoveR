# finished --

#' Read SVM data file
#' @description This functions reads the data file
#' @param filepath the path to the svm data file
#' @return data_frame with the read data
#' @export
tor_read_svm_data_file <- function(filepath) {

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
                 sample = col_character(),
                 series = col_character()
               )
    ) %>%
    as_data_frame()

  return(data)
}


# finished but will have to work with filtering criteria other than pred_cutoff (once we determine which filter criteria are actually most useful)

#' Process the SVM data
#' @description process and filter data (Step 10 in current workflow)
#' @param svm_data the SVM data
#' @param pred_cutoff the probability cutoff from the svm prediction
tor_filter_peptides_by_spectral_fit <- function(svm_data, pred_cutoff = 0.75, quiet = FALSE) {

  if (missing(svm_data)) stop("need to supply the svm data frame", call. =FALSE)
  if (!is.data.frame(svm_data)) {
    glue("wrong data type supplied: {class(svm_data)[1]}") %>% stop(call. = FALSE)
  }

  filtered_data <- svm_data %>%
    # filter out peptides that do have a low probability score
    filter(svmPred >= pred_cutoff)


  # information for user
  if (!quiet) {
    n_original <- nrow(svm_data)
    n_kept <- nrow(filtered_data)
    glue("Info: kept {n_kept} of {n_original} ({round(n_kept/n_original*100, 1)}%) peptide measurements") %>%
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
tor_rename_proteins <- function(data, renaming_protein_map_file, prot_col = "prot", prot_new_col = "protNew", quiet = FALSE) {

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
#' @description Calculate fraclab, add to metadata(replace excel step 11)
#' @param  filtered_data the filtered svm data
#' @export
tor_calculate_labeled_fraction <- function(data) {

  #checking whether data file was supplied, and in correct format
  if (missing(data)) stop("need to supply the data frame", call. =FALSE)
  if (!is.data.frame(data)) {
    glue("wrong data type supplied: {class(data)[1]}") %>% stop(call. = FALSE)
  }

  # adding columns for unlabeled and labeled fraction.
  data <- data %>%
    mutate(
      #create new columns frac_ulab and frac_lab using ampU and ampL values from filtered dataset
      frac_ulab = ampU / (ampU + ampL),
      frac_lab = 1 - frac_ulab)


  return(data)
}

# finished -- except for maybe a few more checks AND allow for overwriting default column names
# Q: are we actuall still using this one or is it going to be filtering after the degradation curves are calculated?

#' Filter data based on number of timepoints where the peptide is identified
#' @description (Step 11 in workflow)
#' @param  data the data to be filtered
#' @param min number of timepoints present
filter_min_timepoints <-function(data, min_timepoint_present = 3, quiet = FALSE) {
    #checking whether data file was supplied, and in correct format
    if (missing(data))
      stop("need to supply the data frame", call. = FALSE)
    if (!is.data.frame(data)) {
      glue("wrong data type supplied: {class(data)[1]}") %>% stop(call. = FALSE)
    }

    # @TODO: maybe add safety checks that the expected columns are present: protein, isopep, etc.

    filtered_data <- data %>%
      group_by(protein, isopep) %>%
      # calculate the number of time points where the peptide is identified
      mutate(n_time_points_identified = n()) %>% ungroup() %>%
      # filter out the ones that are not present in enough time points
      filter(n_time_points_identified > min_timepoint_present)

    # information for user
    if (!quiet) {
      n_original <- nrow(data)
      n_kept <- nrow(filtered_data)
      glue(
        "Info: kept {n_kept} of {n_original} ({round(n_kept/n_original*100, 1)}%) rows"
      ) %>%
        message()
    }
    return(filtered_data)
  }



