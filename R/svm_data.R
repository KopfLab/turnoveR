#' Read SVM data file
#' @description This functions reads the data file
#' @param filepath the path to the svm data file
#' @return data_frame with the read data
#' @export
read_svm_data_file <- function(filepath) {

  if(!file.exists(filepath)) {
    glue("This svm file does not seem to exist '{filepath}'") %>%
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




#' Process the SVM data
#' @description process and filter data (Step 10 in current workflow)
#' @param svm_data the SVM data
#' @param pred_cutoff the probability cutoff from the svm prediction
process_svm_data <- function(svm_data, pred_cutoff = 0.75, quiet = FALSE) {

#add steps to rename proteins to match master list? (r code exists for this - asked Matt)
  #a. read in "rename_prot.txt"
  #b. compare protein names in SVM data to list
  #c. change name depending on...?


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
    n_discarded <- nrow(svm_data) - nrow(filtered_data)
    glue("Info: discarded {n_discarded} of {n_original} ({round(n_discarded/n_original*100, 1)}%) rows") %>%
      message()
  }

  return(filtered_data)
}


#' @param data the data set
#' @param renaming_protein_map_file the filepath to the xlsx mapping file
#' @param prot_col the name of the column that has the protein ID/name
rename_proteins <- function(data, renaming_protein_map_file, prot_col = "prot", prot_new_col = "protNew", quiet = FALSE) {

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


#' Calculate labeled fraction 1
#' @description Calculate fraclab, add to metadata(replace excel step 11)
#' @param  filtered_data the filtered svm data
#' @param meta_data file?? should we write a separate function to read in metadata file?
#' Maybe ask Matt what format is most useful (what else is in metadata that they use?)
calculate_fraculab <- function(data, meta_data) {
  data <- data %>%
    mutate(
      #create new columns frac_ulab and frac_lab using ampU and ampL values from filtered dataset
      frac_ulab = ampU / (ampU + ampL),
      frac_lab = 1 - frac_ulab) %>%
    left_join(meta_data, by = "sample") %>%
    select(
      unishort = uniShort,
      protein = prot,
      gene,
      isopep = seqz,
      sample,
      frac_lab
    ) %>%
    arrange(sample, protein) %>%
    spread(sample, frac_lab, fill = 0)

  return(data)
}


#OR


#' Calculate labeled fraction and clean: added step to filter based on number timepoints where the peptide is identified
#' @description Calculate fraclab, calculate number of occurences, add to metadata (Step 11 in workflow)
#' @param  filtered_data the filtered data
#' @param meta_data file?? should we write a separate function to read in metadata?
#' @param min number of timepoints present
#' another option is to make a separate function for filtering based on number of time points present
calculate_fraculab_clean <- function(filtered_data, meta_data, min_timepoint_present = 3) {
#meta_data <- readxl::read_excel(file.path("test_scripts", "metadata.xlsx"))

  fraculab_data_clean <- filtered_data %>%
  rename(unishort = uniShort, protein = prot, isopep = seqz) %>%
  mutate(
    frac_ulab = ampU / (ampU + ampL),
    frac_lab = 1-frac_ulab
  ) %>%
  left_join(meta_data, by = "sample") %>%
  group_by(protein, isopep) %>%
  # calculate the number of time points where the peptide is identified
  mutate(
    n_time_points_identified = n()
  ) %>% ungroup() %>%
  # filter out the ones that are not present in enough time points
  filter(n_time_points_identified > min_timepoint_present)

  return(fraculab_data_clean)
}




#' Plot degradation curves
#' @description create graphs showing proportion of unlabeled fraction over time to visualize degradation
#' @param fraculab_data or fraculab_data_clean
#' @param number of proteins (or peptides?) to plot?
plot_deg_curve <- function(fraculab_data_clean, prots_to_plot = 10) {

  fraculab_data_clean %>%
    filter(prot %in% unique(prot)[1:prots_to_plot]) %>%
    ggplot() +
    aes(hours - min(hours), frac_lab, color = prot, shape = replicate, size = svmPred) +
    geom_point() +
    theme(legend.position="none") +
    facet_wrap(~prot)
  #library(plotly)
  #ggplotly()

}




#' Calculate degradation rate, chisquared value
#' @description calculate kdeg and chisquared (Step 12a-c in workflow, previously performed in Igor)
#' @param fraculab_data or fraculab_data_clean
calc_pep_degrate <- function(fraculab_data_clean) {

  # fit each peptide to exponential equation y = Ae^dx (timepoint, frac_lab)
  #extract d value for each peptide, store in new column
  #calculate chisquared value, store in new column

  # output example: return(deg_rate_data)
}




#' make master list: combine peptide data to protein level
#' @description Average kdeg values of peptides per protein (Step 12d-e in workflow, previously performed in Igor)
#' @param deg_rate_data
#' @param min_chisq
#' @param min_num_peptides
make_prot_master <- function(deg_rate_data, min_chisq = 3, min_num_peptides=2) {

  # average kdeg values for all peptides that correspond to single protein = kdegavg
  # calculate number of peptides that correspond to sing protein = npep
  # remove proteins that not meet min_chisq and min_num_peptide criteria

  # output example: return(prot_sum_data)
}




#' Calculate real degradation rate, calculate dissipation rate
#' @description Calculate the percent of the protein that was degraded during one generation(Step 12f-g in workflow, previously performed in Igor)
#' @param prot_sum_data
#' where does growth rate value come from?!
calc_prot_dissipation <- function(prot_sum_data) {

  # kdegReal = kdegavg - growth rate
  # dissipation = kdegReal / (kdegReal + growth rate)

  # output example: return(prot_disp_data)
}
