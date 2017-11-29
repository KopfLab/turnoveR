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
#' @description bla
#' @param svm_data the SVM data
#' @param pred_cutoff the probability cutoff from the svm prediction
process_svm_data <- function(svm_data, pred_cutoff = 0.75, quiet = FALSE) {

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

