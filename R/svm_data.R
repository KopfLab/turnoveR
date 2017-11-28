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
