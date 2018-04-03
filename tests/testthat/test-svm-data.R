context("Testing SVM files")

test_that("Testing read functions", {

  expect_error(tor_read_svm_data_file("DNE"), "svm file does not seem to exist")

  # read a good test file
  expect_true(is.data.frame(svm_data <- tor_read_svm_data_file(file.path("test_data", "svm_pred_results.csv"))))
  expect_equal(nrow(svm_data), 101721)

  # check for warnings in a bad test file
  expect_error(
    suppressWarnings(tor_read_svm_data_file(file.path("test_data", "svm_pred_results_bad.csv"))),
    "missing column")
})


test_that("SVM processiong works", {

  expect_error(tor_filter_peptides_by_spectral_fit_quality())
  expect_error(tor_filter_peptides_by_spectral_fit_quality(5), "wrong data type")

})

test_that("Renaming proteins works", {

  # make sure errors are thrown
  expect_error(tor_recode_protein_ids())
  expect_error(tor_recode_protein_ids(5), "wrong data type")
  expect_error(tor_recode_protein_ids(data_frame(x=1)), "column .* does not exist")
  expect_error(tor_recode_protein_ids(data_frame(prot_id=1), prot_col = "my_protein"), "column .* does not exist")
  expect_error(tor_recode_protein_ids(data_frame(prot_id=1)), "mapping text file.* is required")
  expect_error(tor_recode_protein_ids(data_frame(prot_id=1), renaming_protein_map_file = "DNE"))
  expect_error(tor_recode_protein_ids(data_frame(prot_id = 1), file.path("test_data", "rename_prot_bad.xlsx")), "column .* does not exist")

  # test specific data
  my_test_data <- data_frame(prot_id = c("Q1PI83:KAD", "Q1PI83:KAD", "OTHER"))
  expect_message(result <- tor_recode_protein_ids(my_test_data, file.path("test_data", "rename_prot.xlsx")), "2 protein .* 1 different")
  expect_equal(result, data_frame(prot_id = c("KAD", "KAD", "OTHER")))
})

#tests for "tor_calculate_labeled_fraction" function
test_that("tor_calculate_labeled_fraction works", {

# - test that filtered_data file is correct format
  expect_error(tor_calculate_labeled_fraction())
  expect_error(tor_calculate_labeled_fraction(5), "wrong data type")

# test specific data
  # test that correct inputs used to calculate
  # test that fraculab and fraclab calculations are correct
})

#tests for "filter_min_timepoints" function
test_that("filter_min_timepoints works", {

  # - test that data file is correct format
  expect_error(filter_min_timepoints())
  expect_error(filter_min_timepoints(5), "wrong data type")

})






