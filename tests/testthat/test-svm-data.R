context("Testing SVM files")

test_that("Testing read functions", {

  expect_error(read_svm_data_file("DNE"), "svm file does not seem to exist")

  # read a good test file
  expect_true(is.data.frame(svm_data <- read_svm_data_file(file.path("test_data", "svm_pred_results.csv"))))
  expect_equal(nrow(svm_data), 101721)

  # check for warnings in a bad test file
  expect_warning(read_svm_data_file(file.path("test_data", "svm_pred_results_bad.csv")), "don't match the column names")
})


test_that("SVM processiong works", {

  expect_error(process_svm_data())
  expect_error(process_svm_data(5), "wrong data type")

})

test_that("Renaming proteins works", {

  # make sure errors are thrown
  expect_error(rename_proteins())
  expect_error(rename_proteins(5), "wrong data type")
  expect_error(rename_proteins(data_frame(x=1)), "column .* does not exist")
  expect_error(rename_proteins(data_frame(prot=1), prot_col = "my_protein"), "column .* does not exist")
  expect_error(rename_proteins(data_frame(prot=1)), "mapping text file.* is required")
  expect_error(rename_proteins(data_frame(prot=1), renaming_protein_map_file = "DNE"))
  expect_error(rename_proteins(data_frame(prot = 1), file.path("test_data", "rename_prot_bad.xlsx")), "column .* does not exist")

  # test specific data
  my_test_data <- data_frame(prot = c("Q1PI83:KAD", "Q1PI83:KAD", "OTHER"))
  expect_message(result <- rename_proteins(my_test_data, file.path("test_data", "rename_prot.xlsx")), "2 protein .* 1 different")
  expect_equal(result, data_frame(prot = c("KAD", "KAD", "OTHER")))

})

#tests for "calculate_fraculab" function
test_that("calculate_fraculab works", {

# - test that filtered_data file is correct format
  expect_error(calculate_fraculab())
  expect_error(calculate_fraculab(5), "wrong data type")

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



#tests for "calc_pep_degrate" function
# - test that input file is correct format
# - test that exponential curve fits data(?)
# - test that correct value is saved for deg rate & chisq

#tests for "make_prot_master" function
# - test that input file is correct format
# - test that values are summed per proper column (protein)
# - test that filtering of data based on number peptides was succesful
# - test that filtering of data based on chisq value was succesful

#tests for "calc_prot_dissipation" function
# - test that input file is correct format
# - test that kdegReal is calculated with correct inputs (kdegavg and growth rate)
# - test that dissipation is calculated with correct inputs (kdegReal and growth rate)


