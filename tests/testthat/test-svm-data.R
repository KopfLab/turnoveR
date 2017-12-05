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

#tests for "calculate_fraculab" function
# - test that filtered_data file is correct format
# - test that metadata file is correct format
# - test that fraculab and fraclab calculations are correct (and / or used correct inputs to calculate)

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


