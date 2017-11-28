context("Testing SVM files")

test_that("Testing read functions", {

  expect_error(read_svm_data_file("DNE"), "svm file does not seem to exist")

  # read a good test file
  expect_true(is.data.frame(svm_data <- read_svm_data_file(file.path("test_data", "svm_pred_results.csv"))))
  expect_equal(nrow(svm_data), 101721)

  # check for warnings in a bad test file
  expect_warning(read_svm_data_file(file.path("test_data", "svm_pred_results_bad.csv")), "don't match the column names")
})
