context("Testing SVM files")

test_that("Testing read functions", {

  expect_error(read_svm_data_file("DNE"), "svm file does not seem to exist")


})
