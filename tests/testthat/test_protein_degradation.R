context("Testing protein degradation calculation functions")

test_that("Testing file quality", {

  #test that input data file is correct format
  expect_error(calculate_label_rate())
  expect_error(calculate_label_rate(5), "wrong data type")

})

test_that("Testing protein labeling rate", {

  test_data <- data_frame(
    protein = c("6PGD", "7ABC", "6PGD"),
    isopep = c("VLSGPQAQPAGDK", "LLSGPRD", "YAGHMPQFHSLY"),
    frac_lab = c(.5, .6, .5),
    hours = c(1, 0, 0))

  expect_message(calculate_label_rate(test_data), "proteins")
  expect_message(calculate_label_rate(test_data, combine_peptides = FALSE), "peptides")


})
