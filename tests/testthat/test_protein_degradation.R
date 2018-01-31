context("Testing protein degradation calculation functions")

test_that("Testing file quality", {

  #mockdataframe
  mock_data <- data_frame(
    protein = c("6PGD", "7ABC", "6PGD"),
    isopep = c("VLSGPQAQPAGDK", "LLSGPRD", "YAGHMPQFHSLY"),
    frac_lab = c(.5, .6, .5),
    hours = c(1, 0, 3))

  #test that input data file is correct format
  expect_error(calculate_label_rate())
  expect_error(calculate_label_rate(5), "wrong data type")

  #test that input data file is complete
  expect_error(calculate_label_rate(mock_data %>% select(-hours)))
  expect_error(calculate_label_rate(mock_data %>% select(-protein)))
  expect_error(calculate_label_rate(mock_data %>% select(-isopep)))
  expect_error(calculate_label_rate(mock_data %>% select(-frac_lab)))


  #mockdataframe2
  mock_data2 <- data_frame(
    label_rate = c(0.09220716, 0.07712532, 0.08902904),
    label_rate_se = c(0.0024814317, 0.0036241348, 0.0014125349))

  #test that input data file is correct format
  expect_error(calculate_degrate_dissipation())
  expect_error(calculate_degrate_dissipation(5), "wrong data type")

  #test that input data file is complete
  expect_error(calculate_degrate_dissipation(mock_data %>% select(-label_rate)))
  expect_error(calculate_degrate_dissipation(mock_data %>% select(-label_rate_se)))


})

test_that("Testing protein labeling rate", {

  test_data <- data_frame(
    protein = c("6PGD", "7ABC", "6PGD"),
    isopep = c("VLSGPQAQPAGDK", "LLSGPRD", "YAGHMPQFHSLY"),
    frac_lab = c(.5, .6, .5),
    hours = c(1, 0, 3))

  expect_message(output <- calculate_label_rate(test_data), "proteins")
  expect_message(calculate_label_rate(test_data, combine_peptides = FALSE), "peptides")

})
