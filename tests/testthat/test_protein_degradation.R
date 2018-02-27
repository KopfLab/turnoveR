context("Testing protein degradation calculation functions")

test_that("Testing file quality", {

  #mockdataframe
  mock_data <- data_frame(
    protein = c("6PGD", "7ABC", "6PGD"),
    isopep = c("VLSGPQAQPAGDK", "LLSGPRD", "YAGHMPQFHSLY"),
    frac_lab = c(.5, .6, .5),
    hours = c(1, 0, 3))

  #test that input data file is correct format
  expect_error(tor_calculate_label_rate(), "need .* data set")
  expect_error(tor_calculate_label_rate(5), "wrong data type")

  #test that input data file is complete
  expect_error(tor_calculate_label_rate(mock_data %>% select(-hours)), ".* does not exist") #add appropriate error messages
  expect_error(tor_calculate_label_rate(mock_data %>% select(-protein)), ".* does not exist")
  expect_error(tor_calculate_label_rate(mock_data %>% select(-isopep)), ".* does not exist")
  expect_error(tor_calculate_label_rate(mock_data %>% select(-frac_lab)), ".* does not exist")

  #mockdataframe2
  mock_data2 <- data_frame(
    label_rate = c(0.09220716, 0.07712532, 0.08902904),
    label_rate_se = c(0.0024814317, 0.0036241348, 0.0014125349))

  #test that input data file is correct format
  expect_error(tor_calculate_degradation_dissipation(), "need .* data set")
  expect_error(tor_calculate_degradation_dissipation(5), "wrong data type")
  expect_error(tor_calculate_degradation_dissipation(mock_data2), "need .* growth rate")

  #test that input data file is complete
  expect_error(tor_calculate_degradation_dissipation(mock_data2 %>% select(-label_rate), 0.066), ".* does not exist")
  expect_error(tor_calculate_degradation_dissipation(mock_data2 %>% select(-label_rate_se), 0.066), ".* does not exist")


})

test_that("Testing protein labeling rate", {

  test_data <- data_frame(
    protein = c("6PGD", "7ABC", "6PGD"),
    gene = c("6PGD1", "7ABC2", "6PGD3"),
    isopep = c("VLSGPQAQPAGDK", "LLSGPRD", "YAGHMPQFHSLY"),
    frac_lab = c(.5, .6, .5),
    hours = c(1, 0, 3))

  expect_message(output <- tor_calculate_label_rate(test_data), "proteins")
  expect_message(tor_calculate_label_rate(test_data, combine_peptides = FALSE), "peptides")

})
