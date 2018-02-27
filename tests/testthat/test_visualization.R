context("Testing visualization functions")

test_that("Testing file quality", {

  #mockdataframe
  mock_data <- data_frame(
    growth_rate = c(0.066, 0.066, 0.066),
    label_rate = c(0.09, 0.08, 0.09)
  )

  #test that input data file is correct format
  expect_error(tor_plot_label_rate_hist(), "need .* data set")
  expect_error(tor_plot_label_rate_hist(5), "wrong data type")

  #test that input data file is complete
  expect_error(tor_plot_label_rate_hist(mock_data %>% select(-growth_rate)), ".* does not exist") #add appropriate error messages
  expect_error(tor_plot_label_rate_hist(mock_data %>% select(-label_rate)), ".* does not exist")

  mock_data2 <- data_frame(
    label_rate = c(0.09220716), #, 0.07712532, 0.08902904
    label_rate_se = c(0.0024814317), # , 0.0036241348, 0.0014125349
    growth_rate = c(0.066), # 0.066, 0.066
    fit_rse = c(5.2e-02, 2.5e-02, 1.1e-02)
  )

  #test that input data file is correct format
  expect_error(tor_plot_label_rate_error(), "need .* data set")
  expect_error(tor_plot_label_rate_error(5), "wrong data type")

  #test that input data file is complete
  expect_error(tor_plot_label_rate_error(mock_data2 %>% select(-label_rate)), ".* does not exist")
  expect_error(tor_plot_label_rate_error(mock_data2 %>% select(-label_rate_se)), ".* does not exist")
  expect_error(tor_plot_label_rate_error(mock_data2 %>% select(-growth_rate)), ".* does not exist")
  expect_error(tor_plot_label_rate_error(mock_data2 %>% select(-fit_rse)), ".* does not exist")


  mock_data3 <- data_frame(
    label_rate = c(0.09220716, 0.07712532, 0.08902904),
    nested_data = c(1, 1, 1),##
    growth_rate = c(0.066, 0.066, 0.066),
    fit = c(1, 1, 1) ##
  )
  ##nested data??

  #test that input data file is correct format
  expect_error(tor_plot_labeling_curves(), "need .* data set")
  expect_error(tor_plot_labeling_curves(5), "wrong data type")

  #test that input data file is complete
  expect_error(tor_plot_labeling_curves(mock_data3 %>% select(-label_rate)), ".* does not exist")
  expect_error(tor_plot_labeling_curves(mock_data3 %>% select(-nested_data)), ".* does not exist")
  expect_error(tor_plot_labeling_curves(mock_data3 %>% select(-growth_rate)), ".* does not exist")
  expect_error(tor_plot_labeling_curves(mock_data3 %>% select(-fit)), ".* does not exist")
})
