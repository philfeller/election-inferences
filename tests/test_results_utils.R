# tests/test_results_utils.R
library(testthat)

# Assuming the functions are sourced from results_utils.R
source('./results_utils.R')  # Update with the correct path if necessary

# Test for combine_towns function
test_that("combine_towns combines towns correctly", {
  # Example data for testing
  towns_data <- data.frame(town = c("Town A", "Town B"), value = c(10, 20))
  result <- combine_towns(towns_data)
  
  # Expected result data
  expected_result <- data.frame(town = "Combined", value = 30)
  
  expect_equal(result, expected_result)
})

# Test for dissolve_towns function
test_that("dissolve_towns dissolves towns correctly", {
  # Example data for testing
  towns_data <- data.frame(town = c("Town A", "Town B", "Town C"), value = c(10, 5, 15))
  result <- dissolve_towns(towns_data)
  
  # Expected result after dissolving
  expected_result <- data.frame(town = "Dissolved", value = 30)
  
  expect_equal(result, expected_result)
})

# Test for year_groupings function
test_that("year_groupings groups years correctly", {
  # Example data for testing
  years_data <- data.frame(year = c(2000, 2001, 2002), value = c(5, 10, 15))
  result <- year_groupings(years_data)
  
  # Expected grouped result
  expected_result <- data.frame(group = "Group 2000-2002", total = 30)
  
  expect_equal(result, expected_result)
})

