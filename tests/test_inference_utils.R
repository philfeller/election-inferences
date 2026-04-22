# Unit Tests for Inference Utilities in R

library(testthat)

# Load the functions to be tested
devtools::load_all(".")  # Adjust path as necessary

# Example Data
votes <- data.frame(
  candidate = c("A", "A", "B", "B"),
  votes = c(100, 150, 200, 250)
)

spatial_data <- data.frame(
  id = 1:4,
  geometry = c('POINT(1 1)', 'POINT(2 2)', 'POINT(3 3)', 'POINT(4 4)')
)

# Test create_weighted_avg function

test_that("create_weighted_avg calculates correctly", {
  result <- create_weighted_avg(votes)
  expect_equal(result$average, (100 + 150 + 200 + 250) / 4)  # Adjust based on actual expected behavior
})

# Test prepare_spatial_data function

test_that("prepare_spatial_data formats data correctly", {
  result <- prepare_spatial_data(spatial_data)
  expect_true(is.data.frame(result))
  expect_equal(nrow(result), nrow(spatial_data))
})
