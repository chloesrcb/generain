# Description: This script contains tests for the `dependence.R` functions.

### Test the `get_chiq` function -----------------------------------------------
test_that("get_chiq works well", {

  # Example 1: Basic test with a small dataset
  data1 <- data.frame(x = c(1, 2, 3), y = c(4, 5, 6))
  quantile1 <- 0.5
  result1 <- get_chiq(data1, quantile1)
  expect_type(result1, "double")
  expect_length(result1, 1)

  # Example 2: Test with a different quantile
  quantile2 <- 0.75
  result2 <- get_chiq(data1, quantile2)
  expect_type(result2, "double")
  expect_length(result2, 1)

  # Example 3: Test with another dataset
  data2 <- data.frame(x = c(10, 20, 30), y = c(30, 20, 10))
  quantile3 <- 0.25
  result3 <- get_chiq(data2, quantile3)
  expect_type(result3, "double")
  expect_length(result3, 1)

  # Example 4: Test edge case where quantile is very small (close to 0)
  quantile_small <- 0.01
  result4 <- get_chiq(data1, quantile_small)
  expect_type(result4, "double")
  expect_length(result4, 1)

  # Example 5: Test edge case where quantile is very large (close to 1)
  quantile_large <- 0.99
  result5 <- get_chiq(data1, quantile_large)
  expect_type(result5, "double")
  expect_length(result5, 1)

  # Example 6: Test with a dataset containing negative values
  data3 <- data.frame(x = c(-1, -2, -3), y = c(-4, -5, -6))
  quantile4 <- 0.5
  result6 <- get_chiq(data3, quantile4)
  expect_type(result6, "double")
  expect_length(result6, 1)
})

test_that("get_chiq returns correct chi value", {

  # Example 1: Basic test with a small dataset
  data1 <- data.frame(x = c(1, 2, 3), y = c(4, 5, 6))
  quantile1 <- 0.5
  
  # Compute the expected value manually
  n <- nrow(data1)
  data_unif <- cbind(
    rank(data1[, 1]) / (n + 1),
    rank(data1[, 2]) / (n + 1)
  )
  rowmax <- apply(data_unif, 1, max)
  cu <- mean(rowmax < quantile1)
  chiu <- 2 - log(cu) / log(quantile1)
  chiulb <- 2 - log(pmax(2 * quantile1 - 1, 0)) / log(quantile1)
  expected1 <- pmax(chiu, chiulb)
  
  # Call the function and compare
  result1 <- get_chiq(data1, quantile1)
  expect_equal(result1, expected1, tolerance = 1e-6)
  
  # Example 2: Test with a different quantile
  quantile2 <- 0.75
  
  # Compute the expected value manually
  cu <- mean(rowmax < quantile2)
  chiu <- 2 - log(cu) / log(quantile2)
  chiulb <- 2 - log(pmax(2 * quantile2 - 1, 0)) / log(quantile2)
  expected2 <- pmax(chiu, chiulb)
  
  # Call the function and compare
  result2 <- get_chiq(data1, quantile2)
  expect_equal(result2, expected2, tolerance = 1e-6)

  # Example 3: Test with another dataset
  data2 <- data.frame(x = c(10, 20, 30), y = c(30, 20, 10))
  quantile3 <- 0.25

  # Compute the expected value manually
  data_unif2 <- cbind(
    rank(data2[, 1]) / (nrow(data2) + 1),
    rank(data2[, 2]) / (nrow(data2) + 1)
  )
  rowmax2 <- apply(data_unif2, 1, max)
  cu2 <- mean(rowmax2 < quantile3)
  chiu2 <- 2 - log(cu2) / log(quantile3)
  chiulb2 <- 2 - log(pmax(2 * quantile3 - 1, 0)) / log(quantile3)
  expected3 <- pmax(chiu2, chiulb2)

  # Call the function and compare
  result3 <- get_chiq(data2, quantile3)
  expect_equal(result3, expected3, tolerance = 1e-6)
  
  # Example 4: Test edge case where quantile is very small (close to 0)
  quantile_small <- 0.01

  # Compute the expected value manually
  cu_small <- mean(rowmax < quantile_small)
  chiu_small <- 2 - log(cu_small) / log(quantile_small)
  chiulb_small <- 2 - log(pmax(2 * quantile_small - 1, 0)) / log(quantile_small)
  expected4 <- pmax(chiu_small, chiulb_small)

  # Call the function and compare
  result4 <- get_chiq(data1, quantile_small)
  expect_equal(result4, expected4, tolerance = 1e-6)

  # Example 5: Test edge case where quantile is very large (close to 1)
  quantile_large <- 0.99

  # Compute the expected value manually
  cu_large <- mean(rowmax < quantile_large)
  chiu_large <- 2 - log(cu_large) / log(quantile_large)
  chiulb_large <- 2 - log(pmax(2 * quantile_large - 1, 0)) / log(quantile_large)
  expected5 <- pmax(chiu_large, chiulb_large)

  # Call the function and compare
  result5 <- get_chiq(data1, quantile_large)
  expect_equal(result5, expected5, tolerance = 1e-6)

  # Example 6: Test with a dataset containing negative values
  data3 <- data.frame(x = c(-1, -2, -3), y = c(-4, -5, -6))
  quantile4 <- 0.5

  # Compute the expected value manually
  data_unif3 <- cbind(
    rank(data3[, 1]) / (nrow(data3) + 1),
    rank(data3[, 2]) / (nrow(data3) + 1)
  )
  rowmax3 <- apply(data_unif3, 1, max)
  cu3 <- mean(rowmax3 < quantile4)
  chiu3 <- 2 - log(cu3) / log(quantile4)
  chiulb3 <- 2 - log(pmax(2 * quantile4 - 1, 0)) / log(quantile4)
  expected6 <- pmax(chiu3, chiulb3)

  # Call the function and compare
  result6 <- get_chiq(data3, quantile4)
  expect_equal(result6, expected6, tolerance = 1e-6)
})


### Test the `chispatemp_empirical` function -----------------------------------

# Example simple data
set.seed(123)
data_rain <- data.frame(
  site1 = rexp(100),
  site2 = rexp(100),
  site3 = rexp(100)
)

quantile <- 0.8  # Quantile threshold
tmax <- 5  # Tmax value

quantile <- 0.8
tmax <- 5

test_that("Test temporal_chi without lag", {
  chi_temp <- temporal_chi(data_rain, tmax, quantile, zeros = TRUE, mean = FALSE)

  # Test 1: lag 0 -> 1
  chi_nolag <- chi_temp[, 1]
  expect_true(all(chi_nolag == 1 | is.na(chi_nolag)))

  # Petite aide pour reproduire exactement temporal_chi()
  chi_expected <- function(x, q, lag) {
    r <- rank(x) / (length(x) + 1)
    denom <- sum(r > q)
    if (lag >= length(r) || denom == 0) return(NA_real_)
    R1 <- r[1:(length(r) - lag)]
    R2 <- r[(1 + lag):length(r)]
    sum(R1 > q & R2 > q) / denom
  }

  # Test 2: site1, lag 1
  chival <- chi_expected(data_rain$site1, quantile, 1)
  expect_equal(chi_temp[1, 2], chival, tolerance = 1e-6)

  # Test 3: site1, lag 5
  chival <- chi_expected(data_rain$site1, quantile, 5)
  expect_equal(chi_temp[1, 6], chival, tolerance = 1e-6)

  # Test 4: site3, lag 5 (utiliser bien site3 pour le dénominateur)
  chival <- chi_expected(data_rain$site3, quantile, 5)
  expect_equal(chi_temp[3, 6], chival, tolerance = 1e-6)
})
