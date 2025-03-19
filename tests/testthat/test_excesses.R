# Test simple pour vérifier si les excès sont bien comptés
test_that("empirical_excesses counts excesses correctly", {
  quantile <- 6
  threshold <- TRUE

  data_rain <- data.frame(
    s1 = c(1, 3, 5, 7, 9, 10, 11),
    s2 = c(2, 4, 6, 6, 10, 2, 10)
  )

  df_lags <- data.frame(
    s1 = c(1, 1, 2),
    s2 = c(1, 2, 2),
    tau = c(0, 0, 0)
  )

  result <- empirical_excesses(data_rain, quantile, df_lags, threshold,
                               type = "brownresnick")

  expect_equal(result$Tobs[1], length(data_rain$s1))
  expect_equal(result$kij[2], 2)
  expect_equal(result$kij[1], 4)

  quantile <- 0.6
  threshold <- FALSE
  u <- quantile(data_rain$s1, quantile)
  nb_excesses_s1 <- sum(data_rain$s1 > u)
  nb_excesses_s2 <- sum(data_rain$s2 > u)
  nb_joint <- sum(data_rain$s1 > u & data_rain$s2 > u)

  result <- empirical_excesses(data_rain, quantile, df_lags, threshold,
                               type = "brownresnick")

  expect_equal(result$Tobs[1], length(data_rain$s1))
  expect_equal(result$kij[2], nb_joint)
  expect_equal(result$kij[1], nb_excesses_s1)
  expect_equal(result$kij[3], nb_excesses_s2)

})


test_that("empirical_excesses_rpar works correctly", {
  # Create a small rainfall dataset (3 time steps, 2 locations)
  data_rain <- matrix(c(5, 10, 15, 20, 25, 30), nrow = 3, byrow = TRUE)
  colnames(data_rain) <- c("s1", "s2")
  # Define quantile threshold
  quantile_value <- 15

  # Define a simple lag dataframe
  df_lags <- data.frame(
    s1 = c(1, 1),
    s2 = c(2, 2),
    tau = c(1, 2),
    kij = NA,
    Tobs = NA
  )

  # Run the function
  result <- empirical_excesses_rpar(data_rain, quantile_value, df_lags,
                                                    threshold = FALSE, t0 = 0)
  
  # Expected values
  expected_Tobs <- c(1, 1)  # Should be 1 based on function logic
  expected_kij <- c(1, 1)   # Since both exceed the quantile threshold
  
  # Assertions
  expect_equal(result$Tobs, expected_Tobs)
  expect_equal(result$kij, expected_kij)
})
