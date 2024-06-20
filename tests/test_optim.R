library(testthat)

test_that("get_criterion calculates the correct values", {
  # Create a sample dataframe
  df_result <- data.frame(
    beta1 = c(1, 2, 3),
    beta2 = c(4, 5, 6),
    alpha1 = c(7, 8, 9),
    alpha2 = c(10, 11, 12)
  )

  # Define true parameters
  true_param <- c(2, 5, 8, 11)

  # Call the function
  result <- get_criterion(df_result, true_param)

  # Expected results
  expected_mean <- c(mean(df_result$beta1), mean(df_result$beta2),
                     mean(df_result$alpha1), mean(df_result$alpha2))
  expected_rmse <- c(sqrt(mean((true_param[1] - df_result$beta1)^2)),
                    sqrt(mean((true_param[2] - df_result$beta2)^2)),
                     sqrt(mean((true_param[3] - df_result$alpha1)^2)),
                     sqrt(mean((true_param[4] - df_result$alpha2)^2)))
  expected_mae <- c(mean(abs(true_param[1] - df_result$beta1)),
                    mean(abs(true_param[2] - df_result$beta2)),
                    mean(abs(true_param[3] - df_result$alpha1)),
                    mean(abs(true_param[4] - df_result$alpha2)))

  # Compare results
  expect_equal(result$mean, expected_mean)
  expect_equal(result$rmse, expected_rmse)
  expect_equal(result$mae, expected_mae)
})


test_that("neg_ll calculates the correct negative log-likelihood", {
  # Create sample data and parameters
  params <- c(0.5, 0.7, 1.2, 1.5) # Example parameters
  excesses <- list(N_vect = c(10, 20, 30), n_vect = c(2, 4, 6)) # Excesses
  h_vect <- c(1, 2, 3) # Example h_vect
  hmax <- 3 # Example hmax
  tau <- 1 # Example tau
  df_dist <- data.frame(value = c(1, 2, 3)) # Example df_dist (distance matrix
  locations <- data.frame(Latitude = c(1, 2, 7),
                         Longitude = c(3, 4, 7)) # Example locations

  # Call the function
  result <- neg_ll(params, excesses, h_vect, hmax, tau, df_dist, locations)

  # Check if the result is a numeric value
  expect_is(result, "numeric")

  # Check if the result is a Inf for wrong parameters
  params <- c(0.5, 0.7, 1.2, 3)
  result <- neg_ll(params, excesses, h_vect, hmax, tau, df_dist, locations)
  expect_equal(result, 1e+06)

  params <- c(0.5, 0.7, 2, 1)
  result <- neg_ll(params, excesses, h_vect, hmax, tau, df_dist, locations)
  expect_equal(result, 1e+06)

  params <- c(0, 0.7, 1, 1)
  result <- neg_ll(params, excesses, h_vect, hmax, tau, df_dist, locations)
  expect_equal(result, 1e+06)

  params <- c(0.5, 0, 1, 1)
  result <- neg_ll(params, excesses, h_vect, hmax, tau, df_dist, locations)
  expect_equal(result, 1e+06)

  # with advection=c(0, 0)
  params <- c(0.5, 0.5, 1, 1, 0, 0)

  # # with advection
  # df_dist <- data.frame(value = c(1, 2, 3, 4, 5, 6)) # Example df_dist (distance matrix
  params <- c(0.5, 0.5, 1, 1, 0.1, 0.1)
  # result <- neg_ll(params, excesses, h_vect, hmax, tau, df_dist, locations)
  # expect_equal(result, 1e+06)
})

test_that("empirical_excesses works with sufficient data", {
  # Create sample data
  data_rain <- matrix(runif(100 * 25, 0, 100), nrow = 100, ncol = 25)
  data_rain <- as.data.frame(data_rain)
  sites_coords <- generate_grid_coords(5) # Example sites_coords
  dist_mat <- get_dist_mat(sites_coords, adv = c(0, 0), tau = 1:10,
                            latlon = FALSE) # distance matrix
  df_dist <- reshape_distances(dist_mat) # reshape the distance matrix

  # Calculate the h vector
  h_vect <- get_h_vect(df_dist, sqrt(17))
  quantile <- 0.9
  tau <- 1:10

  # Call the function
  result <- empirical_excesses(data_rain, quantile, tau, h_vect, df_dist,
                               nmin = 5)

  # Check if the result as the correct length
  expect_true(length(result$n_vect) == nrow(df_dist))
  expect_true(length(result$N_vect) == nrow(df_dist))

  # Check if the result is numeric or NA
  expect_true(all(is.na(result$n_vect) | is.numeric(result$n_vect)))
  expect_true(all(is.na(result$N_vect) | is.numeric(result$N_vect)))

  # Check if the result is NA for insufficient data
  result <- empirical_excesses(data_rain, quantile, tau, h_vect, df_dist,
                               nmin = 1000)
  expect_true(all(is.na(result$n_vect)))
  expect_true(all(is.na(result$N_vect)))

  # with advection coordinates
  dist_mat <- get_dist_mat(sites_coords, adv = c(0.1, 0.1), tau = 1:10,
                            latlon = FALSE) # distance matrix
  df_dist <- reshape_distances(dist_mat) # reshape the distance matrix

  # Calculate the h vector
  h_vect <- get_h_vect(df_dist, sqrt(17))
  quantile <- 0.9
  tau <- 1:10

  # Call the function
  result <- empirical_excesses(data_rain, quantile, tau, h_vect, df_dist,
                               nmin = 5)

  # Check if the result as the correct length
  expect_true(length(result$n_vect) == nrow(df_dist))
  expect_true(length(result$N_vect) == nrow(df_dist))

  # Check if the result is numeric or NA
  expect_true(all(is.na(result$n_vect) | is.numeric(result$n_vect)))
  expect_true(all(is.na(result$N_vect) | is.numeric(result$N_vect)))
})
