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

test_that("get_lag_vectors", {
  # Create sample data
  data_rain <- matrix(runif(100 * 25, 0, 100), nrow = 100, ncol = 25)
  data_rain <- as.data.frame(data_rain)
  sites_coords <- generate_grid_coords(5) # Example sites_coords

  # Calculate the h vector
  q <- 0.7
  true_param <- c(0.4, 0.2, 1.5, 1)
  hmax <- sqrt(17)
  tau_vect <- 0:10
  df_lags <- get_lag_vectors(sites_coords, true_param,
                          hmax, tau_vect)

  # Check col
  expect_equal(ncol(df_lags), 6)
  # Check taus
  expect_equal(unique(df_lags$tau), tau_vect)
  # Check hmax
  expect_equal(max(df_lags$hnorm), hmax)

})

test_that("empirical_excesses works with sufficient data", {
  # Create sample data
  data_rain <- matrix(runif(100 * 25, 0, 100), nrow = 100, ncol = 25)
  data_rain <- as.data.frame(data_rain)
  sites_coords <- generate_grid_coords(5) # Example sites_coords

  # Calculate the h vector
  q <- 0.7
  true_param <- c(0.4, 0.2, 1.5, 1)
  df_lags <- get_lag_vectors(sites_coords, true_param,
                          hmax = sqrt(17), tau_vect = 0:10)
  excesses <- empirical_excesses(data_rain, quantile = q, df_lags = df_lags)

  # Check the number of columns
  expect_equal(ncol(excesses), ncol(df_lags) + 2)
  # Check the number of rows
  expect_equal(nrow(excesses), nrow(df_lags))

  # for the hnorm == 0  and tau == 0 ie same site (si, si)
  k_h_tau <- excesses[excesses$hnorm == 0 & excesses$tau == 0, ]
  n_marg <- get_marginal_excess(data_rain, q)
  # check if it is the same
  expect_equal(unique(k_h_tau$kij), n_marg)

})

test_that("theorical_chi", {
  tau_vect <- 0:10
  data_rain <- matrix(runif(100 * 25, 0, 100), nrow = 100, ncol = 25)
  data_rain <- as.data.frame(data_rain)
  sites_coords <- generate_grid_coords(5) # Example sites_coords
  true_param <- c(0.4, 0.2, 1.5, 1)
  df_lags <- get_lag_vectors(sites_coords, true_param,
                            hmax = sqrt(17), tau_vect = tau_vect)

  chi_theorical <- theorical_chi(true_param, df_lags)

  # for one hnorm and one tau
  hnorm <- 1
  tau <- 3
  semivar <- true_param[1]*hnorm^true_param[3] + true_param[2]*tau^true_param[4]
  chi_h_t_verif <- 2 * (1 - pnorm(sqrt(semivar)))

  chi_h_t <- chi_theorical$chi[chi_theorical$hnorm == hnorm &
                                    chi_theorical$tau == tau]

  expect_equal(unique(chi_h_t), chi_h_t_verif)
})
