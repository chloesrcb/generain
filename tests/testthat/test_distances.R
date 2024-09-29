test_that("reshape_distances produces correct format", {
    # Calculate the spatial coordinates of the regular grid
    sites_coords <- generate_grid_coords(5)
    adv <- c(0, 0) # Advection coordinates
    # Calculate the distance matrix with advanced parameters
    dist_mat_adv <- get_dist_mat(sites_coords, adv = adv, tau = 1:10,
                                latlon = FALSE)
    df_dist <- reshape_distances(dist_mat_adv) # reshape the distance matrix

    # Check if the required columns are present in the reshaped distance matrix
    expect_true(all(c("Y", "X", "value") %in% colnames(df_dist)))

    # Check the number of unique values in the 'tau', 'Y', and 'X' columns
    expect_equal(length(unique(df_dist$Y)), 25)
    expect_equal(length(unique(df_dist$X)), 25)
})



test_that("get_h_vect without adv", {
  # Calculate the spatial coordinates of the regular grid
  sites_coords <- generate_grid_coords(5)
  # Calculate the distance matrix with advanced parameters
  dist_mat <- get_dist_mat(sites_coords, adv = c(0, 0), tau = 1:10,
                           latlon = FALSE)
  # Reshape the distance matrix
  df_dist <- reshape_distances(dist_mat)

  # Calculate the h vector
  h_vect <- get_h_vect(df_dist)
  d_max <- max(df_dist$value) # Maximum distance in the distance matrix
  expect_true(length(h_vect) == length(unique(h_vect))) # Check for duplicates
  expect_true(all(h_vect > 0)) # Check if all values are positive
  expect_true(all(h_vect <= d_max)) # Check if all values are less than d_max
  expect_true(all(sort(h_vect) == h_vect)) # Check if the values are sorted

  hmax <- 3 # Maximum value for h
  h_vect <- get_h_vect(df_dist, hmax = hmax) # Calculate the h vector with hmax
  # Check if all values are within the range [1, hmax]
  expect_true(all(h_vect >= 1 & h_vect <= hmax)) 
  expect_true(all(sort(h_vect) == h_vect)) # Check if the values are sorted
})

test_that("get_h_vect with adv", {
  # Calculate the spatial coordinates of the regular grid
  sites_coords <- generate_grid_coords(5)
  adv <- c(1, 1) # Advection coordinates
  # Calculate the distance matrix with advanced parameters
  dist_mat <- get_dist_mat(sites_coords, adv = adv, tau = 1:10,
                           latlon = FALSE)
  # Reshape the distance matrix
  df_dist <- reshape_distances(dist_mat)

  # Calculate the h vector
  h_vect <- get_h_vect(df_dist)
  d_max <- max(df_dist$value) # Maximum distance in the distance matrix
  expect_true(length(h_vect) == length(unique(h_vect))) # Check for duplicates
  expect_true(all(h_vect > 0)) # Check if all values are positive
  expect_true(all(h_vect <= d_max)) # Check if all values are less than d_max
  expect_true(all(sort(h_vect) == h_vect)) # Check if the values are sorted

  hmax <- 3 # Maximum value for h
  h_vect <- get_h_vect(df_dist, hmax = hmax) # Calculate the h vector with hmax
  # Check if all values are within the range [1, hmax]
  expect_true(all(h_vect >= 1 & h_vect <= hmax))
  expect_true(all(sort(h_vect) == h_vect)) # Check if the values are sorted
})


test_that("get_lag_vectors", {
  ngrid <- 5
  sites_coords <- generate_grid_coords(ngrid)

  # Calculate the h vector
  true_param <- c(0.4, 0.2, 1.5, 1)
  hmax <- sqrt(17)
  tau_vect <- 0:10

  # Check row without hmax
  df_lags <- get_lag_vectors(sites_coords, true_param,
                          hmax = NA, tau_vect)

  nsites <- nrow(sites_coords)
  pairs_comb <- combn(nsites, 2)
  diag_pairs <- matrix(rep(1:nsites, each = 2), nrow = 2)
  indices <- cbind(pairs_comb, diag_pairs)
  n_pairs <- ncol(indices)
  expect_equal(nrow(df_lags), n_pairs * length(tau_vect))

  # Check col
  expect_equal(ncol(df_lags), 6)
  # Check taus
  expect_equal(unique(df_lags$tau), tau_vect)
  # Check hmax
  df_lags <- get_lag_vectors(sites_coords, true_param,
                          hmax, tau_vect)
  expect_true(max(df_lags$hnorm) <= hmax)

  # Check value in hnorm
  dist_mat <- get_dist_mat(sites_coords, adv = c(0, 0), tau = 0:10,
                           latlon = FALSE)
  df_dist <- reshape_distances(dist_mat)
  h_vect <- get_h_vect(df_dist)
  hnorm_no0 <- unique(df_lags$hnorm[df_lags$hnorm != 0])
  expect_true(all(hnorm_no0 %in% h_vect))

  df_lags_adv0 <- get_lag_vectors(sites_coords, c(true_param, 0, 0),
                          hmax, tau_vect)

  expect_true(all(df_lags$hnorm %in% df_lags_adv0$hnorm))
})


test_that("get_lag_vectors with advection", {
  ngrid <- 5
  sites_coords <- generate_grid_coords(ngrid)

  # Calculate the h vector
  true_param <- c(0.4, 0.2, 1.5, 1)
  hmax <- sqrt(17)
  tau_vect <- 0:10

  # with advection
  adv <- c(1, 1)
  true_param <- c(true_param, adv)

  df_lags <- get_lag_vectors(sites_coords, true_param,
                          hmax, tau_vect)

  expect_equal(max(df_lags$hnorm), hmax)

  # Check value in hnorm
  dist_mat <- get_dist_mat(sites_coords, adv = adv, tau = tau_vect,
                           latlon = FALSE)
  df_dist <- reshape_distances(dist_mat)
  h_vect <- get_h_vect(df_dist, hmax = hmax)
  hnorm_no0 <- sort(unique(df_lags$hnorm[df_lags$hnorm != 0]))
  expect_true(all(hnorm_no0 %in% h_vect))

  # On simulation data
  ngrid <- 5
  sites_coords <- generate_grid_coords(ngrid)
  temp <- 1:300

  # Number of realizations
  nres <- 10

  # Folder
  foldername <- paste0("./data/simulations_BR/sim_25s_300t/")

  list_simu <- list()
  for (i in 1:nres) {
    file_name <- paste0(foldername, "br_", ngrid^2, "s_",
                          length(temp), "t_", i, ".csv")
    list_simu[[i]] <- read.csv(file_name)
  }

  simu_df <- list_simu[[1]]

  true_param <- c(0.4, 0.2, 1.5, 1)

  df_lags <- get_lag_vectors(sites_coords, true_param,
                          hmax = NA, tau_vect = 0:10)

  adv <- c(0.1, 0.)
  df_lags_adv <- get_lag_vectors(sites_coords, c(true_param, adv),
                          hmax = NA, tau_vect = 0:10)

  # check dimension
  expect_equal(nrow(df_lags_adv), nrow(df_lags))

  same_index <- which(df_lags_adv$tau == df_lags$tau[11] &
                      df_lags_adv$s1 == df_lags$s1[11] &
                      df_lags_adv$s2 == df_lags$s2[11])
  hxadv <- df_lags$hx[11] - adv[1] * df_lags$tau[11]
  expect_equal(df_lags_adv$hx[same_index], hxadv)


  adv <- c(0.5, 0.8)
  df_lags_adv <- get_lag_vectors(sites_coords, c(true_param, adv),
                          hmax = NA, tau_vect = 0:10)

  # check dimension
  expect_equal(nrow(df_lags_adv), nrow(df_lags))

  same_index <- which(df_lags_adv$tau == df_lags$tau[11] &
                      df_lags_adv$s1 == df_lags$s1[11] &
                      df_lags_adv$s2 == df_lags$s2[11])
  hxadv <- df_lags$hx[11] - adv[1] * df_lags$tau[11]
  expect_equal(df_lags_adv$hx[same_index], hxadv)

})


test_that("get_conditional_lag_vectors with advection", {
  ngrid <- 5
  sites_coords <- generate_grid_coords(ngrid)
  s0 <- c(1, 1)
  t0 <- 1

  # Calculate the h vector
  true_param <- c(0.4, 0.2, 1.5, 1)
  hmax <- sqrt(17)
  tau_vect <- 0:10

  # with advection
  adv <- c(1, 1)
  true_param <- c(true_param, adv)

  df_lags <- get_conditional_lag_vectors(sites_coords, true_param,
                          hmax = NA, tau_vect = tau_vect, s0 = s0, t0 = t0)

  # check row
  tau_lag <- sort(unique(abs(tau_vect - t0)))
  expect_equal(nrow(df_lags), ngrid^2 * length(tau_lag))

  # check with hmax
  df_lags <- get_conditional_lag_vectors(sites_coords, true_param,
                          hmax = hmax, tau_vect = tau_vect, s0 = s0, t0 = t0)
  expect_true(max(df_lags$hnorm) <= hmax)

  # Check value in hnorm
  dist_mat <- get_dist_mat(sites_coords, adv = adv, tau = tau_vect,
                           latlon = FALSE)
  df_dist <- reshape_distances(dist_mat)
  h_vect <- get_h_vect(df_dist, hmax = hmax)
  hnorm_no0 <- sort(unique(df_lags$hnorm[df_lags$hnorm != 0]))
  expect_true(all(hnorm_no0 %in% h_vect))
})
