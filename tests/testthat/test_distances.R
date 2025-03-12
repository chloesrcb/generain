test_that("reshape_distances produces correct format", {
    # Calculate the spatial coordinates of the regular grid
    sites_coords <- generate_grid_coords(5)
    adv <- c(0, 0) # Advection coordinates
    # Calculate the distance matrix with advanced parameters
    dist_mat_adv <- get_dist_mat(sites_coords, dmax = 100,
                                latlon = FALSE)
    df_dist <- reshape_distances(dist_mat_adv) # reshape the distance matrix

    # Check if the required columns are present in the reshaped distance matrix
    expect_true(all(c("Y", "X", "value") %in% colnames(df_dist)))

    # Check the number of unique values in the 'tau', 'Y', and 'X' columns
    expect_equal(length(unique(df_dist$Y)), 25)
    expect_equal(length(unique(df_dist$X)), 25)
})



test_that("get_h_vect test", {
  # Calculate the spatial coordinates of the regular grid
  sites_coords <- generate_grid_coords(5)
  # Calculate the distance matrix with advanced parameters
  dist_mat <- get_dist_mat(sites_coords, dmax = 100,
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
  df_lags <- get_lag_vectors(sites_coords, tau_max = 10, latlon=FALSE)

  nsites <- nrow(sites_coords)
  pairs_comb <- combn(nsites, 2)
  diag_pairs <- matrix(rep(1:nsites, each = 2), nrow = 2)
  indices <- cbind(pairs_comb, diag_pairs)
  n_pairs <- ncol(indices)
  expect_equal(nrow(df_lags), n_pairs * length(tau_vect))

  # Check col
  expect_equal(ncol(df_lags), 10)
  # Check taus
  expect_equal(unique(df_lags$tau), tau_vect)

  # Check value in hnorm
  dist_mat <- get_dist_mat(sites_coords, latlon = FALSE)
  df_dist <- reshape_distances(dist_mat)
  h_vect <- get_h_vect(df_dist)
  hnorm_no0 <- unique(df_lags$hnorm[df_lags$hnorm != 0])
  expect_true(all(hnorm_no0 %in% h_vect))
})



test_that("get_conditional_lag_vectors", {
  ngrid <- 5
  sites_coords <- generate_grid_coords(ngrid)
  s0 <- c(1, 1)
  t0 <- 1

  # Calculate the h vector
  true_param <- c(0.4, 0.2, 1.5, 1)
  hmax <- sqrt(17)
  tau_vect <- 0:10

  df_lags <- get_conditional_lag_vectors(sites_coords,
                          tau_max = max(tau_vect), s0 = s0, t0 = t0)

  # check row
  expect_equal(nrow(df_lags), ngrid^2 * length(tau_vect))
})


test_that("haversine_distance_with_advection works correctly", {
  lat1 <- 43.62505
  lon1 <- 3.862038
  lat2 <- 43.68333
  lon2 <- 4.133333
  adv <- c(0.0, 0.0)
  tau <- 0
  expected_distance_km <- 22.77

  distance_with_adv <- haversine_distance_with_advection(lat1, lon1, lat2, lon2,
                        adv, tau)$distance
  expect_equal(round(distance_with_adv, 2), expected_distance_km)

  theta_dir <- haversine_distance_with_advection(lat1, lon1, lat2, lon2,
                        adv, tau)$theta
  theta_deg <- theta_dir * 180 / pi
  dir_card <- convert_to_cardinal(theta_deg)
  theta_meteo <- (pi/2 - theta_dir) %% (2 * pi)
  theta_meteo_deg <- theta_meteo * 180 / pi
  dir_card_meteo <- convert_to_cardinal(theta_meteo_deg)
  expect_equal(dir_card_meteo, "E")
})
