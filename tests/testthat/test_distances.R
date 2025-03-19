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
  df_lags <- get_lag_vectors(sites_coords, tau_vect, latlon = FALSE)

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
  tau_vect <- -10:10

  df_lags <- get_conditional_lag_vectors(sites_coords,
                          tau_vect, s0 = s0, t0 = t0)

  # check row
  expect_equal(nrow(df_lags), ngrid^2 * length(tau_vect))
})


test_that("haversine_distance_with_advection works correctly", {
  # Real-world lat/lon coordinates for Montpellier and Lunel
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
                        adv, tau)$theta_meteo
  theta_meteo_deg <- theta_dir * 180 / pi
  dir_card_meteo <- convert_to_cardinal(theta_meteo_deg)
  expect_equal(dir_card_meteo, "E")

  # Real-world lat/lon coordinates for Paris and London
  lat1 <- 48.8566 # Paris
  lon1 <- 2.3522
  lat2 <- 51.5074 # London
  lon2 <- -0.1278
  adv <- c(0.0, 0.0)
  tau <- 0
  expected_distance_km <- 343.51 # Approx distance in km

  # Expected values assuming haversine distance
  # Angle from Paris to London
  expected_theta_PL <- atan2(48.8566 - 51.5074, 2.3522 - (-0.1278))
  expected_theta_deg_PL <- expected_theta_PL * 180 / pi
  expected_x_polar_PL <- expected_distance_km * cos(expected_theta_PL)
  expected_y_polar_PL <- expected_distance_km * sin(expected_theta_PL)
  # Angle from London to Paris
  expected_theta_LP <- atan2(51.5074 - 48.8566, -0.1278 - 2.3522)
  expected_theta_deg_LP <- expected_theta_LP * 180 / pi
  expected_x_polar_LP <- expected_distance_km * cos(expected_theta_LP)
  expected_y_polar_LP <- expected_distance_km * sin(expected_theta_LP)

  distance_LP <- haversine_distance_with_advection(lat1, lon1, lat2, lon2,
                        adv, tau)$distance
  distance_PL <- haversine_distance_with_advection(lat2, lon2, lat1, lon1,
                        adv, tau)$distance
  expect_equal(round(distance_LP, 0), round(expected_distance_km, 0))
  expect_equal(round(distance_PL, 0), round(expected_distance_km, 0))
  expect_equal(round(distance_LP, 0), round(distance_PL, 0))

  theta_LP <- haversine_distance_with_advection(lat1, lon1, lat2, lon2,
                        adv, tau)$theta
  theta_PL <- haversine_distance_with_advection(lat2, lon2, lat1, lon1,
                        adv, tau)$theta
  expect_equal(round(theta_LP, 2), round(expected_theta_LP, 2))
  expect_equal(round(theta_PL, 2), round(expected_theta_PL, 2))
  expect_false(theta_LP == theta_PL)

  x_polar_LP <- distance_LP * cos(theta_LP)
  y_polar_LP <- distance_LP * sin(theta_LP)
  x_polar_PL <- distance_PL * cos(theta_PL)
  y_polar_PL <- distance_PL * sin(theta_PL)

  expect_equal(round(x_polar_LP, 0), round(expected_x_polar_LP, 0))
  expect_equal(round(y_polar_LP, 0), round(expected_y_polar_LP, 0))
  expect_equal(round(x_polar_PL, 0), round(expected_x_polar_PL, 0))
  expect_equal(round(y_polar_PL, 0), round(expected_y_polar_PL, 0))

})


test_that("haversine_distance_with_advection accounts for advection correctly", {
  lat1 <- 48.8566 # Paris
  lon1 <- 2.3522
  lat2 <- 51.5074 # London
  lon2 <- -0.1278

  # No advection case
  adv_0 <- c(0.0, 0.0)
  tau_0 <- 0

  # With advection: 1 m/s East (longitude) and 1 m/s North (latitude), for 100s
  adv <- c(1, 1)
  tau <- 100

  # Compute results
  result_no_adv <- haversine_distance_with_advection(lat1, lon1, lat2, lon2,
                                                                  adv_0, tau_0)
  result_with_adv <- haversine_distance_with_advection(lat1, lon1, lat2, lon2,
                                                                      adv, tau)

  # Expected base distance
  expected_distance_km <- result_no_adv$distance

  # Compute new lat/lon due to advection
  adv_x <- adv[1] * tau  # Advection in X (meters)
  adv_y <- adv[2] * tau  # Advection in Y (meters)

  # Adjusted lat/lon
  lon2_adj <- lon2 - adv_x / (111319 * cos(lat2 * pi / 180))
  lat2_adj <- lat2 - adv_y / 111319

  # Compute expected distance using Haversine again
  expected_adjusted_distance_km <- haversine_distance_with_advection(lat1,
                              lon1, lat2_adj, lon2_adj, adv_0, tau_0)$distance

  # Check that the baseline distance matches expected
  expect_equal(round(result_no_adv$distance, 2), round(expected_distance_km, 2))

  # Check that the adjusted distance is within reasonable bounds of
  # the recomputed expected distance
  expect_equal(round(result_with_adv$distance, 0),
                          round(expected_adjusted_distance_km, 0))
})



test_that("haversine_distance_with_advection  computes correct values with different etas", {
  wind_vect <- c(-9.8, 8)  # Wind effect for a simple test case
  wind_vect_kmh <- wind_vect  # Convert wind to km/h

  eta1 <- 0.5
  eta2 <- 0.1
  adv1 <- (abs(wind_vect_kmh)^eta1) * sign(wind_vect_kmh) * eta2  # in km/h
  adv2 <- (abs(wind_vect_kmh)^eta2) * sign(wind_vect_kmh) * eta1  # in km/h
  # Real-world lat/lon coordinates for Paris and close to Paris
  lat1 <- 48.8566 # Paris
  lon1 <- 2.3522
  lat2 <- 48.8566 # Close to Paris
  lon2 <- 2.3523

  result_no_adv <- haversine_distance_with_advection(lat1, lon1, lat2, lon2,
                                                      c(0, 0), 0)
  # First case: eta1 = 0.5, eta2 = 0.1
  result1 <- haversine_distance_with_advection(lat1, lon1, lat2, lon2,
                                                adv1, 1)
  # Second case: eta1 = 0.1, eta2 = 0.5
  result2 <- haversine_distance_with_advection(lat1, lon1, lat2, lon2,
                                                adv2, 1)

  distance_no_adv <- result_no_adv$distance
  distance1 <- result1$distance
  distance2 <- result2$distance

  # Check that the baseline distance is the same
  expect_false(distance_no_adv == distance1)
  expect_false(distance_no_adv == distance2)
  expect_false(distance1 == distance2)
})
