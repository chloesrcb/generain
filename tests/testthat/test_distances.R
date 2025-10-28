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

# ======================================================================
# TESTS FOR get_dist_mat()
# ======================================================================

test_that("get_dist_mat computes Euclidean distances correctly", {
  # Simple 3-point grid
  coords <- data.frame(
    Longitude = c(0, 1, 2),
    Latitude = c(0, 1, 2)
  )
  rownames(coords) <- c("A", "B", "C")

  dist_mat <- get_dist_mat(coords, latlon = FALSE)

  # Check matrix structure
  expect_true(is.matrix(dist_mat))
  expect_equal(dim(dist_mat), c(3, 3))
  expect_true(all(rownames(dist_mat) == colnames(dist_mat)))

  # Verify Euclidean distances
  expect_equal(dist_mat["A", "B"], sqrt(2), tolerance = 1e-6)
  expect_equal(dist_mat["A", "C"], sqrt(8), tolerance = 1e-6)
  expect_equal(dist_mat["B", "C"], sqrt(2), tolerance = 1e-6)

  # Diagonal must be zero
  expect_true(all(diag(dist_mat) == 0))
})


test_that("get_dist_mat assigns names from Station column if available", {
  coords <- data.frame(
    Station = c("S1", "S2", "S3"),
    Longitude = c(1, 2, 3),
    Latitude = c(4, 5, 6)
  )

  dist_mat <- get_dist_mat(coords, latlon = FALSE)

  expect_true(all(c("S1", "S2", "S3") %in% rownames(dist_mat)))
  expect_true(all(rownames(dist_mat) == colnames(dist_mat)))
})


test_that("get_dist_mat handles missing Station but uses rownames", {
  coords <- data.frame(
    Longitude = c(0, 1, 2),
    Latitude = c(0, 1, 2)
  )
  rownames(coords) <- c("P1", "P2", "P3")

  dist_mat <- get_dist_mat(coords, latlon = FALSE)

  expect_equal(rownames(dist_mat), c("P1", "P2", "P3"))
  expect_equal(colnames(dist_mat), c("P1", "P2", "P3"))
})


test_that("get_dist_mat applies distance threshold (dmax) correctly", {
  coords <- data.frame(
    Longitude = c(0, 1, 2),
    Latitude = c(0, 0, 0)
  )
  rownames(coords) <- c("X", "Y", "Z")

  # Distance X-Z = 2 > 1.5, should be set to 0
  dist_mat <- get_dist_mat(coords, dmax = 1.5, latlon = FALSE)

  expect_equal(dist_mat["X", "Z"], 0)
  expect_equal(dist_mat["X", "Y"], 1)
})


test_that("get_dist_mat computes haversine distances (latlon = TRUE)", {
  coords <- data.frame(
    Station = c("A", "B"),
    Latitude = c(0, 0),
    Longitude = c(0, 1)
  )

  dist_mat <- get_dist_mat(coords, latlon = TRUE)

  # Haversine: distance at equator for 1° longitude ≈ 111.3 km
  expect_true(all(dim(dist_mat) == c(2, 2)))
  expect_true(all(diag(dist_mat) == 0))
  expect_true(dist_mat["A", "B"] > 100000 && dist_mat["A", "B"] < 120000)
})


test_that("get_dist_mat throws error for invalid coordinates", {
  coords <- data.frame(x = 1:3, y = 1:3)

  expect_error(get_dist_mat(coords, latlon = FALSE),
               "must have columns 'Latitude' and 'Longitude'")
})
