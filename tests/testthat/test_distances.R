# Test that reshape_distances produces correct format
test_that("reshape_distances produces correct format", {
    # Calculate the spatial coordinates of the regular grid
    sites_coords <- generate_grid_coords(5)
    avd <- c(0, 0) # Advection coordinates
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
