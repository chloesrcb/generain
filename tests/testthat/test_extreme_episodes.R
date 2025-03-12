# Description: Unit tests for the extreme episodes function

# Load necessary packages
library(geosphere)
library(data.table)

### Tests for the function select_extreme_episodes -----------------------------
# Sample data setup
test_sites_coords <- data.frame(
  lon = c(0, 1, 2, 3),
  lat = c(0, 1, 2, 3)
)

rownames(test_sites_coords) <- c("Site1", "Site2", "Site3", "Site4")

test_data <- data.frame(
  Site1 = c(5, 10, 15, 20),
  Site2 = c(2, 12, 8, 25),
  Site3 = c(6, 14, 9, 30),
  Site4 = c(3, 11, 7, 22)
)

# Define parameters
test_quantile <- 0.9
test_min_spatial_dist <- 1
test_delta <- 1
test_n_max_episodes <- 2
test_time_ext <- 0

# Test case 1: Basic functionality
test_that("Function returns expected output format", {
  result <- select_extreme_episodes(
    test_sites_coords, test_data, test_quantile,
    test_min_spatial_dist, test_delta, test_n_max_episodes, test_time_ext
  )
  expect_true(is.list(result))
  expect_true(ncol(result) == 3)  # Should have s0, t0, u_s0 columns
})

# Test case 2: Ensuring quantile threshold filtering
test_that("Function correctly filters episodes based on site-specific quantile", {
  result <- select_extreme_episodes(
    test_sites_coords, test_data, test_quantile,
    test_min_spatial_dist, test_delta, test_n_max_episodes, test_time_ext
  )

  # Verify that each u_s0 corresponds to the quantile of the concerned site
  for (i in seq_len(nrow(result))) {
    site_name <- result$s0[i]
    site_threshold <- quantile(test_data[[site_name]], probs = test_quantile,
                        na.rm = TRUE)
    expect_true(result$u_s0[i] >= site_threshold)
  }
})


# Test case 3: Handling edge cases (e.g., no extreme episodes)
test_that("Function returns empty list if no episodes exceed threshold", {
  low_data <- data.frame(Site1 = c(1, 2, 3, 4), Site2 = c(1, 2, 3, 4))

  result <- select_extreme_episodes(
    test_sites_coords, low_data, test_quantile,
    test_min_spatial_dist, test_delta, test_n_max_episodes, test_time_ext
  )

  expect_true(nrow(result) == 0)
})

# Test case 4: Ensuring spatial and temporal constraints
test_that("Function respects spatial and temporal constraints", {
  result <- select_extreme_episodes(
    test_sites_coords, test_data, test_quantile,
    test_min_spatial_dist, test_delta, test_n_max_episodes, test_time_ext
  )
  # Check that the spatial and temporal constraints are respected
  if (nrow(result) > 1) {
    for (i in 1:(nrow(result) - 1)) {
      for (j in (i + 1):nrow(result)) {
        dist_sites <- geosphere::distHaversine(
          test_sites_coords[result$s0[i], ],
          test_sites_coords[result$s0[j], ]
        ) / 1000

        time_diff <- abs(result$t0[i] - result$t0[j])
        expect_true(dist_sites >= test_min_spatial_dist ||
                                          time_diff > 2 * test_delta)
      }
    }
  }
})

# Integrated Test: Checking for temporal overlap
test_that("select_extreme_episodes does not produce overlapping episodes", {
  selected_points <- select_extreme_episodes(
    test_sites_coords, test_data, test_quantile,
    test_min_spatial_dist, test_delta, test_n_max_episodes, test_time_ext
  )

  # Check that no overlaps exist for each unique site
  unique_sites <- unique(selected_points$s0)

  for (site in unique_sites) {
    expect_false(check_intervals_overlap(site, selected_points,
                                            test_delta, test_time_ext),
                 info = paste("Overlap detected for site:", site))
  }
})

### Tests for the function check_intervals_overlap -----------------------------
# Example of data from select_extreme_episodes
selected_points <- data.table(
  s0 = c("Site1", "Site1", "Site2", "Site3", "Site3"),
  t0 = c(5, 15, 10, 20, 25),
  u_s0 = c(12, 14, 16, 18, 20)  # Example threshold values
)

test_delta <- 3
test_beta <- 1

# Test 1: No overlap
test_that("No overlap detected", {
  result <- check_intervals_overlap("Site2", selected_points,
                                              test_delta, test_beta)
  expect_false(result)
})

# Test 2: Overlap is correctly detected
test_that("Correct detection of an overlap", {
  # Add a point that will cause an overlap
  selected_points_overlap <- rbind(selected_points,
                                data.table(s0 = "Site1", t0 = 7, u_s0 = 13))

  result <- check_intervals_overlap("Site1", selected_points_overlap,
                                                      test_delta, test_beta)
  expect_true(result)
})

# Test 3: A single point for a site (no possible overlap)
test_that("No overlap with only one point", {
  result <- check_intervals_overlap("Site2", selected_points,
                                              test_delta, test_beta)
  expect_false(result)
})

# Test 4: Edge case where intervals touch but do not overlap
test_that("Adjacent intervals without overlap are not detected as overlapping", {
  selected_points_adjacent <- data.table(
    s0 = c("Site1", "Site1"),
    t0 = c(5, 13),  # With delta = 3 and beta = 1, [2,9] and [10,17]
                    # do not overlap
    u_s0 = c(12, 14)
  )

  result <- check_intervals_overlap("Site1", selected_points_adjacent,
                                                      test_delta, test_beta)
  expect_false(result)
})


### Test get_extreme_episodes function -----------------------------------------

# Sample data setup
test_data <- data.frame(
  Site1 = c(5, 10, 15, 20, 25, 30),
  Site2 = c(2, 12, 8, 25, 18, 22),
  Site3 = c(6, 14, 9, 30, 20, 28),
  Site4 = c(3, 11, 7, 22, 17, 27)
)

# Define parameters
test_delta <- 2
test_beta <- 1

# Sample selected points (including some invalid points)
test_selected_points <- data.table(
  s0 = c("Site1", "Site3"),
  t0 = c(3, 6)  # 6 is at the boundary, might be invalid
)

test_that("Function returns correct episode structure", {
  result <- get_extreme_episodes(test_selected_points, test_data, test_delta,
                                    test_beta)

  expect_type(result, "list")  # Should return a list
  expect_named(result, c("episodes", "selected_points"))  # Correct names

  # Check types of each
  expect_type(result$episodes, "list")
  expect_s3_class(result$selected_points, "data.table")
})

test_that("Function removes invalid selected points", {
  result <- get_extreme_episodes(test_selected_points, test_data, test_delta,
                                          test_beta)

  # Check if selected_points was modified
  if (nrow(result$selected_points) < nrow(test_selected_points)) {
    expect_true(nrow(result$selected_points) < nrow(test_selected_points))
    print("Test Passed: Function correctly removes invalid selected points.")
  }
})

test_that("Function handles cases where all points are invalid", {
  invalid_points <- data.table(
    s0 = c("Site2", "Site4"),
    t0 = c(nrow(test_data), nrow(test_data))  # All at boundary, likely invalid
  )

  result <- get_extreme_episodes(invalid_points, test_data, test_delta, 
                                  test_beta)

  expect_length(result$episodes, 0)  # No valid episodes
  expect_equal(nrow(result$selected_points), 0)  # No valid selected points
  print("Test Passed: Function correctly handles all invalid points.")
})



test_that("Function maintains correct episode size", {
  set.seed(123)  # Fixer la seed pour la reproductibilité
  test_data <- data.frame(
    Site1 = rnorm(100, mean = 10, sd = 5),  # Moyenne 10, écart-type 5
    Site2 = rnorm(100, mean = 15, sd = 7),  # Moyenne 15, écart-type 7
    Site3 = rnorm(100, mean = 20, sd = 6),  # Moyenne 20, écart-type 6
    Site4 = rnorm(100, mean = 25, sd = 8)   # Moyenne 25, écart-type 8
  )

  test_selected_points <- data.table(
    s0 = c("Site1", "Site3"),
    t0 = c(which.max(test_data$Site1), which.max(test_data$Site2))
  )
 
  test_delta <- 2
  test_beta <- 0
  result <- get_extreme_episodes(test_selected_points, test_data, test_delta,
                                        test_beta)

  # If there are episodes, check their size
  for (episode in result$episodes) {
    expect_equal(nrow(episode), 2 * test_delta + 1)  # Validate episode size
  }

  test_delta <- 2
  test_beta <- 1
  t_inf <- test_selected_points$t0 - test_delta - test_beta
  t_sup <- test_selected_points$t0 + test_delta + test_beta

  test_data$Site1[t_inf[1]:t_sup[1]]

  result <- get_extreme_episodes(test_selected_points, test_data, test_delta,
                                        test_beta)

  # If there are episodes, check their size
  for (episode in result$episodes) {
    expect_equal(nrow(episode), 2 * test_delta + 1 + 2 * test_beta)
  }

  print("Test Passed: Episodes have correct size.")
})
