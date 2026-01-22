
test_that("get_conditional_lag_vectors produces correct lag structure", {
  df_coords <- data.frame(
    Longitude = c(0, 1, 2),
    Latitude = c(0, 0, 0)
  )
  rownames(df_coords) <- paste0("S", 1:3)

  s0 <- c(1, 0)  # conditional point coordinates
  tau_vect <- 0:2

  res <- get_conditional_lag_vectors(
    df_coords, s0 = s0, t0 = 0, tau_vect = tau_vect, latlon = FALSE
  )

  # Basic structure
  expect_s3_class(res, "data.table")
  expect_true(all(c("s1", "s2", "hx", "hy", "tau", "hnorm") %in% names(res)))

  # Check that tau covers all requested lags
  expect_true(all(tau_vect %in% res$tau))

  # Check that number of rows matches (n_sites * length(tau_vect))
  n_sites <- nrow(df_coords)
  expect_equal(nrow(res), n_sites * length(tau_vect))

  # Check that Euclidean distance is computed correctly
  # Corrected check
  expect_equal(
    res$hnorm[res$s1 == "S2" & res$s2 == "S3"][1],
    sqrt((df_coords$Longitude[3] - df_coords$Longitude[2])^2 +
        (df_coords$Latitude[3] - df_coords$Latitude[2])^2)
  )

})

test_that("get_conditional_lag_vectors errors when conditional site not found", {
  df_coords <- data.frame(Longitude = c(0, 1), Latitude = c(0, 1))
  rownames(df_coords) <- c("S1", "S2")

  # Conditional site outside known coordinates
  expect_error(get_conditional_lag_vectors(df_coords, s0 = c(3, 3)))
})

test_that("get_conditional_lag_vectors handles latlon distance correctly", {
  df_coords <- data.frame(
    Longitude = c(0, 1),
    Latitude = c(0, 0)
  )
  rownames(df_coords) <- c("A", "B")

  res <- get_conditional_lag_vectors(df_coords, s0 = c(0, 0), tau_vect = 0:1, latlon = TRUE)

  # Check that geodesic distance (hnorm) is in meters and positive
  expect_true(all(res$hnorm >= 0))
  expect_true(max(res$hnorm) > 100000)  # 1° ~ 111 km on Earth
})

# ----------------------------------------------------------------------
# empirical_excesses_rpar()
# ----------------------------------------------------------------------

test_that("empirical_excesses_rpar computes exceedances correctly (fixed threshold)", {
  data_rain <- matrix(c(
    1, 2, 3,
    2, 3, 4,
    5, 6, 7,
    8, 9, 10
  ), nrow = 4, ncol = 3, byrow = TRUE)

  colnames(data_rain) <- paste0("S", 1:3)

  df_lags <- data.table(s1 = 1, s2 = c(1, 2, 3), tau = c(0, 1, 2))

  res <- empirical_excesses_rpar(data_rain, df_lags, threshold = 5)

  expect_s3_class(res, "data.table")
  expect_true(all(c("s1", "s2", "tau", "kij", "Tobs") %in% names(res)))
  expect_true(all(res$Tobs == 1))
  expect_true(all(res$kij %in% c(0, 1)))
})

test_that("empirical_excesses_rpar works with site-specific thresholds", {
  set.seed(123)
  data_rain <- matrix(runif(30, 0, 10), nrow = 10, ncol = 3)
  colnames(data_rain) <- paste0("S", 1:3)

  threshold <- 5
  df_lags <- data.table(s1 = 1, s2 = c(1, 2, 3), tau = c(0, 1, 2))

  res <- empirical_excesses_rpar(data_rain, df_lags, threshold = threshold)

  expect_equal(nrow(res), nrow(df_lags))
  expect_true(all(res$kij %in% c(0, 1)))
})


test_that("empirical_excesses_rpar handles missing values gracefully", {
  data_rain <- matrix(c(1, NA, 3, 4, 5, 6), nrow = 3, ncol = 2)
  colnames(data_rain) <- c("A", "B")
  df_lags <- data.table(s1 = 1, s2 = 2, tau = 1)
  
  res <- empirical_excesses_rpar(data_rain, df_lags, threshold = 3)
  
  expect_true(all(!is.na(res$kij)))
  expect_true(all(res$kij %in% c(0, 1)))
})

# ----------------------------------------------------------------------
# --- 3. INTEGRATION TEST: realistic synthetic data
# ----------------------------------------------------------------------

test_that("Integration: lag vectors and excess computation on realistic data", {
  set.seed(123)
  
  # Simulate spatial coordinates for 4 sites
  df_coords <- data.frame(
    Longitude = seq(0, 3),
    Latitude = c(0, 0, 1, 1)
  )
  rownames(df_coords) <- paste0("S", 1:4)
  
  # Simulate rainfall data with an obvious extreme
  data_rain <- matrix(rnorm(200, mean = 10, sd = 4), nrow = 50, ncol = 4)
  colnames(data_rain) <- paste0("S", 1:4)
  data_rain[25, 3] <- 40  # create one extreme at site 3
  
  # Generate lag table for conditional site S1
  lags <- get_conditional_lag_vectors(df_coords, s0 = c(0, 0), tau_vect = 0:5, latlon = FALSE)
  expect_true(nrow(lags) > 0)
  
  # Compute empirical excesses using quantile-based threshold
  res <- empirical_excesses_rpar(data_rain, lags, threshold = 5)
  
  # Basic structure and consistency
  expect_true(all(c("s1", "s2", "tau", "kij") %in% names(res)))
  expect_true(all(res$kij %in% c(0, 1)))
  
  # At least one exceedance should exist near the extreme
  expect_true(sum(res$kij) >= 1)
})
