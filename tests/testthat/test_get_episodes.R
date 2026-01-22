
# -----------------------------------------------------------------------
# get_spatiotemp_excess()
# -----------------------------------------------------------------------

test_that("get_spatiotemp_excess works with quantile", {
  set.seed(123)
  data <- matrix(runif(100, 0, 10), nrow = 10, ncol = 10)
  colnames(data) <- paste0("S", 1:10)

  res <- get_spatiotemp_excess(data, quantile = 0.9)
  expect_type(res, "list")
  expect_true(all(c("list_s", "list_t", "list_u") %in% names(res)))
  expect_true(length(res$list_s) > 0)
  expect_true(length(res$list_t) == length(res$list_s))
  expect_true(length(res$list_u) == length(res$list_s))
})

test_that("get_spatiotemp_excess works with fixed threshold", {
  data <- matrix(runif(50, 0, 10), nrow = 10, ncol = 5)
  colnames(data) <- paste0("P", 1:5)
  thr <- 5

  res <- get_spatiotemp_excess(data, threshold = thr)
  expect_true(all(unlist(res$list_u) == thr))
})

test_that("get_spatiotemp_excess removes zeros correctly", {
  set.seed(42)
  data <- matrix(runif(100, 0, 10), nrow = 10)
  data[sample(length(data), 20)] <- 0
  colnames(data) <- paste0("S", 1:10)

  res_no_rm <- get_spatiotemp_excess(data, quantile = 0.9, remove_zeros = FALSE)
  res_rm <- get_spatiotemp_excess(data, quantile = 0.9, remove_zeros = TRUE)

  expect_true(length(res_no_rm$list_s) >= length(res_rm$list_s))
})

test_that("get_spatiotemp_excess errors properly on invalid args", {
  data <- matrix(runif(50), nrow = 10, ncol = 5)
  colnames(data) <- paste0("S", 1:5)

  expect_error(get_spatiotemp_excess(data, quantile = 0.9, threshold = 5))
  expect_error(get_spatiotemp_excess(data))
})


# -----------------------------------------------------------------------
# get_s0t0_pairs()
# -----------------------------------------------------------------------

test_that("get_s0t0_pairs selects valid space-time pairs", {
  set.seed(123)
  data <- matrix(runif(50, 0, 10), nrow = 10, ncol = 5)
  colnames(data) <- paste0("P", 1:5)

  sites_coords <- data.frame(
    Longitude = c(0, 1, 2, 3, 4),
    Latitude = c(0, 0, 0, 0, 0)
  )
  rownames(sites_coords) <- paste0("P", 1:5)

  set_st_excess <- get_spatiotemp_excess(data, quantile = 0.9)
  res <- get_s0t0_pairs(sites_coords, data,
                        min_spatial_dist = 1,
                        episode_size = 2,
                        set_st_excess = set_st_excess,
                        n_max_episodes = 100,
                        latlon = FALSE)

  expect_s3_class(res, "data.table")
  expect_true(all(c("s0", "t0", "u_s0") %in% names(res)))
  expect_true(nrow(res) > 0)
  # if same t0, s0 should be different
  expect_true(all(!duplicated(res[, .(s0, t0)])))
  # if same s0, t0 should be different
  expect_true(all(!duplicated(res[, .(t0, s0)])))
  # if same t0, distances between s0 and other sites should be >= min_spatial_dist
  for (i in 1:nrow(res)) {
    s0 <- res$s0[i]
    t0 <- res$t0[i]
    other_sites <- setdiff(rownames(sites_coords), s0)
    dists <- sqrt((sites_coords[other_sites, "Longitude"] - sites_coords[s0, "Longitude"])^2 +
                    (sites_coords[other_sites, "Latitude"] - sites_coords[s0, "Latitude"])^2)
    expect_true(all(dists >= 1))
  }
})

test_that("get_s0t0_pairs respects both spatial and temporal constraints", {
    set.seed(42)
    data <- matrix(runif(50, 0, 10), nrow = 10, ncol = 5)
    colnames(data) <- paste0("P", 1:5)
    sites_coords <- data.frame(Longitude = 1:5, Latitude = rep(0, 5))
    rownames(sites_coords) <- paste0("P", 1:5)

    set_st_excess <- get_spatiotemp_excess(data, quantile = 0.9)

    # Case 1: very strict constraint (large distance + large episode)
    res_strict <- get_s0t0_pairs(
        sites_coords, data,
        min_spatial_dist = 10, episode_size = 10,
        set_st_excess = set_st_excess, latlon = FALSE
    )
    expect_true(nrow(res_strict) <= 1)

    # Case 2: more relaxed temporal constraint
    res_relaxed <- get_s0t0_pairs(
        sites_coords, data,
        min_spatial_dist = 10, episode_size = 5,
        set_st_excess = set_st_excess, latlon = FALSE
    )
    expect_true(nrow(res_relaxed) >= nrow(res_strict))

    # Verify that each pair respects the constraints
    if (nrow(res_relaxed) > 1) {
        for (i in 2:nrow(res_relaxed)) {
                prev <- res_relaxed[1:(i - 1)]
                curr <- res_relaxed[i]
                spatial_distances <- get_dist_mat(sites_coords, latlon = FALSE)[curr$s0, prev$s0]
                temporal_distances <- abs(curr$t0 - prev$t0)
                expect_true(all(spatial_distances >= 10 | temporal_distances >= 5))
        }
    }

})

test_that("get_s0t0_pairs enforces spatial and temporal exclusion", {
  set.seed(456)
  data <- matrix(runif(50, 0, 10), nrow = 10, ncol = 5)
  colnames(data) <- paste0("P", 1:5)

  sites_coords <- data.frame(Longitude = 1:5, Latitude = rep(0, 5))
  rownames(sites_coords) <- paste0("P", 1:5)
  dist_mat <- as.matrix(dist(sites_coords[, c("Longitude", "Latitude")]))
  set_st_excess <- get_spatiotemp_excess(data, quantile = 0.5)
  res <- get_s0t0_pairs(sites_coords, data,
                        min_spatial_dist = 10,  # very large -> almost no pairs
                        episode_size = 10,
                        set_st_excess = set_st_excess,
                        latlon = FALSE)
  expect_true(nrow(res) <= 1)
})

test_that("get_s0t0_pairs detects missing coordinates", {
  data <- matrix(runif(20), nrow = 10, ncol = 2)
  colnames(data) <- c("A", "B")

  coords <- data.frame(Longitude = 1, Latitude = 2)
  rownames(coords) <- "A"

  st_excess <- get_spatiotemp_excess(data, quantile = 0.9)
  expect_error(get_s0t0_pairs(coords, data, 1, 2, st_excess))
})


# -----------------------------------------------------------------------
# Integration test for get_spatiotemp_excess and get_s0t0_pairs
# -----------------------------------------------------------------------

test_that("get_spatiotemp_excess + get_s0t0_pairs integration works on complex data", {
    set.seed(999)
    n_sites <- 5
    n_time <- 50
    site_names <- paste0("S", 1:n_sites)
    coords <- data.frame(
        Longitude = seq(0, 20, length.out = n_sites),
        Latitude = rep(0, n_sites)
    )
    rownames(coords) <- site_names

    # Simulated data with 2 extreme episodes
    data <- matrix(rnorm(n_sites * n_time, 10, 4), nrow = n_time, ncol = n_sites)
    colnames(data) <- site_names
    data[25, 3] <- 40
    data[26, 4] <- 35

    # Calculate spatiotemporal excess
    st_excess <- get_spatiotemp_excess(data, quantile = 0.95, remove_zeros = TRUE)
    expect_true(length(st_excess$list_s) > 0)

    # Select space-time pairs
    pairs <- get_s0t0_pairs(coords, data,
                                                    min_spatial_dist = 2,
                                                    episode_size = 5,
                                                    set_st_excess = st_excess,
                                                    n_max_episodes = 50,
                                                    latlon = FALSE)

    expect_true(nrow(pairs) > 0)
    # For each site present in pairs$s0, compare to local quantile
    for (site in unique(pairs$s0)) {
        q_site <- quantile(data[, site], 0.95, na.rm = TRUE)
        u_vals <- pairs$u_s0[pairs$s0 == site]
        expect_true(all(u_vals >= q_site - 1e-8),
                                info = paste("Threshold mismatch for site", site))
    }
})




