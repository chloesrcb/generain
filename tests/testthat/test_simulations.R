
test_that("sim_rpareto returns expected structure with random s0", {
  skip_if_not_installed("RandomFields")
  skip_if_not_installed("evd")

  ngrid <- 3
  spa <- 1:ngrid
  temp <- 0:2
  nres <- 2
  param <- c(0.4, 0.2, 1.5, 1)
  adv <- c(0.1, 0.2)
  s0 <- c(2, 2)
  t0 <- 0
  random_s0 <- TRUE
  s0_radius <- 1.5

  result <- sim_rpareto(param[1], param[2], param[3], param[4],
                        spa, spa, temp,
                        adv = adv, t0 = t0, nres = nres,
                        random_s0 = random_s0, s0 = s0,
                        s0_radius = s0_radius, seed = 123)

  expect_type(result, "list")
  expect_true("Z" %in% names(result))
  expect_true("s0_used" %in% names(result))

  expect_equal(dim(result$Z), c(ngrid, ngrid, length(temp), nres, 1))
  expect_equal(length(result$s0_used), nres)

  for (i in seq_len(nres)) {
    expect_equal(length(result$s0_used[[i]]), 1)
    s0_curr <- result$s0_used[[i]][[1]]
    expect_equal(names(s0_curr), c("x", "y"))
    s0_vec <- as.numeric(s0_curr[1, c("x", "y")])
    expect_true(all(s0_vec %in% 1:ngrid))
    expect_true(all(s0_vec == as.integer(s0_vec)))

    dist <- sqrt((s0_curr$x - s0[1])^2 + (s0_curr$y - s0[2])^2)
    expect_lte(dist, s0_radius)
  }
})

test_that("sim_rpareto supports multiple fixed s0 points", {
  skip_if_not_installed("RandomFields")
  skip_if_not_installed("evd")

  ngrid <- 3
  spa <- 1:ngrid
  temp <- 0:2
  nres <- 2
  param <- c(0.4, 0.2, 1.5, 1)
  adv <- c(0, 0)
  s0 <- rbind(c(1, 1), c(3, 2))
  t0 <- 0

  result <- sim_rpareto(param[1], param[2], param[3], param[4],
                        spa, spa, temp,
                        adv = adv, t0 = t0, nres = nres,
                        random_s0 = FALSE, s0 = s0, seed = 456)

  expect_equal(dim(result$Z), c(ngrid, ngrid, length(temp), nres, 2))

  for (i in seq_len(nres)) {
    expect_equal(length(result$s0_used[[i]]), 2)
    for (j in seq_len(nrow(s0))) {
      expect_equal(as.numeric(result$s0_used[[i]][[j]][1, ]), s0[j, ])
    }
  }
})

test_that("sim_rpareto returns a value at s0 >= 1 for t0", {
  skip_if_not_installed("RandomFields")
  skip_if_not_installed("evd")

  ngrid <- 3
  spa <- 1:ngrid
  temp <- 0:2
  nres <- 1
  param <- c(0.4, 0.2, 1.5, 1)
  adv <- c(0, 0)
  s0 <- c(1, 3)
  t0 <- 0

  result <- sim_rpareto(param[1], param[2], param[3], param[4],
                        x = spa, y = spa, t = temp,
                        adv = adv, t0 = t0, nres = nres,
                        random_s0 = FALSE, s0 = s0, seed = 999)

  t0_idx <- which(temp == t0)
  val_at_s0 <- result$Z[s0[1], s0[2], t0_idx, 1, 1]
  expect_gt(val_at_s0, 1 - 1e-12)
})




test_that("compute_st_variogram returns gamma_s + gamma_t", {
  # Define a simple 1x1x1 grid (only one point in space and time)
  grid <- data.frame(
    x = 1, y = 1, t = 0,
    shifted_x = 1, shifted_y = 1
  )

  # Dummy variograms
  gamma_space <- list(2)   # gamma_s = 2
  gamma_temp <- list(3)    # gamma_t = 3
  adv <- c(0, 0)           # No advection

  # Run the function
  result <- compute_st_variogram(
    grid = grid,
    gamma_space = gamma_space,
    gamma_temp = gamma_temp,
    adv = adv
  )

  # Expect gamma_s + gamma_t = 2 + 3 = 5
  expect_equal(result[1, 1, 1], 5)
})


test_that("compute_st_variogram returns gamma_s + gamma_t", {
  # Define a simple 1x1x1 grid (only one point in space and time)
  grid <- data.frame(
    x = 1, y = 1, t = 0,
    shifted_x = 1, shifted_y = 1
  )

  # Dummy variograms
  gamma_space_x <- list(2)   # gamma_sx = 2
  gamma_space_y <- list(1)   # gamma_sy = 1
  gamma_temp <- list(3)    # gamma_t = 3
  adv <- c(0, 0)           # No advection

  # Run the function
  result <- compute_st_variogram(
    grid = grid,
    gamma_space_x = gamma_space_x,
    gamma_space_y = gamma_space_y,
    gamma_temp = gamma_temp,
    adv = adv
  )

  # Expect gamma_s + gamma_t = 2 + 1 + 3 = 6
  expect_equal(result[1, 1, 1], 6)
})

test_that("compute_st_variogram handles multi-point grid with zero advection", {
  grid <- expand.grid(x = 1:2, y = 1:2, t = 0:1)
  grid$shifted_x <- grid$x
  grid$shifted_y <- grid$y

  coords <- unique(grid[, c("shifted_x", "shifted_y")])
  nsites <- nrow(coords)

  gamma_space <- as.list(seq_len(nsites))
  gamma_temp <- list(10, 20) # t=0 -> 10, t=1 -> 20

  result <- compute_st_variogram(
    grid = grid,
    gamma_space = gamma_space,
    gamma_temp = gamma_temp,
    adv = c(0, 0)
  )

  expected <- array(NA, dim = c(2, 2, 2))
  for (i in seq_len(nrow(grid))) {
    t_idx <- grid$t[i] + 1
    ind_g_s <- if (i %% nsites == 0) {
      i - nsites * (i %/% nsites - 1)
    } else {
      i - nsites * (i %/% nsites)
    }
    expected[grid$x[i], grid$y[i], t_idx] <- gamma_space[[ind_g_s]] + gamma_temp[[t_idx]]
  }

  expect_equal(result, expected)
})

test_that("compute_st_variogram handles lalpha with advection", {
  grid <- expand.grid(x = 1:2, y = 1:2, t = 0:1)
  grid$shifted_x <- grid$x - grid$t * 0.5
  grid$shifted_y <- grid$y - grid$t * 0.25

  nsites <- nrow(grid)
  gamma_space_x <- as.list(seq_len(nsites))
  gamma_space_y <- as.list(seq_len(nsites) * 2)
  gamma_temp <- list(1, 3)

  result <- compute_st_variogram(
    grid = grid,
    gamma_space_x = gamma_space_x,
    gamma_space_y = gamma_space_y,
    gamma_temp = gamma_temp,
    adv = c(0.5, 0.25)
  )

  expected <- array(NA, dim = c(2, 2, 2))
  for (i in seq_len(nrow(grid))) {
    t_idx <- grid$t[i] + 1
    expected[grid$x[i], grid$y[i], t_idx] <-
      gamma_space_x[[i]] + gamma_space_y[[i]] + gamma_temp[[t_idx]]
  }

  expect_equal(result, expected)
})

test_that("compute_st_variogram errors when spatial variogram inputs are missing", {
  grid <- data.frame(
    x = 1, y = 1, t = 0,
    shifted_x = 1, shifted_y = 1
  )

  expect_error(
    compute_st_variogram(grid = grid, gamma_temp = list(0), adv = c(0, 0)),
    "Provide either gamma_space"
  )
})

test_that("compute_st_gaussian_process handles euclidean with zero advection", {
  grid <- expand.grid(x = 1:2, y = 1:2, t = 0:1)
  grid$shifted_x <- grid$x
  grid$shifted_y <- grid$y

  coords <- cbind(grid$shifted_x, grid$shifted_y)
  key <- paste(coords[, 1], coords[, 2], sep = "_")
  key_u <- unique(key)
  map_idx <- match(key, key_u)

  W_s <- 1:4
  W_t <- c(10, 20)

  result <- compute_st_gaussian_process(
    grid = grid,
    W_s = W_s,
    W_t = W_t,
    adv = c(0, 0)
  )

  expected <- array(NA, dim = c(2, 2, 2))
  t_levels <- sort(unique(grid$t))
  t_index <- match(grid$t, t_levels)

  for (i in seq_len(nrow(grid))) {
    expected[grid$x[i], grid$y[i], t_index[i]] <- W_s[map_idx[i]] + W_t[t_index[i]]
  }

  expect_equal(result, expected)
})

test_that("compute_st_gaussian_process handles lalpha with advection", {
  grid <- expand.grid(x = 1:2, y = 1:2, t = 0:1)
  grid$shifted_x <- grid$x - grid$t * 0.5
  grid$shifted_y <- grid$y - grid$t * 0.25

  W_sx <- seq_len(nrow(grid))
  W_sy <- 100 + seq_len(nrow(grid))
  W_t <- c(5, 6)

  result <- compute_st_gaussian_process(
    grid = grid,
    W_s_x = W_sx,
    W_s_y = W_sy,
    W_t = W_t,
    adv = c(0.5, 0.25)
  )

  expected <- array(NA, dim = c(2, 2, 2))
  t_levels <- sort(unique(grid$t))
  t_index <- match(grid$t, t_levels)

  for (i in seq_len(nrow(grid))) {
    expected[grid$x[i], grid$y[i], t_index[i]] <- W_sx[i] + W_sy[i] + W_t[t_index[i]]
  }

  expect_equal(result, expected)
})

test_that("compute_st_gaussian_process errors when spatial inputs are missing", {
  grid <- data.frame(
    x = 1, y = 1, t = 0,
    shifted_x = 1, shifted_y = 1
  )

  expect_error(
    compute_st_gaussian_process(grid = grid, W_t = 1, adv = c(0, 0)),
    "Provide either W_s"
  )
})

test_that("compute_st_variogram handles advection correctly", {
  # Define a simple grid with two time points
  grid <- data.frame(
    x = c(1, 1), y = c(1, 1), t = c(0, 1),
    shifted_x = c(1, 1), shifted_y = c(1, 1)
  )

  # Dummy variograms
  gamma_space <- list(2, 3)   # gamma_s = 2
  gamma_temp <- list(3, 2)    # gamma_t = 3
  adv <- c(1, 0)           # Advection in x-direction

  # Run the function
  result <- compute_st_variogram(
    grid = grid,
    gamma_space = gamma_space,
    gamma_temp = gamma_temp,
    adv = adv
  )

  # For t=0, no advection effect, expect gamma_s + gamma_t = 2 + 3 = 5
  expect_equal(result[1, 1, 1], 5)

  # For t=1, advection shifts point to (2,1), expect same variogram value
  expect_equal(result[1, 1, 2], 5)
})