test_that("sim_episode returns correct structure and reasonable values", {
  # --- Input setup -------------------------------------------------------------
  coords <- expand.grid(x = 1:2, y = 1:2)   # 4 sites
  times <- 0:2                            # 3 time steps
  params_vario <- list(beta1 = 1, beta2 = 1, alpha1 = 1, alpha2 = 1)
  
  # Marginal parameters for each site (vectors of length = nrow(coords))
  params_margins <- list(
    xi    = rep(0.1, 4),
    sigma = rep(1, 4),
    kappa = rep(0.5, 4)
  )
  
  # Local thresholds
  u_s <- c(10, 12, 15, 20)
  
  # Call the function
  X <- sim_episode(
    params_vario  = params_vario,
    params_margins = params_margins,
    coords = coords,
    times  = times,
    adv    = c(0, 0),
    t0     = 1,
    s0     = c(1, 1),
    u_s    = u_s
  )
  
  # --- Structural checks ------------------------------------------------------
  expect_true(is.array(X))  # X should be an array
  expect_identical(dim(X), c(2L, 2L, 3L))  # Dimensions should match (nx, ny, nt)
  
  # --- Value sanity checks ----------------------------------------------------
  expect_true(all(is.finite(X)))                # No NA or Inf values
  expect_true(all(X >= 0))                      # Rainfall must be non-negative

})
