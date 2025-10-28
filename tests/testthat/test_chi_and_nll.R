# ======================================================================
# Tests for theoretical_chi, neg_ll, neg_ll_composite, neg_ll_composite_fixed_eta
# ======================================================================

context("theoretical chi and likelihoods")

# --- Small helper to build a minimal df_lags with columns required
make_df_lags <- function(tau_vec = c(0, 1, 2),
                         s1x = 0, s1y = 0, s2x = 1, s2y = 0) {
  n <- length(tau_vec)
  df <- data.frame(
    s1   = rep(1L, n),
    s2   = rep(2L, n),
    tau  = tau_vec,
    s1x  = rep(s1x, n),
    s1y  = rep(s1y, n),
    s2x  = rep(s2x, n),
    s2y  = rep(s2y, n),
    hnorm = rep(sqrt((s2x - s1x)^2 + (s2y - s1y)^2), n)
  )
  df
}

# --- A small "excesses" factory
make_excesses <- function(n, kij_vals = NULL, Tobs_vals = NULL) {
  if (is.null(kij_vals)) kij_vals <- rep(0L, n)
  if (is.null(Tobs_vals)) Tobs_vals <- rep(1L, n)
  data.frame(kij = as.integer(kij_vals), Tobs = as.integer(Tobs_vals))
}

# ======================================================================
# theoretical_chi
# ======================================================================

test_that("theoretical_chi basic structure and bounds", {
  set.seed(1)
  df_lags <- make_df_lags(tau_vec = c(0, 1, 2))
  # params: beta1, beta2, alpha1, alpha2, advx, advy
  params <- c(1, 1, 1, 1, 0, 0)

  res <- theoretical_chi(params, df_lags, distance = "euclidean", normalize = FALSE)

  # structure
  expect_true(is.data.frame(res))
  expect_true(all(c("chi", "vario") %in% names(res)))
  expect_equal(nrow(res), nrow(df_lags))

  # chi in (0,1)
  expect_true(all(res$chi > 0 & res$chi < 1))
  # numeric finite
  expect_true(all(is.finite(res$chi)))
  expect_true(all(is.finite(res$vario)))
})

test_that("theoretical_chi: zero-variogram parameters give chi ~ 1", {
  df_lags <- make_df_lags(tau_vec = c(0, 1, 5))
  # beta1=beta2=0 => vario = 0, chi = 2*(1 - Phi(0)) = 1, then clipped to < 1
  params <- c(0, 0, 1, 1, 0, 0)

  res <- theoretical_chi(params, df_lags, distance = "euclidean", normalize = FALSE)

  expect_true(all(res$vario == 0))
  expect_true(all(res$chi < 1 & res$chi > 1e-12))
  # all chi are equal and very close to 1
  expect_lt(max(abs(res$chi - res$chi[1])), 1e-15)
})

test_that("theoretical_chi: normalize TRUE vs FALSE yield identical chi", {
  # multiple distinct lags so medians are meaningful
  df_lags <- rbind(
    make_df_lags(tau_vec = c(0, 2, 4), s2x = 2, s2y = 0),
    make_df_lags(tau_vec = c(1, 3, 5), s2x = 4, s2y = 3)
  )
  params <- c(0.8, 0.5, 1.2, 0.7, 0.0, 0.0)

  res_no  <- theoretical_chi(params, df_lags, distance = "euclidean", normalize = FALSE)
  res_yes <- theoretical_chi(params, df_lags, distance = "euclidean", normalize = TRUE)

  # Normalization rescales lags and betas consistently in the function:
  # chi should match (up to tight numerical tolerance)
  expect_equal(res_yes$chi, res_no$chi, tolerance = 1e-10)
})

test_that("theoretical_chi: lalpha branch matches manual vario for a simple row", {
  df_lags <- make_df_lags(tau_vec = 3) # one row
  params  <- c(2, 3, 1.5, 0.5, 0, 0)   # beta1,beta2,alpha1,alpha2, no advection

  res <- theoretical_chi(params, df_lags, distance = "lalpha", normalize = FALSE)

  # Manual vario: (2*beta1)*(|hx|^alpha1 + |hy|^alpha1) + (2*beta2)*|t|^alpha2
  hx <- df_lags$s2x - df_lags$s1x
  hy <- df_lags$s2y - df_lags$s1y
  tau <- abs(df_lags$tau)

  vario_manual <- (2*params[1])*(abs(hx)^params[3] + abs(hy)^params[3]) +
                  (2*params[2])*(abs(tau)^params[4])

  expect_equal(res$vario[1], vario_manual[1], tolerance = 1e-12)
})

test_that("theoretical_chi: advection changes the advected distance", {
  # If advection is non-zero and tau != 0, advected coordinates shift s2.
  df_lags <- make_df_lags(tau_vec = c(0, 2))  # two rows, tau 0 and 2
  params0 <- c(1, 1, 1, 1, 0, 0)  # no advection
  paramsA <- c(1, 1, 1, 1, 0.5, 0) # advect along +x, so s2xv = s2x - 0.5 * tau

  res0 <- theoretical_chi(params0, df_lags, distance = "euclidean", normalize = FALSE)
  resA <- theoretical_chi(paramsA, df_lags, distance = "euclidean", normalize = FALSE)

  # For tau=0, identical; for tau=2, different (unless by coincidence)
  expect_equal(res0$chi[1], resA$chi[1], tolerance = 1e-15)
  expect_false(isTRUE(all.equal(res0$chi[2], resA$chi[2], tolerance = 1e-12)))
})

# ======================================================================
# neg_ll
# ======================================================================

test_that("neg_ll returns finite numeric and honors hmax filtering", {
  df_lags <- rbind(
    make_df_lags(tau_vec = 0),               # hnorm=1
    make_df_lags(tau_vec = 1, s2x = 3)       # hnorm=3
  )
  params <- c(0.4, 0.3, 1.0, 1.2, 0, 0)

  # Tobs = 1 for each, kij mix
  excesses <- make_excesses(nrow(df_lags), kij_vals = c(0, 1), Tobs_vals = c(1, 1))

  nll_all <- neg_ll(params, df_lags, excesses, hmax = NA, distance="euclidean", normalize = FALSE)
  nll_h1  <- neg_ll(params, df_lags, excesses, hmax = 1.5, distance="euclidean", normalize = FALSE)

  expect_true(is.numeric(nll_all) && length(nll_all) == 1 && is.finite(nll_all))
  expect_true(is.numeric(nll_h1)  && length(nll_h1)  == 1 && is.finite(nll_h1))

  # With hmax small, only the first row (hnorm=1) is kept => different value
  expect_false(isTRUE(all.equal(nll_all, nll_h1, tolerance = 1e-12)))
})

test_that("neg_ll matches closed-form when vario=0 (chi ~ 1)", {
  # beta1=beta2=0 => chi ~= 1 - 1e-12 (due to clipping). p=1 (rpar=TRUE).
  df_lags <- make_df_lags(tau_vec = c(0, 1, 2))
  params  <- c(0, 0, 1, 1, 0, 0)
  excesses <- make_excesses(nrow(df_lags),
                            kij_vals  = c(0, 1, 0),
                            Tobs_vals = c(1, 1, 1))

  # Compute expected analytically with the same eps as in the code
  eps   <- 1e-12
  chi   <- rep(1 - eps, nrow(df_lags)) # from clipping
  pchi  <- pmax(pmin(1 - 1 * chi, 1 - eps), eps)  # => eps
  ll    <- c( log(pchi[1]), log(chi[2]), log(pchi[3]) )
  nll_expected <- -sum(ll)

  nll_got <- neg_ll(params, df_lags, excesses, rpar = TRUE)
  expect_equal(nll_got, nll_expected, tolerance = 1e-10)
})

# ======================================================================
# neg_ll_composite
# ======================================================================

test_that("neg_ll_composite equals the sum of per-episode neg_ll", {
  # Build two tiny episodes that share the same lags structure
  lags1 <- make_df_lags(tau_vec = c(0, 1))
  lags2 <- make_df_lags(tau_vec = c(1, 2), s2x = 2)

  exc1 <- make_excesses(nrow(lags1), kij_vals = c(1, 0), Tobs_vals = c(1, 1))
  exc2 <- make_excesses(nrow(lags2), kij_vals = c(0, 1), Tobs_vals = c(1, 1))

  params <- c(0.6, 0.8, 1.0, 1.0)  # length 4 => adv defaults to (0,0) in composite

  ll1 <- neg_ll(c(params, 0, 0), lags1, exc1)
  ll2 <- neg_ll(c(params, 0, 0), lags2, exc2)

  total_ref <- ll1 + ll2

  nll_comp <- neg_ll_composite(params, list_episodes = list(NA, NA),
                               list_excesses = list(exc1, exc2),
                               list_lags = list(lags1, lags2),
                               wind_df = NA, distance="euclidean",
                               rpar = TRUE, normalize = FALSE)

  expect_equal(nll_comp, total_ref, tolerance = 1e-12)
})

test_that("neg_ll_composite: wind advection changes the result when tau != 0", {
  # Lags with non-zero tau so advection matters
  lags1 <- make_df_lags(tau_vec = c(0, 2), s2x = 3, s2y = 0)
  lags2 <- make_df_lags(tau_vec = c(1, 3), s2x = 3, s2y = 0)

  exc1 <- make_excesses(nrow(lags1), kij_vals = c(1, 0))
  exc2 <- make_excesses(nrow(lags2), kij_vals = c(0, 1))

  # wind_df provides vx, vy per episode; eta1, eta2 set so adv is non-zero
  wind_df <- data.frame(vx = c(5, -3), vy = c(0, 0)) # km/h or consistent units
  params_wind <- c(0.6, 0.8, 1.0, 1.0, 0.1, 1.0)     # eta1=0.1, eta2=1 => adv = 0.1*sign(v)*|v|

  nll_nowind <- neg_ll_composite(c(0.6, 0.8, 1.0, 1.0, 0, 0),
                                 list_episodes = list(NA, NA),
                                 list_excesses = list(exc1, exc2),
                                 list_lags = list(lags1, lags2),
                                 wind_df = NA,
                                 distance="euclidean", rpar = TRUE)

  nll_wind <- neg_ll_composite(params_wind,
                               list_episodes = list(NA, NA),
                               list_excesses = list(exc1, exc2),
                               list_lags = list(lags1, lags2),
                               wind_df = wind_df,
                               distance="euclidean", rpar = TRUE)

  # Expect a difference due to advected coordinates
  expect_false(isTRUE(all.equal(nll_nowind, nll_wind, tolerance = 1e-12)))
})

# ======================================================================
# neg_ll_composite_fixed_eta
# ======================================================================

test_that("neg_ll_composite_fixed_eta reproduces neg_ll_composite when fixing eta", {
  lags <- make_df_lags(tau_vec = c(0, 1, 2), s2x = 2)
  exc  <- make_excesses(nrow(lags), kij_vals = c(1, 0, 1))
  wind_df <- data.frame(vx = 2, vy = -1)

  # Base params: (beta1,beta2,alpha1,alpha2,eta1,eta2)
  base <- c(0.5, 0.7, 1.1, 0.9, 0.2, 1.3)

  nll_free <- neg_ll_composite(base,
                               list_episodes = list(NA),
                               list_excesses = list(exc),
                               list_lags = list(lags),
                               wind_df = wind_df,
                               rpar = TRUE, distance="euclidean")

  # Lock eta1, eta2 to the same values
  nll_fixed <- neg_ll_composite_fixed_eta(params = base[1:4],
                                          list_episodes = list(NA),
                                          list_excesses = list(exc),
                                          list_lags = list(lags),
                                          wind_df = wind_df,
                                          fixed_eta1 = base[5],
                                          fixed_eta2 = base[6],
                                          rpar = TRUE, distance="euclidean")

  expect_equal(nll_free, nll_fixed, tolerance = 1e-12)
})

test_that("neg_ll_composite validates distance argument", {
  lags <- make_df_lags(tau_vec = 0)
  exc  <- make_excesses(1, kij_vals = 1)
  expect_error(
    neg_ll_composite(c(0.5,0.5,1,1,0,0),
                     list_episodes = list(NA),
                     list_excesses = list(exc),
                     list_lags = list(lags),
                     distance = "invalid"),
    "should be 'lalpha' or 'euclidean'"
  )
})
