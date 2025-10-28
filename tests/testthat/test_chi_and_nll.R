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



context("Advanced property-based tests for chi and likelihoods")

# ----------------- Helpers -----------------
make_df_lags <- function(tau_vec, s1x=0, s1y=0, s2x=1, s2y=0) {
  n <- length(tau_vec)
  data.frame(
    s1   = rep(1L, n),
    s2   = rep(2L, n),
    tau  = tau_vec,
    s1x  = rep(s1x, n),
    s1y  = rep(s1y, n),
    s2x  = rep(s2x, n),
    s2y  = rep(s2y, n),
    hnorm = rep(sqrt((s2x - s1x)^2 + (s2y - s1y)^2), n)
  )
}
make_excesses <- function(n, kij = NULL, Tobs = NULL) {
  if (is.null(kij)) kij <- rep(0L, n)
  if (is.null(Tobs)) Tobs <- rep(1L, n)
  data.frame(kij = as.integer(kij), Tobs = as.integer(Tobs))
}

# ================== theoretical_chi ==================

test_that("chi decreases when variogram scales (beta1,beta2) increase", {
  # Deux lignes: même h, deux valeurs de tau distinctes
  df_lags <- make_df_lags(tau_vec = c(0, 3), s2x = 4, s2y = 0)
  # beta1,beta2 petit -> vario plus petit -> chi plus grand
  p_small <- c(0.2, 0.2, 1.0, 1.0, 0, 0)
  p_large <- c(2.0, 2.0, 1.0, 1.0, 0, 0)

  chi_small <- theoretical_chi(p_small, df_lags, distance="euclidean", normalize=FALSE)$chi
  chi_large <- theoretical_chi(p_large, df_lags, distance="euclidean", normalize=FALSE)$chi

  expect_true(all(chi_large < chi_small))
})

test_that("chi is identical for normalize TRUE/FALSE under the intended scaling", {
  # lags hétérogènes pour que les médianes soient bien définies
  df_lags <- rbind(
    make_df_lags(tau_vec = c(0,2,4), s2x=2, s2y=1),
    make_df_lags(tau_vec = c(1,3,5), s2x=5, s2y=4)
  )
  params <- c(0.7, 0.9, 1.3, 0.6, 0.0, 0.0)

  res_no  <- theoretical_chi(params, df_lags, distance="euclidean", normalize=FALSE)
  res_yes <- theoretical_chi(params, df_lags, distance="euclidean", normalize=TRUE)

  expect_equal(res_yes$chi, res_no$chi, tolerance = 1e-10)
})

test_that("advection is neutral when tau==0 and changes chi when tau!=0", {
  df_lags <- make_df_lags(tau_vec = c(0, 2), s2x=3, s2y=0)
  params0 <- c(0.8, 0.8, 1.0, 1.0, 0.0, 0.0)
  paramsA <- c(0.8, 0.8, 1.0, 1.0, 0.5, 0.0)  # advx=0.5

  chi0 <- theoretical_chi(params0, df_lags)$chi
  chiA <- theoretical_chi(paramsA, df_lags)$chi

  expect_equal(chi0[1], chiA[1], tolerance = 1e-15)  # tau=0
  expect_false(isTRUE(all.equal(chi0[2], chiA[2], tolerance = 1e-12))) # tau=2
})

test_that("lalpha branch equals manual formula on random rows", {
  set.seed(123)
  df_lags <- data.frame(
    s1=1, s2=2,
    tau = sample(1:5, 5, replace=TRUE),
    s1x=0, s1y=0,
    s2x = runif(5, -2, 2),
    s2y = runif(5, -2, 2),
    hnorm = NA_real_
  )
  df_lags$hnorm <- sqrt((df_lags$s2x - df_lags$s1x)^2 + (df_lags$s2y - df_lags$s1y)^2)
  params <- c(1.7, 0.9, 1.4, 0.8, 0, 0)

  out <- theoretical_chi(params, df_lags, distance="lalpha", normalize=FALSE)
  hx <- df_lags$s2x - df_lags$s1x
  hy <- df_lags$s2y - df_lags$s1y
  tt <- abs(df_lags$tau)

  vario_manual <- (2*params[1])*(abs(hx)^params[3] + abs(hy)^params[3]) +
                  (2*params[2])*(abs(tt)^params[4])
  expect_equal(out$vario, vario_manual, tolerance=1e-12)
})

test_that("extreme lags produce chi within clipping bounds and finite", {
  # Très grands h et tau
  df_lags <- data.frame(
    s1=1, s2=2, tau = c(0, 50, 200),
    s1x=0, s1y=0, s2x= c(0, 1e3, 1e4), s2y=0,
    hnorm = c(0, 1e3, 1e4)
  )
  params <- c(1, 1, 1.2, 0.7, 0, 0)
  out <- theoretical_chi(params, df_lags)

  expect_true(all(is.finite(out$chi)))
  expect_true(all(out$chi >= 1e-12 & out$chi <= 1 - 1e-12))
})

# ================== neg_ll ==================

test_that("neg_ll is additive over rows and matches manual per-row sum", {
  df_lags <- rbind(
    make_df_lags(tau_vec = 0, s2x=1),
    make_df_lags(tau_vec = 2, s2x=3)
  )
  params <- c(0.5, 0.4, 1.0, 1.0, 0, 0)
  exc <- data.frame(kij=c(1L, 0L), Tobs=c(1L, 2L))

  # calcul interne
  nll <- neg_ll(params, df_lags, exc, rpar=TRUE)

  # recompose manuellement: chi, pchi, ll
  chi <- theoretical_chi(params, df_lags)$chi
  p <- 1
  eps <- 1e-12
  pchi <- pmax(pmin(1 - p*chi, 1 - eps), eps)
  ll_row1 <- exc$kij[1]*log(chi[1]) + (exc$Tobs[1]-exc$kij[1])*log(pchi[1])
  ll_row2 <- exc$kij[2]*log(chi[2]) + (exc$Tobs[2]-exc$kij[2])*log(pchi[2])
  nll_manual <- -(ll_row1 + ll_row2)

  expect_equal(nll, nll_manual, tolerance = 1e-12)
})

test_that("neg_ll increases with beta when all kij=1 (monotonicity check)", {
  # Si tous kij=1, augmenter vario (via beta) diminue chi et augmente -log-lik
  df_lags <- rbind(
    make_df_lags(tau_vec=0, s2x=2),
    make_df_lags(tau_vec=3, s2x=2)
  )
  exc <- make_excesses(nrow(df_lags), kij = c(1,1), Tobs = c(1,1))

  nll1 <- neg_ll(c(0.3, 0.3, 1.0, 1.0, 0,0), df_lags, exc)
  nll2 <- neg_ll(c(1.5, 1.5, 1.0, 1.0, 0,0), df_lags, exc)

  expect_gt(nll2, nll1)
})

test_that("hmax subsetting equals recomputing on the filtered subset", {
  df_lags <- rbind(
    make_df_lags(tau_vec=0, s2x=1),  # h=1
    make_df_lags(tau_vec=1, s2x=3)   # h=3
  )
  params <- c(0.7, 0.5, 1.0, 1.0, 0,0)
  exc <- make_excesses(2, kij=c(1,0), Tobs=c(2,2))

  # nll avec hmax=1.5
  nll_hmax <- neg_ll(params, df_lags, exc, hmax=1.5)

  # filtre manuel puis nll
  keep <- df_lags$hnorm <= 1.5
  nll_ref <- neg_ll(params, df_lags[keep,], exc[keep,])

  expect_equal(nll_hmax, nll_ref, tolerance = 1e-12)
})

# ================== neg_ll_composite ==================

test_that("composite equals sum of per-episode nll with advection per-episode", {
  # 3 épisodes avec adv différents (via wind_df)
  l1 <- make_df_lags(tau_vec=c(0,2), s2x=2)
  l2 <- make_df_lags(tau_vec=c(1,3), s2x=2)
  l3 <- make_df_lags(tau_vec=c(2,4), s2x=2)

  e1 <- make_excesses(nrow(l1), kij=c(1,0))
  e2 <- make_excesses(nrow(l2), kij=c(0,1))
  e3 <- make_excesses(nrow(l3), kij=c(1,1))

  wind_df <- data.frame(vx=c(0, 5, -2), vy=c(0, 0, 0))  # 3 advx différents
  params  <- c(0.6, 0.8, 1.0, 1.0, 0.1, 1.0)           # eta1=0.1, eta2=1

  # somme manuelle
  adv_vec <- 0.1 * sign(wind_df$vx) * abs(wind_df$vx)^1
  nll_sum <- neg_ll(c(0.6,0.8,1.0,1.0, adv_vec[1], 0), l1, e1) +
             neg_ll(c(0.6,0.8,1.0,1.0, adv_vec[2], 0), l2, e2) +
             neg_ll(c(0.6,0.8,1.0,1.0, adv_vec[3], 0), l3, e3)

  nll_comp <- neg_ll_composite(params,
                               list_episodes=list(NA,NA,NA),
                               list_excesses=list(e1,e2,e3),
                               list_lags=list(l1,l2,l3),
                               wind_df=wind_df, rpar=TRUE)

  expect_equal(nll_comp, nll_sum, tolerance=1e-12)
})

test_that("fixed_eta wrapper reproduces composite when fixing both etas", {
  lags <- make_df_lags(tau_vec=c(0,2,4), s2x=3)
  exc  <- make_excesses(nrow(lags), kij=c(1,0,1))
  wind_df <- data.frame(vx=3, vy=-4)
  base <- c(0.5,0.7,1.1,0.9, 0.3, 1.2)

  nll_free <- neg_ll_composite(base,
                               list_episodes=list(NA),
                               list_excesses=list(exc),
                               list_lags=list(lags),
                               wind_df=wind_df, rpar=TRUE)

  nll_fixed <- neg_ll_composite_fixed_eta(params=base[1:4],
                                          list_episodes=list(NA),
                                          list_excesses=list(exc),
                                          list_lags=list(lags),
                                          wind_df=wind_df,
                                          fixed_eta1=base[5],
                                          fixed_eta2=base[6],
                                          rpar=TRUE)
  expect_equal(nll_free, nll_fixed, tolerance=1e-12)
})

# ================== Stress / edge cases ==================

test_that("handles NAs in df_lags rows gracefully (excluded in sums)", {
  df_lags <- make_df_lags(tau_vec=c(0,1))
  df_lags$hnorm[2] <- NA_real_
  params <- c(0.4, 0.6, 1.0, 1.0, 0,0)
  exc <- make_excesses(2, kij=c(1,0), Tobs=c(1,1))
  # theoretical_chi->chi has NA in vario for row2; neg_ll sums with na.rm=TRUE
  nll <- neg_ll(params, df_lags, exc)
  expect_true(is.finite(nll))
})

test_that("very small chi values are clipped -> finite nll", {
  # Gros betas et lags => chi très petit, mais pas 0 grâce au clipping
  df_lags <- make_df_lags(tau_vec=c(10, 20), s2x=100, s2y=0)
  params <- c(50, 50, 1.5, 1.5, 0,0)
  exc <- make_excesses(2, kij=c(0,0), Tobs=c(10,10))
  nll <- neg_ll(params, df_lags, exc)
  expect_true(is.finite(nll))
})





context("End-to-end optimization of neg_ll_composite")

# --- Helpers to build synthetic lags and excesses -----------------------------

make_random_df_lags <- function(n_rows, seed=1) {
  set.seed(seed)
  # s1 toujours 1, s2 aléatoire parmi {2,3,4} pour divers hx/hy
  s1  <- rep(1L, n_rows)
  s2  <- sample(2:4, n_rows, replace=TRUE)
  s1x <- rep(0, n_rows); s1y <- rep(0, n_rows)
  # positions s2 dans un petit plan euclidien
  s2x <- rnorm(n_rows, mean=3, sd=1.0)
  s2y <- rnorm(n_rows, mean=0, sd=1.0)
  tau <- sample(0:6, n_rows, replace=TRUE)
  hnorm <- sqrt((s2x - s1x)^2 + (s2y - s1y)^2)

  data.frame(s1,s2,tau,s1x,s1y,s2x,s2y,hnorm, check.names = FALSE)
}

simulate_excesses_from_params <- function(params6, df_lags, Tobs = 200L) {
  # params6 = c(beta1,beta2,alpha1,alpha2, advx, advy)
  chi <- theoretical_chi(params6, df_lags, distance="euclidean", normalize=FALSE)$chi
  # Simule kij ~ Binomial(Tobs, chi) indépendamment pour chaque ligne
  kij <- stats::rbinom(n = length(chi), size = Tobs, prob = chi)
  data.frame(kij = as.integer(kij), Tobs = as.integer(rep(Tobs, length(chi))))
}

# --- Test 1: optimisation avec vent (eta1, eta2) -----------------------------

test_that("optim on neg_ll_composite recovers parameters with wind (eta1,eta2)", {
  set.seed(42)

  # VRAIS paramètres (modérés pour éviter des chi trop extrêmes)
  true <- c(
    beta1  = 0.6,
    beta2  = 0.8,
    alpha1 = 1.2,
    alpha2 = 0.9,
    eta1   = 0.20,  # échelle advection
    eta2   = 1.10   # exposant vitesse
  )

  # 3 épisodes, vents différents pour que l’advection soit identifiable
  wind_df <- data.frame(
    vx = c( 0.0,  4.0, -2.0),
    vy = c( 0.0, -1.0,  3.0)
  )

  # On fabrique pour chaque épisode un df_lags distinct, puis les excès simulés
  list_lags <- list(
    make_random_df_lags(30, seed=101),
    make_random_df_lags(35, seed=202),
    make_random_df_lags(40, seed=303)
  )

  # Construire adv par épisode pour simuler (comme le fait neg_ll_composite)
  eta1 <- true["eta1"]; eta2 <- true["eta2"]
  adv_x <- eta1 * sign(wind_df$vx) * abs(wind_df$vx)^eta2
  adv_y <- eta1 * sign(wind_df$vy) * abs(wind_df$vy)^eta2

  # Simule les excès à partir des "vrais" paramètres épisode par épisode
  list_excesses <- vector("list", length(list_lags))
  for (i in seq_along(list_lags)) {
    params6 <- c(true[1:4], adv_x[i], adv_y[i])
    list_excesses[[i]] <- simulate_excesses_from_params(params6, list_lags[[i]], Tobs=250L)
  }

  # Pas utilisé quand rpar=TRUE, mais l’API demande un list_episodes
  list_episodes <- vector("list", length(list_lags))
  list_episodes[] <- list(NA)

  # Point de départ (bruité mais raisonnable)
  start <- true
  start[1:4] <- pmax(1e-3, start[1:4] * runif(4, 0.6, 1.4)) # beta/alpha > 0
  start["eta1"] <- start["eta1"] * runif(1, 0.5, 1.5)
  start["eta2"] <- pmax(0.5, start["eta2"] * runif(1, 0.7, 1.3))

  # NLL au départ
  nll0 <- neg_ll_composite(start, list_episodes, list_excesses, list_lags,
                           wind_df=wind_df, rpar=TRUE, distance="euclidean")

  # Optimisation L-BFGS-B avec bornes logiques (positivité pour {beta, alpha, eta2})
  lower <- c(1e-6, 1e-6, 0.5, 0.5, -Inf, 0.5)
  upper <- c(  5.0,   5.0, 2.0, 2.0,  Inf, 2.5)

  opt <- optim(
    par     = start,
    fn      = function(par) neg_ll_composite(par, list_episodes, list_excesses, list_lags,
                                             wind_df=wind_df, rpar=TRUE, distance="euclidean"),
    method  = "L-BFGS-B",
    lower   = lower,
    upper   = upper,
    control = list(maxit = 300, factr = 1e7) # convergence tolérante et rapide
  )

  # Vérifs: l’optim baisse la NLL
  expect_true(is.finite(opt$value))
  expect_lt(opt$value, nll0)

  # Paramètres estimés ~ vrais (tolérances réalistes, les eta peuvent être plus durs)
  est <- opt$par
  # Tolerances: betas et alphas ~10-20%, etas un peu plus large
  rel_err <- function(est, tru) abs(est - tru) / pmax(1e-8, abs(tru))

  re_ba <- rel_err(est[1:4], true[1:4])
  expect_true(all(re_ba < 0.25))  # 25% de relatif : robuste sur petits jeux

  re_eta <- rel_err(est[5:6], true[5:6])
  expect_true(all(re_eta < 0.40)) # etas un peu plus souples

  # Sanity: remettre les estims dans la NLL doit être <= au point de départ
  nll_hat <- neg_ll_composite(est, list_episodes, list_excesses, list_lags,
                              wind_df=wind_df, rpar=TRUE)
  expect_lte(nll_hat, nll0)
})

# --- Test 2: cas sans vent (params à 4) --------------------------------------

test_that("optim works without wind (params length 4) and reduces NLL", {
  set.seed(777)

  true4 <- c(beta1=0.9, beta2=0.7, alpha1=1.1, alpha2=0.8) # pas d’advection
  # 2 épisodes, pas de wind_df
  list_lags <- list(
    make_random_df_lags(25, seed=11),
    make_random_df_lags(28, seed=22)
  )

  # Simule excès directement avec params_adv = (true4, adv=0)
  list_excesses <- lapply(list_lags, function(L) {
    simulate_excesses_from_params(c(true4, 0, 0), L, Tobs=200L)
  })
  list_episodes <- list(NA, NA)

  start <- c(beta1=0.6, beta2=1.0, alpha1=0.9, alpha2=0.9) # point de départ

  nll0 <- neg_ll_composite(start, list_episodes, list_excesses, list_lags,
                           wind_df=NA, rpar=TRUE)

  opt <- optim(
    par     = start,
    fn      = function(par) neg_ll_composite(par, list_episodes, list_excesses, list_lags,
                                             wind_df=NA, rpar=TRUE),
    method  = "L-BFGS-B",
    lower   = c(1e-6, 1e-6, 0.5, 0.5),
    upper   = c(  5.0,   5.0, 2.0, 2.0),
    control = list(maxit=200)
  )

  expect_true(is.finite(opt$value))
  expect_lt(opt$value, nll0)

  # paramètres proches (tolérance 25% réaliste)
  rel_err <- function(est, tru) abs(est - tru) / pmax(1e-8, abs(tru))
  expect_true(all(rel_err(opt$par, true4) < 0.25))
})

# --- Test 3: normalisation TRUE vs FALSE n’empêche pas l’optim de converger ---

test_that("optimization also works with normalize=TRUE", {
  set.seed(999)
  true <- c(beta1=0.5, beta2=0.6, alpha1=1.3, alpha2=0.7, eta1=0.15, eta2=1.2)
  wind_df <- data.frame(vx=c(2, -3), vy=c(0, 1))

  list_lags <- list(
    make_random_df_lags(40, seed=5),
    make_random_df_lags(32, seed=6)
  )

  eta1 <- true["eta1"]; eta2 <- true["eta2"]
  adv_x <- eta1 * sign(wind_df$vx) * abs(wind_df$vx)^eta2
  adv_y <- eta1 * sign(wind_df$vy) * abs(wind_df$vy)^eta2

  list_excesses <- list(
    simulate_excesses_from_params(c(true[1:4], adv_x[1], adv_y[1]), list_lags[[1]], Tobs=220L),
    simulate_excesses_from_params(c(true[1:4], adv_x[2], adv_y[2]), list_lags[[2]], Tobs=220L)
  )
  list_episodes <- list(NA, NA)

  start <- c(beta1=0.4, beta2=0.8, alpha1=1.1, alpha2=0.8, eta1=0.1, eta2=1.0)

  nll0 <- neg_ll_composite(start, list_episodes, list_excesses, list_lags,
                           wind_df=wind_df, rpar=TRUE, normalize=TRUE)

  opt <- optim(
    par     = start,
    fn      = function(par) neg_ll_composite(par, list_episodes, list_excesses, list_lags,
                                             wind_df=wind_df, rpar=TRUE, normalize=TRUE),
    method  = "L-BFGS-B",
    lower   = c(1e-6,1e-6,0.5,0.5, -Inf, 0.5),
    upper   = c(  5.0,  5.0,2.0,2.0,  Inf, 2.5),
    control = list(maxit=250)
  )

  expect_lt(opt$value, nll0)
  # tolérances plus larges car normalize change l’échelle interne
  rel_err <- function(est, tru) abs(est - tru) / pmax(1e-8, abs(tru))
  expect_true(all(rel_err(opt$par[1:4], true[1:4]) < 0.35))
})

