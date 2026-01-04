gamma_st <- function(dx, dy, dt, beta1, beta2, alpha1, alpha2, adv) {
  hx <- dx - dt * adv[1]
  hy <- dy - dt * adv[2]
  beta1 * (sqrt(hx^2 + hy^2))^alpha1 + beta2 * abs(dt)^alpha2
}

sim_W_from_variogram <- function(coords, times,
                                 beta1, beta2, alpha1, alpha2,
                                 adv = c(0,0),
                                 ref_site = 1, ref_time = 1,
                                 seed = NULL, jitter = 1e-10) {
  if (!is.null(seed)) set.seed(seed)

  n_sites <- nrow(coords)
  lt <- length(times)
  n <- n_sites * lt

  grid <- expand.grid(site = seq_len(n_sites), it = seq_len(lt))
  x  <- coords$Longitude[grid$site]
  y  <- coords$Latitude [grid$site]
  tt <- times[grid$it]

  ref <- (ref_time - 1) * n_sites + ref_site
  x0 <- x[ref]; y0 <- y[ref]; t0 <- tt[ref]

  g_ref <- gamma_st(x - x0, y - y0, tt - t0, beta1, beta2, alpha1, alpha2, adv)

  C <- matrix(0, n, n)
  for (i in 1:n) {
    for (j in i:n) {
      gij <- gamma_st(x[i]-x[j], y[i]-y[j], tt[i]-tt[j],
                      beta1, beta2, alpha1, alpha2, adv)
      Cij <- g_ref[i] + g_ref[j] - gij
      C[i,j] <- Cij
      C[j,i] <- Cij
    }
  }

  diag(C) <- diag(C) + jitter
  R <- chol(C)              # R'R = C (R upper)
  z <- rnorm(n)
  W_vec <- as.vector(t(R) %*% z)  # cov = C

  W <- matrix(W_vec, nrow = n_sites, ncol = lt, byrow = FALSE)
  W[ref_site, ref_time] <- 0  # numeric safety
  list(W = W, C = C, ref = ref)
}


sim_rpareto_coords <- function(coords, times,
                               beta1, beta2, alpha1, alpha2,
                               adv = c(0,0),
                               threshold = 1,
                               s0_index = 1, t0_index = 1,
                               seed = NULL) {

  simW <- sim_W_from_variogram(coords, times,
                               beta1, beta2, alpha1, alpha2,
                               adv = adv,
                               ref_site = s0_index, ref_time = t0_index,
                               seed = seed)

  W <- simW$W
  n_sites <- nrow(coords)
  lt <- length(times)

  x0 <- coords$Longitude[s0_index]
  y0 <- coords$Latitude [s0_index]
  t0 <- times[t0_index]

  gamma0_vec <- gamma_st(
    dx = rep(coords$Longitude, times = lt) - x0,
    dy = rep(coords$Latitude,  times = lt) - y0,
    dt = rep(times, each = n_sites) - t0,
    beta1, beta2, alpha1, alpha2, adv
  )
  gamma0 <- matrix(gamma0_vec, nrow = n_sites, ncol = lt, byrow = FALSE)

  W0 <- W[s0_index, t0_index]
  Y  <- exp(W - W0 - gamma0)

  R <- 1 / runif(1)      # Pareto(1)
  Z <- threshold * R * Y
  rownames(Z) <- rownames(coords)

  list(Z = Z, W = W, gamma0 = gamma0)
}


# theoretical_chi <- function(params, df_lags, latlon = FALSE) {

#   beta1  <- params[1]
#   beta2  <- params[2]
#   alpha1 <- params[3]
#   alpha2 <- params[4]
#   adv_x  <- params[5]
#   adv_y  <- params[6]

#   gamma_vals <- gamma_st(
#     dx = df_lags$dx,
#     dy = df_lags$dy,
#     dt = df_lags$tau,
#     beta1 = beta1,
#     beta2 = beta2,
#     alpha1 = alpha1,
#     alpha2 = alpha2,
#     adv = c(adv_x, adv_y)
#   )

#   chi <- 2 * (1 - pnorm(sqrt(1/2 * gamma_vals)))

#   return(data.frame(chi = chi))
# }



# coords <- data.frame(
#   Longitude = c(0, 5),
#   Latitude  = c(0, 0)
# )
# times <- 0:10

# beta1 <- 0.2
# beta2 <- 0.1
# alpha1 <- 1.5
# alpha2 <- 1.0
# adv <- c(2, 10)   # advection forte


# iA <- 1; tA <- 2
# iB <- 2; tB <- 7


# dx <- coords$Longitude[iB] - coords$Longitude[iA]
# dy <- coords$Latitude [iB] - coords$Latitude [iA]
# dt <- times[tB] - times[tA]

# gamma_th <- gamma_st(dx, dy, dt,
#                      beta1, beta2, alpha1, alpha2, adv)

# M <- 2000
# g_emp <- numeric(M)

# for (m in 1:M) {
#   W <- sim_W_from_variogram(
#     coords, times,
#     beta1, beta2, alpha1, alpha2,
#     adv = adv,
#     seed = 1000 + m
#   )$W

#   g_emp[m] <- 0.5 * (W[iA, tA] - W[iB, tB])^2
# }

# c(empirical = mean(g_emp), theoretical = gamma_th)


# taus <- c(1, 3, 5)
# res <- data.frame()

# for (tau in taus) {
#   dx1 <- adv[1] * tau      # aligné
#   dx2 <- adv[1] * tau + 5  # désaligné

#   g1 <- gamma_st(dx1, 0, tau,
#                  beta1, beta2, alpha1, alpha2, adv)
#   g2 <- gamma_st(dx2, 0, tau,
#                  beta1, beta2, alpha1, alpha2, adv)

#   res <- rbind(res, data.frame(
#     tau = tau,
#     gamma_aligned = g1,
#     gamma_shifted = g2
#   ))
# }

# res

# chi_theo <- function(gamma) {
#   2 * (1 - pnorm(sqrt(gamma / 2)))
# }

# chi_th <- chi_theo(gamma_th)


# M <- 3000
# chi_emp <- numeric(M)

# for (m in 1:M) {
#   Z <- sim_rpareto_coords(
#     coords, times,
#     beta1, beta2, alpha1, alpha2,
#     adv = adv,
#     threshold = 1,
#     s0_index = iA,
#     t0_index = tA,
#     seed = 5000 + m
#   )$Z

#   chi_emp[m] <- as.numeric(Z[iB, tB] > 1)
# }

# c(empirical = mean(chi_emp), theoretical = chi_th)


# adv <- c(10, 0)
# beta1 <- 0.05
# beta2 <- 0.02

# coords <- data.frame(
#   Longitude = c(0, 50),
#   Latitude  = c(0, 0)
# )
# times <- 0:20

# iA <- 1; tA <- 2
# iB <- 2; tB <- 12

# dx <- coords$Longitude[iB] - coords$Longitude[iA]
# dt <- times[tB] - times[tA]

# gamma_th <- gamma_st(dx, 0, dt,
#                      beta1, beta2, alpha1, alpha2, adv)

# chi_th <- chi_theo(gamma_th)


# taus <- seq(0, 10, by = 1)
# dxs  <- seq(-20, 20, by = 0.5)

# grid <- expand.grid(dx = dxs, tau = taus)
# grid$gamma <- with(grid,
#   gamma_st(dx, 0, tau,
#            beta1, beta2, alpha1, alpha2, adv)
# )
# grid$chi <- chi_theo(grid$gamma)

# library(ggplot2)
# ggplot(grid, aes(dx, chi, color = factor(tau))) +
#   geom_line() +
#   labs(
#     x = expression(h),
#     y = expression(chi(h,tau)),
#     color = expression(tau),
#     title = "Directional r-extremogram with strong advection"
#   ) +
#   theme_minimal()
