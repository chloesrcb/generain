
dt_hours <- 5/60
tau_lags <- 0:10
h_breaks_m <- seq(0, 1500, by = 150)

dist_mat <- get_dist_mat(location_gauges) / 1000
df_dist <- reshape_distances(dist_mat)
n_hbins <- 10
h_all <- df_dist$value

h_breaks_omsev <- quantile(
  h_all,
  probs = seq(0, 1, length.out = n_hbins + 1),
  na.rm = TRUE
)

h_breaks_omsev <- unique(as.numeric(h_breaks_omsev))
h_breaks_omsev[length(h_breaks_omsev)] <- 1.6
h_breaks_m <- h_breaks_omsev * 1000  # convert to meters

taus_minutes_to_plot <- seq(0, 50, by = 5)
taus_hours_to_plot   <- taus_minutes_to_plot / 60


make_dist_matrix <- function(coords_m) {
  coords_m <- as.data.frame(coords_m)
  xy <- as.matrix(coords_m[, c("Longitude", "Latitude")])
  D <- as.matrix(dist(xy))
  rownames(D) <- rownames(coords_m)
  colnames(D) <- rownames(coords_m)
  D
}


chi_one_episode_idx <- function(exceed_ep, D, s0, t0 = 1,
                                tau_lags, h_breaks) {
  sites <- colnames(exceed_ep)
  if (is.null(sites)) stop("exceed_ep must have colnames = site names")
  if (!(s0 %in% sites)) return(NULL)
  if (!(s0 %in% rownames(D))) return(NULL)

  out <- expand.grid(
    tau_lag = tau_lags,
    hbin    = seq_len(length(h_breaks) - 1)
  )
  out$K <- 0L
  out$N <- 0L

  nT <- nrow(exceed_ep)

  # Distances from conditioning site to all sites present in episode
  d0 <- D[s0, sites]
  hb_all <- findInterval(d0, h_breaks, rightmost.closed = TRUE)

  for (lag in tau_lags) {
    tt <- t0 + lag
    if (tt < 1 || tt > nT) next

    y <- exceed_ep[tt, ]  # 0/1 exceedances at time t0+lag

    for (b in seq_len(length(h_breaks) - 1)) {
      idx_sites <- which(hb_all == b)
      if (!length(idx_sites)) next

      rowi <- which(out$tau_lag == lag & out$hbin == b)
      out$N[rowi] <- out$N[rowi] + length(idx_sites)
      out$K[rowi] <- out$K[rowi] + sum(y[idx_sites], na.rm = TRUE)
    }
  }

  out
}


compute_chi_obs <- function(list_episodes, s0_list, u_list, coords_m,
                            tau_lags, h_breaks, t0 = 1, dt_hours = 5/60) {

  D <- make_dist_matrix(coords_m)

  agg <- expand.grid(tau_lag = tau_lags, hbin = seq_len(length(h_breaks)-1))
  agg$K <- 0L
  agg$N <- 0L

  for (e in seq_along(list_episodes)) {
    ep <- as.matrix(list_episodes[[e]])
    sites <- colnames(ep)
    s0 <- s0_list[e]
    u0 <- u_list[e]

    if (is.null(sites) || !(s0 %in% sites)) next
    if (!is.finite(u0)) next

    # Episode exceedances X > u0
    exceed_ep <- (ep > u0) * 1L

    tmp <- chi_one_episode_idx(exceed_ep, D, s0 = s0, t0 = t0,
                               tau_lags = tau_lags, h_breaks = h_breaks)
    if (is.null(tmp)) next

    keyA <- paste(agg$tau_lag, agg$hbin)
    keyT <- paste(tmp$tau_lag, tmp$hbin)
    m <- match(keyT, keyA)

    agg$K[m] <- agg$K[m] + tmp$K
    agg$N[m] <- agg$N[m] + tmp$N
  }

  agg$chi_hat   <- with(agg, ifelse(N > 0, K / N, NA_real_))
  agg$h_mid_km   <- (h_breaks[agg$hbin] + h_breaks[agg$hbin + 1]) / 2
#   agg$h_mid_km  <- agg$h_mid_m / 1000
  agg$tau_hours <- agg$tau_lag * dt_hours

  agg
}


compute_chi_sim <- function(sims_by_ep, s0_list, u_list, coords_m,
                            tau_lags, h_breaks, t0 = 1, dt_hours = 5/60) {

  D <- make_dist_matrix(coords_m)

  agg <- expand.grid(tau_lag = tau_lags, hbin = seq_len(length(h_breaks)-1))
  agg$K <- 0L
  agg$N <- 0L

  for (j in seq_along(sims_by_ep)) {
    s0 <- s0_list[j]
    u0 <- u_list[j]
    if (!is.finite(u0)) next

    reps <- sims_by_ep[[j]]
    if (!length(reps)) next

    for (m in seq_along(reps)) {
      X <- reps[[m]]$X
      if (is.null(X)) next

      ep <- as.matrix(X) # time x site
      sites <- colnames(ep)
      if (is.null(sites) || !(s0 %in% sites)) next

      exceed_ep <- (ep > u0) * 1L

      tmp <- chi_one_episode_idx(exceed_ep, D, s0 = s0, t0 = t0,
                                 tau_lags = tau_lags, h_breaks = h_breaks)
      if (is.null(tmp)) next

      keyA <- paste(agg$tau_lag, agg$hbin)
      keyT <- paste(tmp$tau_lag, tmp$hbin)
      mm <- match(keyT, keyA)

      agg$K[mm] <- agg$K[mm] + tmp$K
      agg$N[mm] <- agg$N[mm] + tmp$N
    }
  }

  agg$chi_hat   <- with(agg, ifelse(N > 0, K / N, NA_real_))
  agg$h_mid_km   <- (h_breaks[agg$hbin] + h_breaks[agg$hbin + 1]) / 2
#   agg$h_mid_km  <- agg$h_mid_m / 1000
  agg$tau_hours <- agg$tau_lag * dt_hours

  agg
}


plot_chi_curves_obs_vs_sim <- function(chi_obs_df, chi_sim_df,
                                       taus_hours_to_plot,
                                       min_N_obs = 10, min_K_obs = 1,
                                       min_N_sim = 10) {

  obs <- subset(chi_obs_df, N >= min_N_obs & K >= min_K_obs & is.finite(chi_hat))
  sim <- subset(chi_sim_df, N >= min_N_sim & is.finite(chi_hat))

  taus_hours_to_plot <- sort(unique(obs$tau_hours))

  op <- par(no.readonly = TRUE)
  nT <- length(taus_hours_to_plot)
  par(mfrow = c(ceiling(nT/3), min(3, nT)), mar = c(4,4,2,1))

  for (th in taus_hours_to_plot) {
    o <- obs[obs$tau_hours == th, ]
    s <- sim[sim$tau_hours == th, ]

    hs <- sort(intersect(o$h_mid_km, s$h_mid_km))
    o <- o[match(hs, o$h_mid_km), ]
    s <- s[match(hs, s$h_mid_km), ]
    ss <- smooth_sim_curve(s)
    so <- smooth_sim_curve(o)
    plot(o$h_mid_km, o$chi_hat, pch = 16,
         xlab = "Distance h (km)",
         ylab = expression(hat(chi)[r](h, tau)),
         main = paste0("tau = ", round(th*60), " min"),
         ylim = range(c(o$chi_hat, s$chi_hat), na.rm = TRUE))
        # ylim = c(0, 0.8))
        

    points(s$h_mid_km, s$chi_hat, lwd = 2)
    lines(ss$h_mid_km, ss$chi_smooth, lwd = 2)
    lines(so$h_mid_km, so$chi_smooth, lwd = 2, lty = 2)
    legend("topright",
           legend = c("Observation", "Simulation"),
           pch = c(16, NA), lty = c(NA, 1), lwd = c(NA, 2),
           bty = "n")
  }
}



plot_emp_vs_emp_scatter <- function(chi_obs_df, chi_sim_df,
                                    min_N_obs = 10, min_K_obs = 1,
                                    min_N_sim = 10) {

  obs <- subset(chi_obs_df, N >= min_N_obs & K >= min_K_obs & is.finite(chi_hat))
  sim <- subset(chi_sim_df, N >= min_N_sim & is.finite(chi_hat))

  m <- merge(obs, sim, by = c("tau_hours", "hbin"), suffixes = c("_obs", "_sim"))
  if (!nrow(m)) stop("No common bins after filtering. Lower thresholds or adjust bins.")

  plot(m$chi_hat_sim, m$chi_hat_obs, pch = 16,
       xlab = expression(hat(chi)[sim]),
       ylab = expression(hat(chi)[obs]))

  abline(0, 1, lty = 2, lwd = 2)
  invisible(m)
}

h_breaks_m <- seq(0, 1600, by = 100)
h_breaks <- h_breaks_m 
chi_obs_df <- compute_chi_obs(
  list_episodes = list_episodes,
  s0_list       = s0_list,
  u_list        = u_list,
  coords_m      = df_coords,
  tau_lags      = tau_lags,
  h_breaks      = h_breaks,
  t0            = 1,
  dt_hours      = dt_hours
)

chi_sim_df <- compute_chi_sim(
  sims_by_ep = sims_by_ep,
  s0_list    = s0_list,
  u_list     = u_list,
  coords_m   = df_coords,
  tau_lags   = tau_lags,
  h_breaks   = h_breaks,
  t0         = 1,
  dt_hours   = dt_hours
)

plot_chi_curves_obs_vs_sim(
  chi_obs_df, chi_sim_df,
  taus_hours_to_plot = taus_hours_to_plot,
  min_N_obs = 20, min_K_obs = min_K_obs,
  min_N_sim = 100
)

m_emp_emp <- plot_emp_vs_emp_scatter(
  chi_obs_df, chi_sim_df,
  min_N_obs = 20, min_K_obs = min_K_obs,
  min_N_sim = 100
)


smooth_sim_curve <- function(df_sim, span = 0.6) {
  ord <- order(df_sim$h_mid_km)
  lo <- loess(chi_hat ~ h_mid_km,
              data = df_sim[ord, ],
              span = span,
              degree = 1)
  data.frame(
    h_mid_km = df_sim$h_mid_km[ord],
    chi_smooth = predict(lo)
  )
}
