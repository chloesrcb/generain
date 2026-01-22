#' pEGPD_full
#' 
#' CDF of the entire distribution with a mass at zero and an EGPD for
#' positive values
#' @param x numeric vector of quantiles
#' @param p0 probability of zero
#' @param xi shape parameter of the EGPD
#' @param sigma scale parameter of the EGPD
#' @param kappa second shape parameter of the EGPD
#' @return numeric vector of CDF values
#' @export
pEGPD_full <- function(x, p0, xi, sigma, kappa) {
    v <- numeric(length(x))
    v[x <= 0] <- 0
    v[x == 0] <- p0
    idx <- which(x > 0)
    if(length(idx) > 0) {
      v[idx] <- p0 + (1 - p0) * pextgp(x[idx], xi = xi, sigma = sigma, kappa = kappa)
    }
    v <- pmin(pmax(v, 1e-12), 1 - 1e-12)
    return(v)
}

#' qEGPD_full
#' 
#' Quantile function of the entire distribution with a mass at zero and an EGPD for
#' positive values
#' @param v numeric vector of probabilities
#' @param p0 probability of zero
#' @param xi shape parameter of the EGPD
#' @param sigma scale parameter of the EGPD
#' @param kappa second shape parameter of the EGPD
#' @return numeric vector of quantiles
#' @export
qEGPD_full <- function(v, p0, xi, sigma, kappa) {
  v <- pmin(pmax(v, 1e-12), 1 - 1e-12)

  out <- numeric(length(v))
  
  na_idx <- is.na(v)
  out[na_idx] <- NA
  
  dry  <- (v <= p0) & !na_idx
  wet  <- (v > p0) & !na_idx

  # rain = 0 if U ≤ p0
  out[dry] <- 0

  # if > p0
  if (any(wet)) {
    x <- (v[wet] - p0) / (1 - p0)  # tail prob
    out[wet] <- qextgp(x, type = 1, xi = xi, sigma = sigma, kappa = kappa)
  }

  return(out)
}

#' G_std
#' 
#' Standardized CDF transformation function to get rainfall scale from latent 
#' pareto scale
#' @param z numeric vector of latent variables
#' @param p0 probability of zero
#' @param u threshold in the latent space
#' @return numeric vector of transformed variables
#' @export
G_std <- function(z, p0, u) {
  v <- numeric(length(z))
  v[] <- NA_real_
  z <- as.numeric(z)
  x0 <- 2 / (1 - p0)
  
  # z < 0
  idx_neg <- which(z < 0)
  if (length(idx_neg)) v[idx_neg] <- 0
  
  # 0 <= z <= x0 : G(z) = p0 + a*z
  a <- (1 - p0)^2 / 4
  idx1 <- which(z >= 0 & z <= x0)
  if (length(idx1)) v[idx1] <- p0 + a * z[idx1]
  
  # x0 < z <= u : linear interpolation between v_x0 and 1 - 1/u
  idx2 <- which(z > x0 & z <= u)
  if (length(idx2)) {
    v_x0 <- p0 + a * x0   # = (1 + p0)/2
    v_u  <- 1 - 1 / u
    # if u == x0 then no interval; treat numerically
    if (abs(u - x0) < .Machine$double.eps) {
      v[idx2] <- v_u
    } else {
      slope2 <- (v_u - v_x0) / (u - x0)
      v[idx2] <- v_x0 + slope2 * (z[idx2] - x0)
    }
  }
  
  # region z > u : tail
  idx3 <- which(z > u)
  if (length(idx3)) v[idx3] <- 1 - 1 / z[idx3]
  
  # clamps numeriques
  v <- pmin(pmax(v, 1e-12), 1 - 1e-12)
  return(v)
}


#' G_std_inv
#' 
#' Inverse of the standardized CDF transformation function to get latent 
#' pareto scale from rainfall scale
#' @param v numeric vector of transformed variables
#' @param p0 probability of zero
#' @param u threshold in the latent space
#' @return numeric vector of latent variables
#' @export
G_std_inv <- function(v, p0, u) {
  v <- pmin(pmax(v, 1e-12), 1 - 1e-12)
  z <- numeric(length(v))
  x0 <- 2 / (1 - p0)
  a <- (1 - p0)^2 / 4
  v_x0 <- p0 + a * x0  # = (1 + p0)/2
  v_switch <- 1 - 1 / u
  
  z[v <= 0] <- -Inf
  
  idx_mass0 <- which(v > 0 & v <= p0)
  if (length(idx_mass0)) z[idx_mass0] <- 0
  
  idx_low <- which(v > p0 & v <= v_x0)
  if (length(idx_low)) z[idx_low] <- (v[idx_low] - p0) / a
  
  idx_mid <- which(v > v_x0 & v <= v_switch)
  if (length(idx_mid)) {
    if (abs(u - x0) < .Machine$double.eps) {
      z[idx_mid] <- u
    } else {
      slope2 <- (v_switch - v_x0) / (u - x0)
      z[idx_mid] <- x0 + (v[idx_mid] - v_x0) / slope2
    }
  }
  
  idx_high <- which(v > v_switch)
  if (length(idx_high)) z[idx_high] <- 1 / (1 - v[idx_high])
  
  return(z)
}


#' sim_episode_coords
#' 
#' Simulate a single episode over given coordinates and times
#' with r-Pareto process and transform to observed scale
#' @param params_vario list of variogram parameters
#' @param params_margins list of marginal parameters
#' @param coords matrix of coordinates (sites x 2)
#' @param times vector of times
#' @param adv vector of advection (vx, vy)
#' @param t0 time index of the event origin
#' @param s0 site name of the event origin
#' @param u threshold in the latent space
#' @param u_emp empirical threshold in the observed space
#' @param plot_debug logical, whether to plot debug plots
#' @param filename filename prefix to save debug plots if plot_debug is TRUE
#' @return list with latent Z, observed X, and latent threshold u_latent
#' @export
sim_episode_coords <- function(params_vario, params_margins,
                               coords, times, adv, t0, s0,
                               u, u_emp,
                               plot_debug = FALSE, filename = NULL) {
  idx_s0 <- which(names(params_margins$p0) == s0)
  p0_s0 <- params_margins$p0[idx_s0]
  xi_s0 <- params_margins$xi[idx_s0]
  sigma_s0 <- params_margins$sigma[idx_s0]
  kappa_s0 <- params_margins$kappa[idx_s0]
  # value at site s0 corresponding to u_emp
  x_s0 <- pEGPD_full(u_emp,
                      p0 = p0_s0,
                      xi = xi_s0,
                      sigma = sigma_s0,
                      kappa = kappa_s0)
  # latent threshold
  u <- G_std_inv(x_s0, p0 = p0_s0, u = u)
  s0_index <- which(rownames(coords) == s0)
  # simulate r-Pareto process
  sim <- sim_rpareto_coords(
    beta1 = params_vario$beta1,
    beta2 = params_vario$beta2,
    alpha1 = params_vario$alpha1,
    alpha2 = params_vario$alpha2,
    adv    = adv,
    coords = coords,
    times      = times,
    t0_index     = t0 + 1,
    s0_index     = s0_index,
    threshold = u
  )
  
  Z <- sim$Z # latent r-Pareto processs
  nS <- nrow(coords)
  nT <- length(times)
  
  X <- matrix(NA_real_, nS, nT, dimnames = list(rownames(coords), NULL))
  V <- matrix(NA_real_, nS, nT, dimnames = list(rownames(coords), NULL))
  for (k in seq_len(nS)) {
    # marginal parameters for site k
    p0    <- params_margins$p0[k]
    xi    <- params_margins$xi[k]
    sigma <- params_margins$sigma[k]
    kappa <- params_margins$kappa[k]

    Zk <- Z[k, ] # site k time series
    
    # transformed to uniform scale
    V[k, ] <- G_std(Zk, p0 = p0, u = u)

    # Final rainfall values
    X[k, ] <- qEGPD_full(V[k, ], p0, xi, sigma, kappa)

  }

  if (plot_debug) {
    for (s in rownames(Z)) {
      if(is.null(filename)) {
        filename <- "plot_transformation_"
      } 
      plot_transformation_gg(Z, X, u, site_name = s, save_plot = TRUE,
              filename = paste0(im_folder, "swg/omsev/", filename, s, ".png"))
    }
  }
  
  return(list(Z = Z, X = X, u_latent = u))
}



#' compute_p0_episode
#' 
#' Compute the probability of zero for each site in a single episode
#' @param episode matrix of observed values (sites x time)
#' @return named numeric vector of p0 values for each site
#' @export
compute_p0_episode <- function(episode) {
  out <- numeric(ncol(episode))
  names(out) <- colnames(episode)

  for(s in colnames(episode)) {

    x <- episode[, s]

    if (all(is.na(x))) {
      out[s] <- NA
    } else {
      out[s] <- mean(x == 0, na.rm = TRUE)
    }
  }

  return(out)
}

#' compute_p0_all_episodes
#' 
#' Compute the probability of zero for each site across multiple episodes
#' @param list_episodes list of matrices of observed values (sites x time)
#' @return list with matrix of p0 values by episode and mean p0 values
#' @export
compute_p0_all_episodes <- function(list_episodes) {

  list_p0 <- lapply(list_episodes, compute_p0_episode)
  mat_p0  <- do.call(rbind, list_p0)

  p0_mean <- colMeans(mat_p0, na.rm = TRUE)

  return(list(
    p0_by_episode = mat_p0,
    p0_mean = p0_mean
  ))
}

#' simulate_many_episodes
#' 
#' Simulate multiple episodes of rainfall data
#' @param N number of episodes to simulate
#' @param u latent threshold
#' @param u_emp empirical threshold
#' @param params_vario variogram parameters
#' @param params_margins marginal parameters
#' @param coords coordinates of sites
#' @param times time points
#' @param adv advection parameters
#' @param t0 initial time index
#' @param s0 initial site
#' @param plot_debug whether to plot debug information
#' @return list of simulated episodes
#' @export
simulate_many_episodes <- function(N, u, u_emp, params_vario, params_margins,
                                   coords, times, adv, t0, s0, plot_debug = FALSE) {
  sims <- vector("list", N)
  for (i in seq_len(N)) {
    # print progress every 50
    if (i %% 50 == 0) cat("Sim", i, " / ", N, "\n")
    Xsim <- sim_episode_coords(
      params_vario   = params_vario,
      params_margins = params_margins,
      coords         = coords,
      times          = times,
      adv            = adv,
      t0             = t0,
      s0             = s0,
      u              = u,
      u_emp           = u_emp,
      plot_debug     = FALSE
    )
    sims[[i]] <- Xsim
  }
  return(sims)
}


#' sim_episode_grid
#' 
#' Simulate a single episode over given regular grid coordinates and times
#' with r-Pareto process and transform to observed scale
#' @param params_vario list of variogram parameters
#' @param params_margins_common list of common marginal parameters
#' @param coords matrix of coordinates (sites x 2)
#' @param times vector of times
#' @param adv vector of advection (vx, vy)
#' @param t0 time index of the event origin
#' @param s0_pixel_id site name of the event origin
#' @param u threshold in the latent space
#' @param u_emp empirical threshold in the observed space
#' @param plot_debug logical, whether to plot debug plots
#' @param filename filename prefix to save debug plots if plot_debug is TRUE
#' @return matrix of observed X values (sites x time)
#' @export
sim_episode_grid <- function(params_vario, params_margins_common,
                        coords, times, adv, t0, s0_pixel_id,
                        u, u_emp,
                        plot_debug = FALSE, filename = NULL) {

  s0_coords <- as.numeric(coords[rownames(coords) == s0_pixel_id, ])

  x_s0 <- pEGPD_full(u_emp,
                     p0    = params_margins_common$p0,
                     xi    = params_margins_common$xi,
                     sigma = params_margins_common$sigma,
                     kappa = params_margins_common$kappa)
  u <- G_std_inv(x_s0, p0 = params_margins_common$p0, u = u)
  s0_index <- which(rownames(coords) == s0_pixel_id)

  sim <- sim_rpareto_coords(
    beta1 = params_vario$beta1,
    beta2 = params_vario$beta2,
    alpha1 = params_vario$alpha1,
    alpha2 = params_vario$alpha2,
    adv    = adv,
    coords = coords,
    t      = times,
    t0_index     = t0 + 1,
    s0_index     = s0_index,
    threshold = u
  )

  Z <- sim$Z
  nS <- nrow(coords)
  nT <- length(times)

  X <- matrix(NA_real_, nS, nT, dimnames = list(rownames(coords), NULL))
  V <- matrix(NA_real_, nS, nT, dimnames = list(rownames(coords), NULL))

  for (k in seq_len(nS)) {
    Zk <- Z[k, ]
    V[k, ] <- G_std(Zk, p0 = params_margins_common$p0, u = u)
    X[k, ] <- qEGPD_full(V[k, ],
                         p0    = params_margins_common$p0,
                         xi    = params_margins_common$xi,
                         sigma = params_margins_common$sigma,
                         kappa = params_margins_common$kappa)
  }

  return(X)
}


#' plot_th_emp_chi
#' 
#' Plot empirical vs theoretical chi values for multiple episodes
#' @param list_lags list of data frames with lag information for each episode
#' @param list_excesses list of data frames with excess information for each episode
#' @param list_adv data frame with advection parameters for each episode
#' @param params_estimates numeric vector of estimated model parameters
#' @param tau_min minimum tau value to consider
#' @param tau_fixed fixed tau value for specific plot
#' @param h_breaks breaks for h binning
#' @param latlon logical, whether coordinates are in lat/lon
#' @param adv_transform logical, whether to apply advection transformation
#' @return invisible list with results and plots
#' @export
plot_th_emp_chi <- function(list_lags,
                            list_excesses,
                            list_adv,
                            params_estimates,
                            tau_min = 0,
                            tau_fixed = 0,
                            h_breaks = seq(0, 10, by = 1),
                            latlon = FALSE,
                            adv_transform = TRUE) {
  stopifnot(
    length(list_lags) == length(list_excesses),
    length(list_lags) <= nrow(list_adv)
  )

  if (!requireNamespace("data.table", quietly = TRUE)) {
    stop("Package 'data.table' is required.")
  }
  if (!requireNamespace("dplyr", quietly = TRUE)) {
    stop("Package 'dplyr' is required.")
  }
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required.")
  }

  params <- params_estimates
  eta1 <- params[5]
  eta2 <- params[6]

  res_list <- vector("list", length(list_lags))

  for (i in seq_along(list_lags)) {
    lags_i    <- list_lags[[i]]
    excess_i  <- list_excesses[[i]]
    adv_x     <- list_adv$vx[i]
    adv_y     <- list_adv$vy[i]

    if (adv_transform) {
      adv_norm <- sqrt(adv_x^2 + adv_y^2)
      adv_norm_transformed <- eta1 * adv_norm^eta2
      if (adv_norm > 0) {
        adv_x <- adv_x / adv_norm * adv_norm_transformed
        adv_y <- adv_y / adv_norm * adv_norm_transformed
      } else {
        adv_x <- 0
        adv_y <- 0
      }
    }
    
    params_adv <- c(params[1:4], adv_x, adv_y)
    chi_th_i <- theoretical_chi(
      params   = params_adv,
      df_lags  = lags_i,
      latlon   = latlon
    )

    res_list[[i]] <- data.table::data.table(
      episode  = i,
      s2       = lags_i$s2,
      tau      = lags_i$tau,
      h        = lags_i$hnorm,
      hx      = lags_i$hx,
      hy      = lags_i$hy,
      chi_emp  = excess_i$kij,
      chi_theo = chi_th_i$chi,
      adv_x    = adv_x,
      adv_y    = adv_y
    )
  }

  res <- data.table::rbindlist(res_list, fill = TRUE)
  res$hbin <- cut(res$h, breaks = h_breaks, include.lowest = TRUE)

  grouped <- dplyr::group_by(
    dplyr::filter(res, .data$tau >= tau_min),
    .data$tau,
    .data$hbin
  )
  res_cmp <- dplyr::summarise(
    grouped,
    chi_emp_bar  = mean(.data$chi_emp, na.rm = TRUE),
    chi_theo_bar = mean(.data$chi_theo, na.rm = TRUE),
    n_pairs = dplyr::n(),
    .groups = "drop"
  )

  res_tau <- dplyr::filter(res_cmp, .data$tau == tau_fixed)

  plot_all <- ggplot2::ggplot(res_cmp, ggplot2::aes(
    x = .data$chi_theo_bar,
    y = .data$chi_emp_bar
  )) +
    ggplot2::geom_point(alpha = 0.6, color = btfgreen) +
    ggplot2::geom_abline(
      slope = 1,
      intercept = 0,
      color = "red",
      linetype = "dashed"
    ) +
    ggplot2::labs(
      title = "Empirical vs Theoretical Chi",
      x = "Theoretical Chi",
      y = "Empirical Chi"
    ) +
    ggplot2::theme_minimal() +
    btf_theme

  plot_tau <- ggplot2::ggplot(res_tau, ggplot2::aes(
    x = .data$chi_theo_bar,
    y = .data$chi_emp_bar
  )) +
    ggplot2::geom_point(alpha = 0.6, color = btfgreen) +
    ggplot2::geom_abline(
      slope = 1,
      intercept = 0,
      color = "red",
      linetype = "dashed"
    ) +
    ggplot2::labs(
      title = paste0("Empirical vs Theoretical Chi (tau = ", tau_fixed, ")"),
      x = "Theoretical Chi",
      y = "Empirical Chi"
    ) +
    ggplot2::theme_minimal() +
    btf_theme

  invisible(list(
    res = res,
    res_cmp = res_cmp,
    res_tau = res_tau,
    plots = list(all = plot_all, tau = plot_tau)
  ))
}


#' plot_transformation_gg
#' 
#' Plot the transformation from latent Z to observed X for a given site
#' @param Z matrix of latent variables (sites x time)
#' @param X matrix of observed variables (sites x time)
#' @param u threshold in the latent space
#' @param site_name name of the site to plot
#' @param save_plot logical, whether to save the plot
#' @param filename filename to save the plot if save_plot is TRUE
#' @return ggplot object
#' @export
plot_transformation_gg <- function(Z, X, u, site_name,
                              save_plot = FALSE, filename = NULL) {
  
  s <- site_name
  Zs <- Z[s, ]
  Xs <- X[s, ]
  
  ext     <- Zs > u
  interm <- (Zs > 0) & (Zs <= u)
  low  <- !ext & !interm
  
  df <- data.frame(
    time = seq_along(Zs),
    Z    = Zs,
    X    = Xs,
    type = case_when(
      ext      ~ "extreme",
      interm   ~ "intermediate",
      low      ~ "low"
    )
  )
  
  df_long <- df %>%
    pivot_longer(cols = c(X), names_to = "variable", values_to = "value")
  
  df_long <- df_long %>%
    mutate(plot_type = case_when(
      variable == "X" & type == "extreme" ~ "Z extreme",
      variable == "X" & type == "intermediate" ~ "Z intermediate",
      variable == "X" & type == "low" ~ "Z low"
    ))
  

  gg <- ggplot(df_long, aes(x = time, y = value, color = plot_type, shape = plot_type)) +
    geom_point(data = df_long %>% filter(variable == "X"), size = 2) +
    scale_color_manual(values = c("Z extreme" = "#a72909", "Z intermediate" = "#f4a261", "Z low" = btfgreen)) +
    scale_shape_manual(values = c("Z extreme" = 19, "Z intermediate" = 19, "Z low" = 19)) +
    labs(
        x = "Time", y = "X value") +
    theme_minimal() +
    theme(legend.title = element_blank()) +
    btf_theme
  
  print(gg)
  if (save_plot & !is.null(filename)) {
    ggsave(filename, plot = gg, width = 12, height = 6, dpi = 300)
  }

}



#' sum_over_time_by_site
#' 
#' Sum values over time for each site in an episode
#' @param ep_mat matrix of episode values (time x site) or (site x time)
#' @param sites_ref vector of site names
#' @return named numeric vector of summed values by site
#' @export
sum_over_time_by_site <- function(ep_mat, sites_ref) { 
  # ep_mat can be [time x site] or [site x time] 
  x <- as.matrix(ep_mat) 
  if (!is.null(colnames(x)) && all(colnames(x) %in% sites_ref)) { 
    # time x site 
    s <- colSums(x, na.rm = TRUE) 
    return(s) 
  } 
  if (!is.null(rownames(x)) && all(rownames(x) %in% sites_ref)) { 
    # site x time 
    s <- rowSums(x, na.rm = TRUE) 
    return(s) 
  } 
  if (ncol(x) == length(sites_ref)) { 
      s <- colSums(x, na.rm = TRUE) 
      names(s) <- sites_ref
   } else if (nrow(x) == length(sites_ref)) { 
    s <- rowSums(x, na.rm = TRUE) 
    names(s) <- sites_ref 
    } 
  return(s) 
} 
  
# sum_over_time_by_site <- function(ep_mat, sites_ref) {

#   x <- as.matrix(ep_mat)

#   # ensure orientation: time x site
#   if (ncol(x) != length(sites_ref) && nrow(x) == length(sites_ref)) {
#     x <- t(x)
#   }

#   stopifnot(ncol(x) == length(sites_ref))

#   apply(x, 2, function(v) {
#     if (all(is.na(v))) {
#       NA_real_
#     } else {
#       sum(v, na.rm = TRUE)
#     }
#   }) |> setNames(sites_ref)
# }


#' make_qq_df
#' 
#' Create a data frame for QQ plot from observed and simulated values
#' @param x_obs numeric vector of observed values
#' @param x_sim numeric vector of simulated values
#' @param min_n minimum number of values required to create the QQ data frame
#' @return data frame with quantiles of observed and simulated values
#' @export
make_qq_df <- function(x_obs, x_sim, min_n = 50) { 
  x_obs <- x_obs[is.finite(x_obs)] 
  x_sim <- x_sim[is.finite(x_sim)] 
  n <- min(length(x_obs), length(x_sim)) 
  if (n < min_n) return(NULL) 
  p <- (seq_len(n) - 0.5) / n 
  data.frame( 
    q_obs = as.numeric(quantile(x_obs, probs = p, names = FALSE)), 
    q_sim = as.numeric(quantile(x_sim, probs = p, names = FALSE)) 
  ) 
}
