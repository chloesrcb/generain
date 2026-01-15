
qEGPD_marg <- function(u, p0, xi, sigma, kappa) {
  u <- pmin(pmax(u, 1e-12), 1 - 1e-12)
  out <- numeric(length(u))
  idx_pos <- (u > p0)
  out[!idx_pos] <- 0.0
  if (any(idx_pos)) {
    v <- (u[idx_pos] - p0) / (1 - p0)
    out[idx_pos] <- qextgp(
      v, type = 1,
      xi = xi, sigma = sigma, kappa = kappa
    )
  }
  out
}

F_s <- function(x, p0, xi, sigma, kappa) {
  p0 + (1 - p0) * pextgp(x, type = 1, xi = xi, sigma = sigma, kappa = kappa)
}


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

# qEGPD_full <- function(v, p0, xi, sigma, kappa) {
#   # v : matrix n_sites x n_times
#   nS <- nrow(v)
#   nT <- ncol(v)
#   out <- matrix(0, nS, nT)
  
#   for (k in seq_len(nS)) {
#     vk <- v[k, ]
#     pk <- p0[k]
#     xi_k <- xi[k]
#     sigma_k <- sigma[k]
#     kappa_k <- kappa[k]
    
#     dry  <- vk <= pk
#     wet  <- vk > pk
    
#     out[k, dry] <- 0
#     if(any(wet)) {
#       x <- (vk[wet] - pk) / (1 - pk)
#       out[k, wet] <- qextgp(x, type = 1, xi = xi_k, sigma = sigma_k, kappa = kappa_k)
#     }
#   }
  
#   return(out)
# }

# pEGPD_full <- function(x, p0, xi, sigma, kappa) {
#   nS <- nrow(x)
#   nT <- ncol(x)
#   v <- matrix(0, nS, nT)
  
#   for (k in seq_len(nS)) {
#     xk <- x[k, ]
#     pk <- p0[k]
#     xi_k <- xi[k]
#     sigma_k <- sigma[k]
#     kappa_k <- kappa[k]
    
#     idx_pos <- xk > 0
#     v[k, !idx_pos] <- 0
#     if(any(idx_pos)) {
#       v[k, idx_pos] <- pk + (1 - pk) * pextgp(xk[idx_pos], type = 1, xi = xi_k, sigma = sigma_k, kappa = kappa_k)
#     }
#   }
  
#   v <- pmin(pmax(v, 1e-12), 1 - 1e-12)
#   return(v)
# }


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

# G_std <- function(z, p0, u) {
#   v <- numeric(length(z))

#   v[z < 0] <- 0
#   x0 <- 2/(1 - p0)

#   idx1 <- which(z > 0 & z <= x0)
#   v[idx1] <- p0

#   idx2 <- which(z > x0 & z <= u)
#   v[idx2] <- p0 + ( (1 - 1/u) - p0 ) / (u - x0) * (z[idx2] - x0)

#   idx3 <- which(z > u)
#   v[idx3] <- 1 - 1/z[idx3]

#   v <- pmin(pmax(v, 1e-12), 1 - 1e-12)

#   return(v)
# }

# G_std_inv <- function(v, p0, u) {
#   v <- pmin(pmax(v, 1e-12), 1 - 1e-12)
  
#   z <- numeric(length(v))
#   z[v <= 0] <- -Inf
  
#   idx0 <- which(v > 0 & v <= p0)
#   if (length(idx0) > 0) z[idx0] <- 0
  
#   v_switch <- 1 - 1 / u
  
#   x0 <- 2/(1 - p0)
  
#   idx_mid <- which(v > p0 & v <= v_switch)
#   if (length(idx_mid) > 0) {
#     z[idx_mid] <- (4 / (1 - p0)^2) * (v[idx_mid] - p0)
#   }
  
#   idx_high <- which(v > v_switch)
#   if (length(idx_high) > 0) z[idx_high] <- 1 / (1 - v[idx_high])
  
#   return(z)
# }



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


sim_episode_coords <- function(params_vario, params_margins,
                               coords, times, adv, t0, s0,
                               u, u_emp,
                               plot_debug = FALSE, filename = NULL) {
  idx_s0 <- which(names(params_margins$p0) == s0)
  p0_s0 <- params_margins$p0[idx_s0]
  xi_s0 <- params_margins$xi[idx_s0]
  sigma_s0 <- params_margins$sigma[idx_s0]
  kappa_s0 <- params_margins$kappa[idx_s0]
  x_s0 <- pEGPD_full(u_emp,
                      p0 = p0_s0,
                      xi = xi_s0,
                      sigma = sigma_s0,
                      kappa = kappa_s0)
  u <- G_std_inv(x_s0, p0 = p0_s0, u = u)
  print(u)
  s0_coords <- as.numeric(coords[rownames(coords) == s0, ])
  s0_index <- which(rownames(coords) == s0)
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
  
  Z <- sim$Z
  nS <- nrow(coords)
  nT <- length(times)
  
  X <- matrix(NA_real_, nS, nT,
              dimnames = list(rownames(coords), NULL))
  V <- matrix(NA_real_, nS, nT,
              dimnames = list(rownames(coords), NULL))
  extreme_idx <- matrix(FALSE, nS, nT)
  for (k in seq_len(nS)) {

    p0    <- params_margins$p0[k]
    xi    <- params_margins$xi[k]
    sigma <- params_margins$sigma[k]
    kappa <- params_margins$kappa[k]

    Zk <- Z[k, ] # site k time series
    
    V[k, ] <- G_std(Zk, p0 = p0, u = u)

    X[k, ] <- qEGPD_full(V[k, ], p0, xi, sigma, kappa)

  }

  if (plot_debug) {
    for (s in rownames(Z)) {
      if(is.null(filename)) {
        filename <- "plot_transformation_"
      } 
      plot_transformation_gg(Z, X, u, site_name = s,
                             save_plot = TRUE,
                             filename = paste0(im_folder, "swg/omsev/", filename, s, ".png"))
    }
  }
  
  return(list(
    Z = Z,
    X = X,
    u_latent = u
  ))
}




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

compute_p0_all_episodes <- function(list_episodes) {

  list_p0 <- lapply(list_episodes, compute_p0_episode)
  mat_p0  <- do.call(rbind, list_p0)

  p0_mean <- colMeans(mat_p0, na.rm = TRUE)

  return(list(
    p0_by_episode = mat_p0,
    p0_mean = p0_mean
  ))
}


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

  sim <- sim_rpareto_coords(
    beta1 = params_vario$beta1,
    beta2 = params_vario$beta2,
    alpha1 = params_vario$alpha1,
    alpha2 = params_vario$alpha2,
    adv    = adv,
    coords = coords,
    t      = times,
    t0     = t0,
    s0     = s0_coords,
    threshold = u
  )

  Z <- sim$Z[,,1, drop = TRUE]
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



gamma_st <- function(dx, dy, dt, beta1, beta2, alpha1, alpha2, adv) {
  # dx,dy: differences in space (same units as adv*time)
  # dt: time difference
  hx <- dx - dt * adv[1]
  hy <- dy - dt * adv[2]
  2 * beta1 * (sqrt(hx^2 + hy^2))^alpha1 + 2 * beta2 * abs(dt)^alpha2
}


sim_W_from_variogram <- function(coords, times,
                                 beta1, beta2, alpha1, alpha2,
                                 adv = c(0,0), seed = NULL, jitter = 1e-10) {
  if (!is.null(seed)) set.seed(seed)

  n_sites <- nrow(coords)
  lt <- length(times)
  n <- n_sites * lt

  # points p = (x, y, t)
  grid <- expand.grid(site = seq_len(n_sites), it = seq_len(lt))
  x <- coords$Longitude[grid$site]
  y <- coords$Latitude [grid$site]
  tt <- times[grid$it]

  # choose a reference point (conditioning point is convenient)
  ref <- 1L
  x0 <- x[ref]; y0 <- y[ref]; t0 <- tt[ref]

  # gamma(p, ref)
  g_ref <- gamma_st(x - x0, y - y0, tt - t0, beta1, beta2, alpha1, alpha2, adv)

  # build covariance matrix of increments: C_ij = (g_i + g_j - gamma(p_i - p_j))/2
  C <- matrix(0, n, n)
  for (i in 1:n) {
    dx_i <- x[i] - x0; dy_i <- y[i] - y0; dt_i <- tt[i] - t0
    for (j in i:n) {
      dx <- x[i] - x[j]
      dy <- y[i] - y[j]
      dt <- tt[i] - tt[j]
      gij <- gamma_st(dx, dy, dt, beta1, beta2, alpha1, alpha2, adv)
      Cij <- (g_ref[i] + g_ref[j] - gij)
      C[i,j] <- Cij
      C[j,i] <- Cij
    }
  }

  # numerical stabilization
  diag(C) <- diag(C) + jitter

  # sample Gaussian increments: W(ref)=0, others are increments
  L <- chol(C)
  z <- rnorm(n)
  W_inc <- as.vector(t(L) %*% z)  # mean 0, cov C

  # reshape to (sites x times)
  W <- matrix(W_inc, nrow = n_sites, ncol = lt, byrow = FALSE)

  return(list(W = W, C = C, grid = grid, x = x, y = y, tt = tt, ref = ref))
}


sim_rpareto_coords <- function(coords, times,
                                beta1, beta2, alpha1, alpha2,
                                adv = c(0,0),
                                threshold = 1,
                                s0_index = 1, t0_index = 1,
                                seed = NULL) {

  n_sites <- nrow(coords)
  lt <- length(times)

  # simulate W increments with reference = (s0,t0) to simplify
  simW <- sim_W_from_variogram(coords, times, beta1, beta2, alpha1, alpha2, adv, seed = seed)
  W <- simW$W

  x0 <- coords$Longitude[s0_index]
  y0 <- coords$Latitude [s0_index]
  t0 <- times[t0_index]

  # gamma0(s,t) = gamma( (s,t) - (s0,t0) )
  gamma0_vec <- gamma_st(
    dx = rep(coords$Longitude, times = lt) - rep(x0, n_sites*lt),
    dy = rep(coords$Latitude,  times = lt) - rep(y0, n_sites*lt),
    dt = rep(times, each = n_sites) - rep(t0, n_sites*lt),
    beta1, beta2, alpha1, alpha2, adv
  )
  gamma0 <- matrix(gamma0_vec, nrow = n_sites, ncol = lt, byrow = FALSE)

  W0 <- W[s0_index, t0_index]

  Y <- exp(W - W0 - gamma0)
  R <- 1 / runif(1)  # Pareto(1)
  Z <- threshold * R * Y
  rownames(Z) <- rownames(coords)

  return(list(Z = Z, W = W, gamma0 = gamma0))
}


bootstrap_ci <- function(x, B = 500, probs = c(0.025, 0.975)) {
  boot_means <- replicate(
    B,
    mean(sample(x, replace = TRUE), na.rm = TRUE)
  )
  quantile(boot_means, probs = probs, na.rm = TRUE)
}

