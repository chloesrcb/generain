
plot_th_emp_chi <- function(list_lags_filtered,
                            list_excesses_filtered,
                            wind_df_filtered,
                            params_estimates,
                            tau_min = 0,
                            tau_fixed = 0,
                            h_breaks = seq(0, 10, by = 1),
                            latlon = FALSE,
                            adv_transform = TRUE) {
  stopifnot(
    length(list_lags_filtered) == length(list_excesses_filtered),
    length(list_lags_filtered) <= nrow(wind_df_filtered)
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

  res_list <- vector("list", length(list_lags_filtered))

  for (i in seq_along(list_lags_filtered)) {
    lags_i    <- list_lags_filtered[[i]]
    excess_i  <- list_excesses_filtered[[i]]
    adv_x     <- wind_df_filtered$vx[i]
    adv_y     <- wind_df_filtered$vy[i]

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

    chi_th_i <- theoretical_chi(
      params   = params,
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

params_estimates <- params_est

# params_estimates <- param_estim_kmh
wind_df_filtered <- V_episodes_filtered
colnames(wind_df_filtered) <- c("vx", "vy")
df_p = plot_th_emp_chi(list_lags_filtered,
                list_excesses_filtered,
                wind_df_filtered,
                params_estimates,
                tau_min = 0,
                tau_fixed = 0,
                h_breaks = seq(0, 2, by = 0.2),
                latlon = FALSE,
                adv_transform = TRUE)

print(df_p$plots$all)
res_cmp <- df_p$res_cmp

p <- ggplot(res_cmp, aes(x = chi_theo_bar, y = chi_emp_bar)) +
  geom_jitter(width = 0, height = 0.005, alpha = 0.6) +
  geom_abline(linetype = "dashed") +
  theme_minimal()
# cor(res_cmp$chi_emp_bar, res_cmp$chi_theo_bar, method = "spearman")
ggplot(res_cmp, aes(x = chi_theo_bar, y = chi_emp_bar, color = factor(tau))) +
  geom_point(alpha = 0.7) +
  geom_abline(linetype = "dashed", color = "black") +
  theme_bw() +
  labs(color = "tau")

cor(res_cmp$chi_emp_bar, res_cmp$chi_theo_bar, method = "spearman")
p <- ggplot(res_cmp, aes(
      x = chi_theo_bar,
      y = chi_emp_bar,
      size = n_pairs
  )) +
  labs(
      x = "Theoretical Chi",
      y = "Empirical Chi"
    ) +
  geom_abline(linetype = "dashed", color = "red") +
  geom_point(alpha = 0.8, color = btfgreen) +
  theme_minimal() +
  btf_theme
p

# save plot
foldername <- paste0(im_folder, "/optim/omsev/th_vs_emp_chi/")
if (!dir.exists(foldername)) {
  dir.create(foldername, recursive = TRUE)
}

filename <- paste0(
  foldername,
  "th_vs_emp_chi_alltaus_npairs_adv_still03.png"
)
ggsave(filename,
       plot = p,
       width = 6,
       height = 5)


# facet by tau
p_facet <- ggplot(res_cmp, aes(
  x = chi_theo_bar,
  y = chi_emp_bar
)) +
  geom_point(alpha = 0.6, color = btfgreen) +
  geom_abline(linetype = "dashed", color = "red") +
  theme_bw() +
  facet_wrap(~ tau) +
  labs(title = "Empirical vs Theoretical Chi by tau")
p_facet 

