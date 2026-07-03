library(dplyr)
library(ggplot2)
library(ggrepel)
library(latex2exp)


# omsev_params <- c(1.101, 3.716, 0.103, 0.676, 3.896299, 2.221052)
omsev_params <- c(1.076, 3.911, 0.116, 0.654)

dist_mat <- get_dist_mat(location_gauges) / 1000
df_dist <- reshape_distances(dist_mat)

h_all <- df_dist$value
h_pos <- h_all[is.finite(h_all) & h_all > 0]
sort(unique(h_pos))

make_h_breaks <- function(h, n_bins = 6, min_width = 0.05) {
  h_pos <- h[is.finite(h) & h > 0]

  br <- quantile(
    h_pos,
    probs = seq(0, 1, length.out = n_bins + 1),
    na.rm = TRUE
  )

  br <- unique(as.numeric(br))

  # Force coverage
  br[1] <- 0
  br[length(br)] <- max(h_pos, na.rm = TRUE) + 0.1

  # Remove bins that are too narrow
  repeat {
    widths <- diff(br)
    if (all(widths >= min_width) || length(br) <= 3) break
    bad <- which(widths < min_width)[1]
    br <- br[-(bad + 1)]
  }

  br
}

h_breaks_omsev <- make_h_breaks(h_all, n_bins = 12, min_width = 0.05)

print(h_breaks_omsev)
h_all <- df_dist$value
h_pos <- h_all[is.finite(h_all) & h_all > 0]

print(table(
  cut(h_pos, breaks = h_breaks_omsev, include.lowest = FALSE),
  useNA = "ifany"
))

#  Compute empirical and theoretical chi

chi_omsev <- compute_group_chi(
  list_lags = list_lags_filtered,
  list_excesses = list_excesses_filtered,
  wind_df = V_episodes_filtered,
  params = omsev_params,
  tau_fixed = 0,
  h_breaks = h_breaks_omsev,
  adv_transform = TRUE
)

res_om_cmp <- chi_omsev$res_cmp


res_chi <- res_om_cmp %>%
  mutate(
    diff = chi_emp_bar - chi_theo_bar,
    abs_diff = abs(diff),
    rel_diff = diff / chi_theo_bar,
    tau_lab = paste0(round(tau * 60), " min"),
    hbin_chr = as.character(hbin),
    h_type = ifelse(
      grepl("^\\[0", hbin_chr) | grepl("^\\(0", hbin_chr),
      "Short spatial lags",
      "Positive spatial lags"
    ),
    direction = case_when(
      diff > 0 ~ "empirical > theoretical",
      diff < 0 ~ "empirical < theoretical",
      TRUE ~ "on diagonal"
    )
  )


corr_omsev_all <- summary_correlation(res_chi)
message(sprintf("OMSEV chi Spearman correlation, all bins: %.3f", corr_omsev_all))


p_chi_scatter <- ggplot(
  res_chi,
  aes(x = chi_theo_bar, y = chi_emp_bar, size = n_pairs)
) +
  geom_point(alpha = 0.7, color = btfgreen) +
  geom_abline(slope = 1, intercept = 0,
              linetype = "dashed", color = "red") +
  scale_size_continuous(range = c(1.5, 4)) +
  labs(
    x = TeX("Theoretical $\\chi$"),
    y = TeX("Empirical $\\chi$"),
    size = "Number of pairs"
  ) +
  btf_theme +
  theme(
    legend.position = "bottom",
    legend.title = element_text(size = 16),
    legend.text = element_text(size = 14)
  )

print(p_chi_scatter)


p_chi_scatter_tau <- ggplot(
  res_chi,
  aes(x = chi_theo_bar, y = chi_emp_bar,
      color = factor(tau_lab), size = n_pairs)
) +
  geom_point(alpha = 0.75) +
  geom_abline(slope = 1, intercept = 0,
              linetype = "dashed", color = "red") +
  scale_size_continuous(range = c(1.5, 4)) +
  labs(
    x = TeX("Theoretical $\\chi$"),
    y = TeX("Empirical $\\chi$"),
    color = "Temporal lag",
    size = "Number of pairs"
  ) +
  btf_theme +
  theme(legend.position = "bottom",
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 14))

print(p_chi_scatter_tau)

# Tau > 0 only

res_taupos <- res_chi %>%
  filter(tau > 0)

p_chi_taupos <- ggplot(
  res_taupos,
  aes(x = chi_theo_bar, y = chi_emp_bar, size = n_pairs)
) +
  geom_point(alpha = 0.7, color = btfgreen) +
  geom_abline(slope = 1, intercept = 0,
              linetype = "dashed", color = "red") +
  scale_size_continuous(range = c(1.5, 4)) +
  labs(
    x = TeX("Theoretical $\\chi$"),
    y = TeX("Empirical $\\chi$"),
    size = "Number of pairs"
  ) +
  btf_theme +
  theme(legend.position = "bottom",
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 14))

print(p_chi_taupos)

# -------------------------------------------------------------------------
# Tau = 0 only
# -------------------------------------------------------------------------

res_tau0 <- res_chi %>%
  filter(tau == 0)

p_chi_tau0 <- ggplot(
  res_tau0,
  aes(x = chi_theo_bar, y = chi_emp_bar, size = n_pairs)
) +
  geom_point(alpha = 0.7, color = btfgreen) +
  geom_abline(slope = 1, intercept = 0,
              linetype = "dashed", color = "red") +
  geom_text_repel(
    aes(label = hbin),
    size = 3,
    max.overlaps = 20
  ) +
  scale_size_continuous(range = c(1.5, 4)) +
  labs(
    x = TeX("Theoretical $\\chi$"),
    y = TeX("Empirical $\\chi$"),
    size = "Number of pairs"
  ) +
  btf_theme +
  theme(legend.position = "bottom",
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 14))

print(p_chi_tau0)

# -------------------------------------------------------------------------
# 9. Residuals as a function of theoretical chi
# -------------------------------------------------------------------------

p_chi_residual <- ggplot(
  res_chi,
  aes(x = chi_theo_bar, y = diff,
      color = factor(tau_lab), size = n_pairs)
) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  geom_point(alpha = 0.75) +
  scale_size_continuous(range = c(1.5, 4)) +
  labs(
    x = TeX("Theoretical $\\chi$"),
    y = TeX("$\\chi_{emp} - \\chi_{theo}$"),
    color = "Temporal lag",
    size = "Number of pairs"
  ) +
  btf_theme +
  theme(legend.position = "right")

print(p_chi_residual)

# -------------------------------------------------------------------------
# 10. Residuals by spatial distance bin
# -------------------------------------------------------------------------

p_chi_residual_hbin <- ggplot(
  res_chi,
  aes(x = hbin, y = diff, fill = factor(tau_lab))
) +
  geom_col(position = "dodge") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  labs(
    x = "Spatial distance bin",
    y = TeX("$\\chi_{emp} - \\chi_{theo}$"),
    fill = "Temporal lag"
  ) +
  btf_theme +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "right"
  ) +
  coord_cartesian(ylim = c(-0.25, 0.25))

print(p_chi_residual_hbin)

# -------------------------------------------------------------------------
# 11. Chi as function of spatial distance, one panel per temporal lag
# -------------------------------------------------------------------------
get_h_mid <- function(hbin) {
  x <- as.character(hbin)

  sapply(x, function(z) {
    if (z == "h = 0") return(0)

    nums <- as.numeric(
      unlist(regmatches(z, gregexpr("[0-9.]+", z)))
    )

    mean(nums, na.rm = TRUE)
  })
}

res_chi <- res_chi %>%
  mutate(
    h_mid = get_h_mid(hbin)
  )
# reorder time lags for plotting
res_chi$tau_lab <- factor(res_chi$tau_lab,
                              levels = unique(res_chi$tau_lab[order(res_chi$tau)]))

res_chi <- res_chi %>%
  arrange(tau, h_mid) %>%
  group_by(tau) %>%
  mutate(
    chi_emp_smooth = as.numeric(stats::filter(
      chi_emp_bar,
      filter = c(1/3, 1/3, 1/3),
      sides = 2
    ))
  ) %>%
  ungroup()

# remove tau > 0.6
res_chi_plot <- res_chi %>%
  filter(tau <= 0.7)

p_chi_by_h_tau_smooth <- ggplot(
  res_chi_plot,
  aes(x = h_mid)
) +
  geom_point(aes(y = chi_emp_bar, color = "Empirical"),
             size = 1.7, alpha = 0.45) +
  geom_smooth(aes(y = chi_emp_bar, color = "Empirical smoothed"),
              method = "loess", se = FALSE,
              span = 0.5, linewidth = 0.9) +
  geom_line(aes(y = chi_theo_bar, color = "Theoretical"),
            linewidth = 0.9, linetype = "dashed") +
  facet_wrap(~ tau_lab) +
  labs(
    x = "Spatial distance h (km)",
    y = TeX("$\\chi$"),
    color = ""
  ) +
  btf_theme +
  theme(legend.position = "bottom")
print(p_chi_by_h_tau_smooth)


res_chi_ci <- res_chi %>%
  mutate(
    hbin_id = as.numeric(hbin),
    hbin_mid = case_when(
      hbin == "h = 0" ~ 0,
      TRUE ~ sapply(as.character(hbin), function(x) {
        nums <- as.numeric(unlist(regmatches(x, gregexpr("[0-9.]+", x))))
        mean(nums)
      })
    )
  )

x_labs <- res_chi_ci %>%
  distinct(hbin_id, hbin_mid) %>%
  arrange(hbin_id) %>%
  mutate(
    label = ifelse(row_number() %% 2 == 1, sprintf("%.2f", hbin_mid), "")
  )

# -------------------------------------------------------------------------
# 12. Heatmap of signed residuals
# -------------------------------------------------------------------------

p_chi_heatmap <- ggplot(
  res_chi_plot,
  aes(x = hbin, y = factor(tau_lab), fill = diff)
) +
  geom_tile(color = "white") +
  scale_fill_gradient2(
    low = "#9d52b4",
    mid = "white",
    high = "#d18181",
    midpoint = 0,
    limits = c(-0.25, 0.25),
    name = TeX("$\\chi_{emp} - \\chi_{theo}$")
  ) +
  geom_text(aes(label = round(diff, 2)), size = 3) +
  labs(
    x = "Spatial distance bin",
    y = "Temporal lag"
  ) +
  btf_theme +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "right"
  )

print(p_chi_heatmap)

# -------------------------------------------------------------------------
# 13. Save plots
# -------------------------------------------------------------------------

foldername <- paste0(im_folder, "optim/omsev/2025_results/")
if (!dir.exists(foldername)) {
  dir.create(foldername, recursive = TRUE)
}

ggsave(paste0(foldername, "omsev_chi_scatter.png"),
       p_chi_scatter, width = 7, height = 7)

ggsave(paste0(foldername, "omsev_chi_scatter_by_tau.png"),
       p_chi_scatter_tau, width = 7, height = 7)

ggsave(paste0(foldername, "omsev_chi_scatter_taupos.png"),
       p_chi_taupos, width = 7, height = 7)

ggsave(paste0(foldername, "omsev_chi_scatter_tau0.png"),
       p_chi_tau0, width = 7, height = 7)

ggsave(paste0(foldername, "omsev_chi_residuals.png"),
       p_chi_residual, width = 7, height = 6)

ggsave(paste0(foldername, "omsev_chi_residuals_by_hbin.png"),
       p_chi_residual_hbin, width = 9, height = 6)

ggsave(paste0(foldername, "omsev_chi_by_hbin_tau.png"),
       p_chi_by_h_tau, width = 10, height = 7)

ggsave(paste0(foldername, "omsev_chi_residual_heatmap.png"),
       p_chi_heatmap, width = 9, height = 6)

ggsave(paste0(foldername, "omsev_chi_by_hbin_tau_smooth.png"),
       p_chi_by_h_tau_smooth, width = 10, height = 7)

foldername_jk <- paste0(data_folder, "omsev/optim_results/jackknife_estimates/")

n_eff <-71
filename <- paste0(foldername_jk, "all_results_jk_by_monthyear_n",
           n_eff, "_q", q*100, "_delta", delta, "_dmin", min_spatial_dist, ".csv")

jack_estimates <- read.csv(filename, header = TRUE)
colnames(jack_estimates) <- c("beta1", "beta2", "alpha1", "alpha2")
       # jack_estimates must contain beta1, beta2, alpha1, alpha2
# If eta fixed, add eta1 eta2 if compute_group_chi expects 6 params
eta_fixed <- c(eta1 = 3.896299, eta2 = 2.221052)


chi_jk_list <- lapply(seq_len(nrow(jack_estimates)), function(i) {

  params_i <- as.numeric(c(
    jack_estimates[i, c("beta1", "beta2", "alpha1", "alpha2")],
    eta_fixed
  ))

  names(params_i) <- c("beta1", "beta2", "alpha1", "alpha2", "eta1", "eta2")

  chi_i <- compute_group_chi(
    list_lags = list_lags_filtered,
    list_excesses = list_excesses_filtered,
    wind_df = V_episodes_filtered,
    params = params_i,
    tau_fixed = 0,
    h_breaks = h_breaks_omsev,
    adv_transform = TRUE
  )

  chi_i$res_cmp %>%
    mutate(jk_id = i)
})

chi_jk_all <- bind_rows(chi_jk_list)
chi_band <- chi_jk_all %>%
  group_by(tau, hbin) %>%
  summarise(
    chi_dot = mean(chi_theo_bar, na.rm = TRUE),
    se_chi = sqrt(
      ((n() - 1) / n()) *
        sum((chi_theo_bar - mean(chi_theo_bar, na.rm = TRUE))^2, na.rm = TRUE)
    ),
    .groups = "drop"
  )

res_chi_ci <- res_chi %>%
  left_join(chi_band, by = c("tau", "hbin")) %>%
  mutate(
    chi_theo_low  = pmax(0, chi_theo_bar - qnorm(0.975) * se_chi),
    chi_theo_high = pmin(1, chi_theo_bar + qnorm(0.975) * se_chi)
  )

ggplot(res_chi_ci,
       aes(x = chi_theo_bar, y = chi_emp_bar, size = n_pairs)) +
  geom_errorbarh(
    aes(xmin = chi_theo_low, xmax = chi_theo_high),
    height = 0,
    alpha = 0.4, color = btfgreen
  ) +
  geom_abline(slope = 1, intercept = 0,
              linetype = "dashed", color = "red") +
  scale_size_continuous(range = c(1.5, 4)) +
  labs(
    x = TeX("Theoretical $\\chi$"),
    y = TeX("Empirical $\\chi$"),
    size = "Number of pairs"
  ) +
  btf_theme +
  ylim(0, 0.6) +
  xlim(0, 0.6) +
  theme(
    legend.position = "bottom",
    legend.title = element_text(size = 16),
    legend.text = element_text(size = 14)
  )

# save 
foldername <- paste0(im_folder, "optim/omsev/2025_results/")
if (!dir.exists(foldername)) {
  dir.create(foldername, recursive = TRUE)
}
filename <- paste0(foldername, "omsev_chi_band_scatter.png")
ggsave(filename, width = 7, height = 7)


p_chi_band_by_h_tau <- ggplot(
  res_chi_ci,
  aes(x = hbin)
) +
  geom_ribbon(
    aes(
      ymin = chi_theo_low,
      ymax = chi_theo_high,
      group = tau_lab
    ),
    alpha = 0.25
  ) +
  geom_line(
    aes(y = chi_theo_bar, group = tau_lab),
    linetype = "dashed",
    linewidth = 0.9
  ) +
  geom_point(
    aes(y = chi_emp_bar, size = n_pairs),
    color = btfgreen,
    alpha = 0.8
  ) +
  facet_wrap(~ tau_lab) +
  scale_size_continuous(range = c(1.5, 4)) +
  labs(
    x = "Spatial distance bin",
    y = TeX("$\\chi$"),
    size = "Number of pairs"
  ) +
  btf_theme +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "bottom"
  )

print(p_chi_band_by_h_tau)


get_h_mid <- function(hbin) {
  x <- as.character(hbin)
  sapply(x, function(z) {
    if (z == "h = 0") return(0)
    nums <- as.numeric(unlist(regmatches(z, gregexpr("[0-9.]+", z))))
    mean(nums, na.rm = TRUE)
  })
}

res_chi_ci <- res_chi_ci %>%
  mutate(
    h_mid = get_h_mid(hbin),
    tau_lab = factor(tau_lab, levels = paste0(sort(unique(tau)) * 60, " min"))
  ) %>%
  arrange(tau, h_mid)

p_chi_band_by_h_tau <- ggplot(
  res_chi_ci,
  aes(x = h_mid)
) +
  geom_ribbon(
    aes(
      ymin = chi_theo_low,
      ymax = chi_theo_high,
      group = tau_lab
    ),
    alpha = 0.25,
    fill = "grey70"
  ) +
  geom_line(
    aes(y = chi_theo_bar),
    linetype = "dashed",
    linewidth = 0.9
  ) +
  geom_point(
    aes(y = chi_emp_bar, size = n_pairs),
    color = btfgreen,
    alpha = 0.8
  ) +
  facet_wrap(~ tau_lab) +
  scale_size_continuous(range = c(0.5, 2)) +
  labs(
    x = "Spatial distance h (km)",
    y = TeX("$\\chi$"),
    size = "Number of pairs"
  ) +
  btf_theme +
  theme(
    legend.position = "bottom"
  )

print(p_chi_band_by_h_tau)

# remove tau > 0.7 for better visualization
res_chi_ci_plot <- res_chi_ci %>%
  filter(tau <= 0.6)

p_chi_band_by_h_tau <- ggplot(
  res_chi_ci_plot,
  aes(x = h_mid)
) +
  geom_ribbon(
    aes(
      ymin = chi_theo_low,
      ymax = chi_theo_high,
      group = tau_lab
    ),
    alpha = 0.25,
    fill = "grey70"
  ) +
  geom_line(
    aes(y = chi_theo_bar),
    linetype = "dashed",
    linewidth = 0.9
  ) +
  geom_point(
    aes(y = chi_emp_bar, size = n_pairs),
    color = btfgreen,
    alpha = 0.8
  ) +
  facet_wrap(~ tau_lab, ncol = 4) +
  scale_size_continuous(range = c(0.5, 2)) +
  labs(
    x = TeX("Spatial distance $h$ (km)"),
    y = TeX("$\\widehat{\\chi}(h,\\tau)$"),
    size = "Number of pairs"
  ) +
  btf_theme +
  theme(
    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.background = element_rect(fill = "white", color = NA),
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 9),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

print(p_chi_band_by_h_tau)


# save
foldername <- paste0(im_folder, "optim/omsev/2025_results/")
if (!dir.exists(foldername)) {
  dir.create(foldername, recursive = TRUE)
}

ggsave(paste0(foldername, "omsev_chi_band_by_h_tau.pdf"),
       p_chi_band_by_h_tau, width = 10, height = 10)


# same but remove h=0
res_chi_ci_no_h0 <- res_chi_ci_plot %>% filter(h_mid > 0)

ggplot(
  res_chi_ci_no_h0,
  aes(x = h_mid)
) +
  geom_ribbon(
    aes(
      ymin = chi_theo_low,
      ymax = chi_theo_high,
      group = tau_lab
    ),
    alpha = 0.25,
    fill = "grey70"
  ) +
  geom_line(
    aes(y = chi_theo_bar),
    linetype = "dashed",
    linewidth = 0.9
  ) +
  geom_point(
    aes(y = chi_emp_bar, size = n_pairs),
    color = btfgreen,
    alpha = 0.8
  ) +
  facet_wrap(~ tau_lab, ncol = 4) +
  scale_size_continuous(range = c(0.5, 2)) +
  labs(
    x = TeX("Spatial distance $h$ (km)"),
    y = TeX("$\\widehat{\\chi}(h,\\tau)$"),
    size = "Number of pairs"
  ) +
  btf_theme +
  theme(
    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.background = element_rect(fill = "white", color = NA),
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 9),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )


# save
foldername <- paste0(im_folder, "optim/omsev/2025_results/")
if (!dir.exists(foldername)) {
  dir.create(foldername, recursive = TRUE)
}

ggsave(paste0(foldername, "omsev_chi_band_by_h_tau_no_h0.pdf"),
       width = 10, height = 10)


# do the same on comephore data

com_params <- c(0.1579922, 0.99232, 0.5800247, 0.6627991, 3.896299, 2.221052)

# # get distance matrix and breaks for chi estimation
df_dist <- reshape_distances(dist_mat)
# round distances to 1 decimal place
df_dist$value <- round(df_dist$value, 0)
n_hbins <- 13
h_all <- df_dist$value

h_breaks_com <- quantile(
  h_all,
  probs = seq(0, 1, length.out = n_hbins + 1),
  na.rm = TRUE
)
h_breaks_com <- seq(0, 13, by=1)

length(h_breaks_com)
h_breaks_com

h_breaks_com <- unique(as.numeric(h_breaks_com))
h_breaks_com[length(h_breaks_com)] <- 13



chi_com <- compute_group_chi(
  list_lags = list_lags,
  list_excesses = list_excesses,
  wind_df = V_adv,
  params = com_params,
  tau_fixed = 0,
  h_breaks = h_breaks_com,
  adv_transform = TRUE
)

corr_com <- summary_correlation(chi_com$res_cmp)
message(sprintf("COM chi Spearman correlation: %.3f", corr_com))

# plot chi_come emp vs theo
print(chi_com$plots$all)
res_com_cmp <- chi_com$res_cmp

ggplot(res_com_cmp, aes(x = chi_theo_bar, y = chi_emp_bar, size = n_pairs)) +
  geom_point(alpha = 0.7, color = btfgreen) +
  geom_abline(linetype = "dashed", color = "red") +
  theme_bw() +
  scale_size_continuous(range = c(1, 4)) +
  labs(x = TeX("Theoretical $\\chi$"), y = TeX("Empirical $\\chi$"),
      size = "Number of pairs") +
  btf_theme +
  theme(legend.position = "bottom",
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 14))
# save plot
foldername <- paste0(im_folder, "workflows/full_pipeline/")
if (!dir.exists(foldername)) {
  dir.create(foldername, recursive = TRUE)
}
filename <- paste0(foldername, "com_chi_th_emp.png")
ggsave(filename,
       width = 7,
       height = 7)



# same plot as before with omsev

# without tau = 0

res_com_cmp_taupos <- res_com_cmp %>%
  filter(tau > 0)

ggplot(res_com_cmp_taupos, aes(x = chi_theo_bar, y = chi_emp_bar, size = n_pairs)) +
  geom_point(alpha = 0.7, color = btfgreen) +
  geom_abline(linetype = "dashed", color = "red") +
  theme_bw() +
  scale_size_continuous(range = c(1, 4)) +
  labs(x = TeX("Theoretical $\\chi$"), y = TeX("Empirical $\\chi$"),
       size = "Number of pairs") +
  btf_theme +
  theme(legend.position = "bottom",
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 14))


# heatmap of residuals for com

diff_com <- res_com_cmp %>%
  mutate(
    diff = chi_emp_bar - chi_theo_bar,
    tau_lab = paste0(tau, " h"),
    hbin_chr = as.character(hbin)
  )

# reorder time lags for plotting
diff_com$tau_lab <- factor(diff_com$tau_lab,
                           levels = unique(diff_com$tau_lab[order(diff_com$tau)]))

# keep only tau <= 6
diff_com <- diff_com %>%
  filter(tau <= 6)

# with jk estimates 

filename_jk <- "/home/cserreco/Documents/These/phd_extremes/data/comephore/optim_results/jackknife_estimates/jackknife_annual_estimates/all_annual_jk_q95_delta24_dmin5.csv"
df_jk_all <- read.csv(filename_jk, sep = ",")
jack_estimates <- df_jk_all %>%
  select(beta1, beta2, alpha1, alpha2, eta1, eta2)


chi_jk_list <- lapply(seq_len(nrow(jack_estimates)), function(i) {

  params_i <- as.numeric(jack_estimates[i, c("beta1", "beta2", "alpha1", "alpha2", "eta1", "eta2")])

  names(params_i) <- c("beta1", "beta2", "alpha1", "alpha2", "eta1", "eta2")

  chi_i <- compute_group_chi(
    list_lags = list_lags,
    list_excesses = list_excesses,
    wind_df = V_adv,
    params = params_i,
    tau_fixed = 0,
    h_breaks = h_breaks_com,
    adv_transform = TRUE
  )

  chi_i$res_cmp %>%
    mutate(jk_id = i)
})

chi_jk_all <- bind_rows(chi_jk_list)

# do scatter plot with band
chi_band <- chi_jk_all %>%
  group_by(tau, hbin) %>%
  summarise(
    chi_dot = mean(chi_theo_bar, na.rm = TRUE),
    se_chi = sqrt(
      ((n() - 1) / n()) *
        sum((chi_theo_bar - mean(chi_theo_bar, na.rm = TRUE))^2, na.rm = TRUE)
    ),
    .groups = "drop"
  )


# plot
res_chi <- res_com_cmp %>%
  mutate(
    diff = chi_emp_bar - chi_theo_bar,
    abs_diff = abs(diff),
    rel_diff = diff / chi_theo_bar,
    tau_lab = paste0(tau, " hours"),
    hbin_chr = as.character(hbin),
    h_type = ifelse(
      grepl("^\\[0", hbin_chr) | grepl("^\\(0", hbin_chr),
      "Short spatial lags",
      "Positive spatial lags"
    ),
    direction = case_when(
      diff > 0 ~ "empirical > theoretical",
      diff < 0 ~ "empirical < theoretical",
      TRUE ~ "on diagonal"
    )
  )

res_chi_ci <- res_chi %>%
  left_join(chi_band, by = c("tau", "hbin")) %>%
  mutate(
    chi_theo_low  = pmax(0, chi_theo_bar - qnorm(0.975) * se_chi),
    chi_theo_high = pmin(1, chi_theo_bar + qnorm(0.975) * se_chi)
  )

ggplot(res_chi_ci,
       aes(x = chi_theo_bar, y = chi_emp_bar, size = n_pairs)) +
  geom_errorbarh(
    aes(xmin = chi_theo_low, xmax = chi_theo_high),
    height = 0,
    alpha = 0.4, color = btfgreen
  ) +
  geom_abline(slope = 1, intercept = 0,
              linetype = "dashed", color = "red") +
  scale_size_continuous(range = c(1.5, 4)) +
  labs(
    x = TeX("Theoretical $\\chi$"),
    y = TeX("Empirical $\\chi$"),
    size = "Number of pairs"
  ) +
  btf_theme +
  ylim(0, 0.6) +
  xlim(0, 0.6) +
  theme(
    legend.position = "bottom",
    legend.title = element_text(size = 16),
    legend.text = element_text(size = 14)
  )

# save plot 
filename <- paste0(foldername, "com_chi_band_scatter.png")
ggsave(filename, width = 7, height = 7)

# reorder time lags for plotting
res_chi_ci$tau_lab <- factor(res_chi_ci$tau_lab,
                             levels = unique(res_chi_ci$tau_lab[order(res_chi_ci$tau)]))

# keep only tau <= 7
res_chi_ci <- res_chi_ci %>%
  filter(tau <= 7)
p_chi_band_by_h_tau <- ggplot(
  res_chi_ci,
  aes(x = hbin)
) +
  geom_ribbon(
    aes(
      ymin = chi_theo_low,
      ymax = chi_theo_high,
      group = tau_lab
    ),
    alpha = 0.1
  ) +
  geom_line(
    aes(y = chi_theo_bar, group = tau_lab),
    linetype = "dashed",
    linewidth = 0.9
  ) +
  geom_point(
    aes(y = chi_emp_bar, size = n_pairs),
    color = btfgreen,
    alpha = 0.8
  ) +
  facet_wrap(~ tau_lab, nrow = 2, ncol = 4) +
  scale_size_continuous(range = c(1.5, 4)) +
  scale_x_discrete(
  labels = function(x) {
    labs <- sapply(x, function(z) {
      nums <- as.numeric(unlist(regmatches(z, gregexpr("[0-9.]+", z))))
      if (length(nums) == 0) return(z)
      sprintf("%.1f", mean(nums))
    })

    labs[seq_along(labs) %% 2 == 0] <- ""
    labs
  }
) +
  labs(
    x = "Spatial distance bin",
    y = expression(hat(chi)(bold(h), tau)),
    size = "Number of pairs"
  ) +
  btf_theme +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
    legend.position = "bottom"
  )

print(p_chi_band_by_h_tau)

# save
foldername <- paste0(im_folder, "optim/comephore/2025_results/")
if (!dir.exists(foldername)) {
  dir.create(foldername, recursive = TRUE)
}

ggsave(paste0(foldername, "comephore_chi_band_by_h_tau.pdf"),
       p_chi_band_by_h_tau, width = 10, height = 10)
