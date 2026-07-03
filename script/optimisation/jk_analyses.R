library(ggplot2)
library(ggrepel)
library(dplyr)
# Load libraries and set theme
source("./script/load_libraries.R")


filename_jk <- "/home/cserreco/Documents/These/phd_extremes/data/comephore/optim_results/jackknife_estimates/jackknife_annual_estimates/all_annual_jk_q95_delta24_dmin5.csv"
filename_jk <- "/home/cserreco/Documents/These/phd_extremes/data/omsev/optim_results/jackknife_estimates/all_results_jk_by_monthyear_n71_q95_delta12_dmin1200.csv"
df_jk_all <- read.csv(filename_jk, sep = ",")
# df_jk_all <- jack_estimates
param_names <- c("beta1", "beta2", "alpha1", "alpha2", "eta1", "eta2")
param_names <- c("beta1", "beta2", "alpha1", "alpha2")

theta_full <- c(
  beta1  = 0.1579947,
  beta2  = 0.9923211,
  alpha1 = 0.5800001,
  alpha2 = 0.6627909,
  eta1   = 3.896299,
  eta2   = 2.221052
)

theta_full <- c(
  beta1  = 1.076,
  beta2  = 3.911,
  alpha1 = 0.116,
  alpha2 = 0.654
)


#omsev
theta_full <- c(
  beta1  = init_param_jk[1],
  beta2  = init_param_jk[2],
  alpha1 = init_param_jk[3],
  alpha2 = init_param_jk[4]
)

colnames(df_jk_all) <- param_names
df_plot <- df_jk_all 

# df_plot$block <- unique_month_years

df_plot <- df_plot %>%
  mutate(
    is_outlier = beta1 < quantile(beta1, 0.05, na.rm = TRUE) |
                 beta1 > quantile(beta1, 0.95, na.rm = TRUE) |
                 alpha1 < quantile(alpha1, 0.05, na.rm = TRUE) |
                 alpha1 > quantile(alpha1, 0.95, na.rm = TRUE)
  )



ggplot(df_plot, aes(x = eta2, y = beta2)) +
  geom_point(color = "#1f77b4", size = 3, alpha = 0.8) +
  geom_smooth(
    method = "lm",
    formula = y ~ poly(x, 2),
    color = "darkred",
    se = FALSE,
    linetype = "dashed",
    linewidth = 0.8
  ) +
  geom_text_repel(
  data = df_plot,
  aes(label = block),
  size = 3.5,
  max.overlaps = Inf,
  box.padding = 0.4,
  point.padding = 0.3,
  segment.color = "grey50",
  segment.size = 0.4
) +
  geom_point(
    aes(x = theta_full["eta2"], y = theta_full["beta2"]),
    color = "red",
    size = 4,
    inherit.aes = FALSE
  ) +
  labs(
    x = expression(hat(eta)[2]),
    y = expression(hat(beta)[2])
  ) +
  btf_theme +
  theme(
    plot.subtitle = element_text(face = "italic", color = "grey40", size = 11),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "grey80")
  )

ggsave(
  paste0(im_folder, "/optim/comephore/jackknife_eta2_beta2_relationship_q",
         q * 100, "_delta", delta, "_dmin", min_spatial_dist, ".pdf"),
  width = 8,
  height = 6
)

df_jk_valid <- df_jk_all %>%
  filter(convergence == 0, !is.na(convergence))

df_jk_valid <- df_jk_all %>%
  filter(!(block %in% c("2024-03", "2025-05", "2024-05", "2022-08", "2022-09")))

G <- nrow(df_jk_valid)

param_names <- c("beta1", "beta2", "alpha1", "alpha2")
jack_estimates <- as.matrix(df_jk_valid[, param_names])
storage.mode(jack_estimates) <- "double"

log_jack_estimates <- log(jack_estimates)
log_theta_full <- log(theta_full)
log_theta_dot <- colMeans(log_jack_estimates)

jack_se_log <- sqrt(
  ((G - 1) / G) *
    colSums(
      (
        log_jack_estimates -
          matrix(log_theta_dot, G, length(param_names), byrow = TRUE)
      )^2
    )
)

theta_dot <- colMeans(jack_estimates)

jack_se_nat <- sqrt(
  ((G - 1) / G) *
    colSums(
      (
        jack_estimates -
          matrix(theta_dot, G, length(param_names), byrow = TRUE)
      )^2
    )
)

# Pseudo-values only as diagnostics

pseudo_values_log <- matrix(
  NA_real_,
  nrow = G,
  ncol = length(param_names),
  # dimnames = list(df_jk_valid$block, param_names)
)

pseudo_values_natural <- matrix(
  NA_real_,
  nrow = G,
  ncol = length(param_names),
  # dimnames = list(df_jk_valid$block, param_names)
)

for (i in seq_len(G)) {
  pseudo_values_log[i, ] <-
    G * log_theta_full - (G - 1) * log_jack_estimates[i, ]

  pseudo_values_natural[i, ] <-
    G * theta_full - (G - 1) * jack_estimates[i, ]
}

theta_jack_log_corr <- exp(colMeans(pseudo_values_log))
theta_jack_nat_corr <- colMeans(pseudo_values_natural)


z <- qnorm(0.975)

ci_log_lower <- exp(log_theta_full - z * jack_se_log)
ci_log_upper <- exp(log_theta_full + z * jack_se_log)

ci_nat_lower <- theta_full - z * jack_se_nat
ci_nat_upper <- theta_full + z * jack_se_nat

# ci_nat_lower <- pmax(ci_nat_lower, 1e-8)

jk_intervals_df <- data.frame(
  Parameter              = param_names,
  Estimate_Full          = as.numeric(theta_full),
  StdError_LogScale      = as.numeric(jack_se_log),
  CI_Lower_LogScale      = as.numeric(ci_log_lower),
  CI_Upper_LogScale      = as.numeric(ci_log_upper),
  StdError_NaturalScale  = as.numeric(jack_se_nat),
  CI_Lower_NaturalScale  = as.numeric(ci_nat_lower),
  CI_Upper_NaturalScale  = as.numeric(ci_nat_upper),
  Jack_Corr_LogScale     = as.numeric(theta_jack_log_corr),
  Jack_Corr_NaturalScale = as.numeric(theta_jack_nat_corr),
  n_blocks_effective     = G
)

print(jk_intervals_df)




# ============================================================
# Conversion km/h -> m/5min
# ============================================================

c_x_m <- 1000
c_t_5min <- 12

# Estimate full converted
b1_full_m5 <- theta_full["beta1"] / (c_x_m ^ theta_full["alpha1"])
b2_full_m5 <- theta_full["beta2"] / (c_t_5min ^ theta_full["alpha2"])

# Leave-one-block estimates converted
beta1_m5_i <- jack_estimates[, "beta1"] / (c_x_m ^ jack_estimates[, "alpha1"])
beta2_m5_i <- jack_estimates[, "beta2"] / (c_t_5min ^ jack_estimates[, "alpha2"])

# Jackknife SE directly on converted scale
G <- nrow(jack_estimates)

se_b1_m5 <- sqrt(
  ((G - 1) / G) * sum((beta1_m5_i - mean(beta1_m5_i))^2)
)

se_b2_m5 <- sqrt(
  ((G - 1) / G) * sum((beta2_m5_i - mean(beta2_m5_i))^2)
)

z <- qnorm(0.975)

ci_b1_m5 <- c(
  b1_full_m5 - z * se_b1_m5,
  b1_full_m5 + z * se_b1_m5
)

ci_b2_m5 <- c(
  b2_full_m5 - z * se_b2_m5,
  b2_full_m5 + z * se_b2_m5
)

ci_b1_m5[1] <- max(ci_b1_m5[1], 1e-8)
ci_b2_m5[1] <- max(ci_b2_m5[1], 1e-8)

# alpha unchanged by unit conversion
param_table_m5min <- data.frame(
  Parameter = c("beta1", "beta2", "alpha1", "alpha2"),
  Estimate_full_m5min = c(
    b1_full_m5,
    b2_full_m5,
    theta_full["alpha1"],
    theta_full["alpha2"]
  ),
  StdError_m5min = c(
    se_b1_m5,
    se_b2_m5,
    jk_intervals_df$StdError_NaturalScale[3],
    jk_intervals_df$StdError_NaturalScale[4]
  ),
  CI_lower_m5min = c(
    ci_b1_m5[1],
    ci_b2_m5[1],
    jk_intervals_df$CI_Lower_NaturalScale[3],
    jk_intervals_df$CI_Lower_NaturalScale[4]
  ),
  CI_upper_m5min = c(
    ci_b1_m5[2],
    ci_b2_m5[2],
    jk_intervals_df$CI_Upper_NaturalScale[3],
    jk_intervals_df$CI_Upper_NaturalScale[4]
  )
)

print(param_table_m5min)
