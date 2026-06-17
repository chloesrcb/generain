library(dplyr)
library(tidyr)
library(ggplot2)
library(tibble)
library(purrr)

M <- 100
p_tip <- 0.2152
timesteps <- 0:11

sites_names <- colnames(rain)
grid_omsev <- grid_coords_m

folder_margins <- paste0(im_folder, "swg/omsev/2025/margins/above0_disc/")
folder_cumuls  <- paste0(im_folder, "swg/omsev/2025/cumuls/")

if (!dir.exists(folder_margins)) dir.create(folder_margins, recursive = TRUE)
if (!dir.exists(folder_cumuls))  dir.create(folder_cumuls, recursive = TRUE)

apply_discretization <- function(x, p = 0.2152) {
#   x[!is.finite(x)] <- NA_real_
  x[x < 1e-2] <- 0
  x[x > 0 & x < p] <- p
  x[x > p & x < 2 * p] <- 2 * p
  x
}

# Prepare advection and filter valid episodes

adv_matrix <- as.matrix(adv_df_transfo[, c("vx_t", "vy_t")])

# km/h -> m/5min
adv_matrix <- adv_matrix * (1000 * 5 / 60)

id_not_invalid <- which(adv_class$group != "invalid")

adv_matrix_filtered    <- adv_matrix[id_not_invalid, , drop = FALSE]
list_episodes_filtered <- list_episodes[id_not_invalid]
s0_list_filtered       <- s0_list[id_not_invalid]
u_list_filtered        <- u_list[id_not_invalid]

Nep <- length(list_episodes_filtered)

# Simulate M repetitions for each episode

set.seed(123)

sims_by_ep <- vector("list", Nep)
u_sim <- numeric(Nep)

for (j in seq_len(Nep)) {
  s0_j  <- s0_list_filtered[j]
  u_j   <- u_list_filtered[j]
  adv_j <- adv_matrix_filtered[j, ]

  u_sim[j] <- u_j
  sims_by_ep[[j]] <- vector("list", M)

  for (m in seq_len(M)) {
    sims_by_ep[[j]][[m]] <- simulate_many_episodes(
      N = 1,
      u_emp = u_j,
      params_vario = params_m5min,
      params_margins = params_margins,
      coords = grid_omsev,
      times = timesteps,
      adv = adv_j,
      t0 = 0,
      s0 = s0_j
    )[[1]]
  }
}


sims_all_raw <- unlist(sims_by_ep, recursive = FALSE)
sims_all_disc <- sims_all_raw

for (k in seq_along(sims_all_disc)) {
  sims_all_disc[[k]]$X <- apply_discretization(
    sims_all_disc[[k]]$X,
    p = p_tip
  )
}
# Marginal plots by site: density + boxplot
for (site in sites_names) {

  Xobs_site <- unlist(
    lapply(list_episodes_filtered,
           function(ep) as.numeric(ep[, site]))
  )

  Xsim_raw <- unlist(
    lapply(sims_all_raw,
           function(sim) as.numeric(sim$X[, site]))
  )

  Xsim_disc <- unlist(
    lapply(sims_all_disc,
           function(sim) as.numeric(sim$X[, site]))
  )

  Xobs_site  <- Xobs_site [is.finite(Xobs_site)  & Xobs_site  > 0]
  Xsim_raw   <- Xsim_raw  [is.finite(Xsim_raw)   & Xsim_raw   > 1e-2]
  Xsim_disc  <- Xsim_disc [is.finite(Xsim_disc)  & Xsim_disc  > 1e-2]

  df_site <- bind_rows(
    tibble(source = "Observations", value = Xobs_site),
    tibble(source = "Simulations raw", value = Xsim_raw),
    tibble(source = "Simulations corrected", value = Xsim_disc)
  )

  color_legend <- c(
    "Observations" = "#ee8686",
    "Simulations raw" = "#86a8e7",
    "Simulations corrected" = "#86e7a8"
  )

  p_site <- ggplot(df_site, aes(x = value, fill = source)) +
    geom_density(
      alpha = 0.35,
      bw = 0.2
    ) +
    geom_boxplot(
      aes(y = -0.1, group = source),
      width = 0.15,
      alpha = 0.6,
      outlier.size = 0.5
    ) +
    labs(
      x = "Rainfall (mm / 5 min)",
      y = "Density",
      fill = "Source"
    ) +
    btf_theme +
    xlim(0,10)

    
  ggsave(
    paste0(folder_margins, "density_box_above0_disc_", site, ".png"),
    p_site,
    width = 8,
    height = 6
  )
}

# Episode cumulative rainfall

sims_all <- unlist(sims_by_ep, recursive = FALSE)

#  Match observation NA pattern and discretize simulations

for (j in seq_len(Nep)) {
  obs_mat <- as.matrix(list_episodes_filtered[[j]][, sites_names, drop = FALSE])

  for (m in seq_len(M)) {
    sim_mat <- as.matrix(sims_by_ep[[j]][[m]]$X[, sites_names, drop = FALSE])

    # same missingness as observations
    sim_mat[is.na(obs_mat)] <- NA

    # pseudo tipping-bucket measurement
    sim_mat <- apply_discretization(sim_mat, p = p_tip)

    sims_by_ep[[j]][[m]]$X[, sites_names] <- sim_mat
  }
}

cumul_obs <- numeric(Nep)
cumul_sim <- matrix(NA_real_, nrow = Nep, ncol = M)

for (j in seq_len(Nep)) {
  obs_mat <- as.matrix(list_episodes_filtered[[j]][, sites_names, drop = FALSE])
  cumul_obs[j] <- sum(obs_mat, na.rm = TRUE)

  for (m in seq_len(M)) {
    sim_mat <- as.matrix(sims_by_ep[[j]][[m]]$X[, sites_names, drop = FALSE])
    cumul_sim[j, m] <- sum(sim_mat, na.rm = TRUE)
  }
}

df_cumul <- bind_rows(
  tibble(
    ep_id = seq_len(Nep),
    rep = NA_integer_,
    source = "Observations",
    cumul = cumul_obs
  ),
  as.data.frame(cumul_sim) %>%
    mutate(ep_id = seq_len(n())) %>%
    pivot_longer(
      cols = -ep_id,
      names_to = "rep",
      values_to = "cumul"
    ) %>%
    mutate(
      rep = as.integer(gsub("V", "", rep)),
      source = "Simulations"
    ) %>%
    select(ep_id, rep, source, cumul)
) %>%
  filter(is.finite(cumul), cumul >= 0) %>%
  mutate(
    source = factor(source, levels = c("Observations", "Simulations")),
    log_cumul = log1p(cumul)
  )

# Cumul density + boxplot, raw scale


p_cumul_raw <- ggplot(df_cumul, aes(x = cumul, fill = source)) +
  geom_density(
    alpha = 0.35,
    bw = 40,
    linewidth = 0.4
  ) +
  geom_boxplot(
    aes(y = -0.0008, group = source),
    width = 0.0006,
    alpha = 0.6,
    outlier.size = 0.12,
    colour = "black"
  ) +
  labs(
    x = "Episode cumulative rainfall",
    y = "Density",
    fill = "Source"
  ) +
  btf_theme +
  xlim(0, 700)

ggsave(
  paste0(folder_cumuls, "cumul_density_box_raw.png"),
  p_cumul_raw,
  width = 10,
  height = 6
)

# Cumul density + boxplot, log scale

p_cumul_log <- ggplot(df_cumul, aes(x = log_cumul, fill = source)) +
  geom_density(
    alpha = 0.35,
    adjust = 1.1,
    linewidth = 0.4
  ) +
  geom_boxplot(
    aes(y = -0.08, group = source),
    width = 0.04,
    alpha = 0.6,
    outlier.size = 0.12,
    colour = "black"
  ) +
  labs(
    x = "log(1 + episode cumulative rainfall)",
    y = "Density",
    fill = "Source"
  ) +
  btf_theme

ggsave(
  paste0(folder_cumuls, "cumul_density_box_log.png"),
  p_cumul_log,
  width = 10,
  height = 6
)
