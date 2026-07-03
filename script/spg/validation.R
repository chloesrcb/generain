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



library(reshape2)
library(ggplot2)

# seuils globaux par site
u_by_s0_from_list <- tapply(u_list, s0_list, median, na.rm = TRUE)
u_global <- u_by_s0_from_list

sites <- colnames(list_episodes_filtered[[1]])

get_X_matrix <- function(obj) {
  if (is.null(obj)) return(NULL)
  if (is.matrix(obj) || is.data.frame(obj)) return(as.matrix(obj))
  if (is.list(obj) && !is.null(obj$X)) return(as.matrix(obj$X))
  NULL
}

apply_obs_na_mask <- function(Xsim, Xobs, sites) {
  Xsim2 <- Xsim
  mask_na <- is.na(as.matrix(Xobs)[, sites, drop = FALSE])
  Xsim2[mask_na] <- NA
  Xsim2
}

cond_probs_episode_q_na <- function(X, u_vec, sites = colnames(X)) {
  X <- as.matrix(X)[, sites, drop = FALSE]
  p <- length(sites)

  if (length(u_vec) == 1) u_vec <- rep(u_vec, p)
  if (is.null(names(u_vec))) names(u_vec) <- sites
  u <- u_vec[sites]

  P <- array(NA_real_, dim = c(p, p, p),
             dimnames = list(s1 = sites, s2 = sites, s3 = sites))
  denom <- rep(0, p); names(denom) <- sites

  for (k in seq_len(p)) {
    xk <- X[, k]
    ind_k <- !is.na(xk) & (xk > u[k])
    nk <- sum(ind_k)
    denom[k] <- nk
    if (nk == 0) next

    Xk <- X[ind_k, , drop = FALSE]
    E <- sweep(Xk, 2, u, `>`)

    for (i in seq_len(p)) {
      for (j in seq_len(p)) {
        ok_ij <- !is.na(E[, i]) & !is.na(E[, j])
        if (any(ok_ij)) {
          P[i, j, k] <- sum(E[ok_ij, i] & E[ok_ij, j]) / nk
        }
      }
    }
  }

  list(P = P, denom = denom)
}

weighted_avg_cond_probs <- function(res_list) {
  p <- dim(res_list[[1]]$P)[1]
  sites <- dimnames(res_list[[1]]$P)$s1

  num <- array(0, dim = c(p, p, p),
               dimnames = list(s1 = sites, s2 = sites, s3 = sites))
  den <- rep(0, p); names(den) <- sites

  for (r in res_list) {
    for (k in seq_len(p)) {
      d <- r$denom[k]
      if (!is.finite(d) || d <= 0) next

      tmp <- r$P[, , k]
      tmp[is.na(tmp)] <- 0

      num[, , k] <- num[, , k] + tmp * d
      den[k] <- den[k] + d
    }
  }

  out <- array(NA_real_, dim = c(p, p, p),
               dimnames = list(s1 = sites, s2 = sites, s3 = sites))

  for (k in seq_len(p)) {
    if (den[k] > 0) out[, , k] <- num[, , k] / den[k]
  }

  out
}

# observed conditional probabilities
obs_res <- lapply(seq_along(list_episodes_filtered), function(j) {
  cond_probs_episode_q_na(
    X = list_episodes_filtered[[j]],
    u_vec = u_global,
    sites = sites
  )
})

P_obs_avg <- weighted_avg_cond_probs(obs_res)

# simulated conditional probabilities, masked like observations
sim_res_by_ep <- lapply(seq_along(sims_by_ep), function(j) {
  Xobs_j <- list_episodes_filtered[[j]]

  reps <- lapply(seq_along(sims_by_ep[[j]]), function(m) {
    Xsim <- get_X_matrix(sims_by_ep[[j]][[m]])
    Xsim <- Xsim[, sites, drop = FALSE]
    Xsim_masked <- apply_obs_na_mask(Xsim, Xobs_j, sites)

    cond_probs_episode_q_na(
      X = Xsim_masked,
      u_vec = u_global,
      sites = sites
    )
  })

  P_list <- lapply(reps, `[[`, "P")
  den_list <- lapply(reps, `[[`, "denom")

  p <- dim(P_list[[1]])[1]
  num <- array(0, dim = dim(P_list[[1]]), dimnames = dimnames(P_list[[1]]))
  den <- rep(0, p); names(den) <- dimnames(P_list[[1]])$s3

  for (m in seq_along(P_list)) {
    for (k in seq_len(p)) {
      d <- den_list[[m]][k]
      if (!is.finite(d) || d <= 0) next

      tmp <- P_list[[m]][, , k]
      tmp[is.na(tmp)] <- 0

      num[, , k] <- num[, , k] + tmp * d
      den[k] <- den[k] + d
    }
  }

  Pbar <- array(NA_real_, dim = dim(num), dimnames = dimnames(num))
  for (k in seq_len(p)) {
    if (den[k] > 0) Pbar[, , k] <- num[, , k] / den[k]
  }

  list(P = Pbar, denom = Reduce(`+`, den_list) / length(den_list))
})

P_sim_avg <- weighted_avg_cond_probs(sim_res_by_ep)

# reorder sites names by alphabetical order for plotting
sites_names <- sort(sites_names)
for(s3_pick in sites_names) {
  P_obs_s3 <- P_obs_avg[, , s3_pick]
  P_sim_s3 <- P_sim_avg[, , s3_pick]

  df_plot <- rbind(
    cond_probs_to_df2d(P_obs_s3, "Observed"),
    cond_probs_to_df2d(P_sim_s3, "Simulated")
  )

  p <- ggplot(df_plot, aes(x = s1, y = s2, fill = prob)) +
    geom_tile() +
    facet_wrap(~ source, nrow = 1) +
    coord_equal() +
    scale_fill_gradient(
      low = "white",
      high = "#357470",
      limits = c(0, 1),
      na.value = "grey90"
    ) +
    labs(fill = "Probability") +
    theme_minimal() +
    xlab(expression(paste("Site s"[1]))) +
    ylab(expression(paste("Site s"[2]))) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 15),
      axis.text.y = element_text(size = 15),
      panel.grid = element_blank(),
          panel.spacing.x = unit(1, "cm"),
      legend.text = element_text(size = 14),
      legend.title = element_text(size = 16),
      strip.text = element_text(size = 16, face = "bold"),
      axis.title = element_text(size = 16)
    )

  print(p)

  foldername <- paste0(im_folder, "/swg/omsev/2025/cond_probs_all/")

  if (!dir.exists(foldername)) dir.create(foldername, recursive = TRUE)

  ggsave(
    paste0(foldername, "cond_probs_", s3_pick, ".png"),
    plot = p,
    width = 15,
    height = 7
  )
}



s3_pick <- "archie"

P_obs_s3 <- P_obs_avg[, , s3_pick]
P_sim_s3 <- P_sim_avg[, , s3_pick]



rel_diff_s3 <- (P_sim_s3 - P_obs_s3)

# éviter les divisions par zéro
rel_diff_s3[!is.finite(rel_diff_s3)] <- NA

df_plot <- cond_probs_to_df2d(rel_diff_s3, "Relative difference")



p <- ggplot(df_plot, aes(x = s1, y = s2, fill = prob)) +
  geom_tile() +
  coord_equal() +
  scale_fill_gradient2(
    low = "#570000",
    mid = "white",
    high = "#256460",
    midpoint = 0,
    limits = c(-1, 1),
    oob = scales::squish,
    na.value = "grey90"
  ) +
  labs(fill = expression(p[sim] - p[obs])) +
  theme_minimal() +
  xlab(expression(paste("Site s"[1]))) +
  ylab(expression(paste("Site s"[2]))) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 15),
    axis.text.y = element_text(size = 15),
    panel.grid = element_blank(),
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 16),
    axis.title = element_text(size = 16)
  )
p









library(ggplot2)

df_u <- data.frame(
  site = c("archie", "archiw", "brives", "cefe", "chu1", "chu2", "chu3", "chu4",
           "chu5", "chu6", "chu7", "cines", "cnrs", "crbm", "hydro", "iem",
           "mse", "poly", "um", "um35"),
  u95 = c(1.291576, 1.076260, 1.076324, 1.291576, 1.291576, 1.076260,
          1.076260, 1.076452, 1.291576, 1.076388, 1.076388, 1.076324,
          0.861136, 1.291576, 1.291576, 1.076388, 1.076535, 1.076324,
          1.291640, 1.291640)
)

df_u$site <- factor(df_u$site, levels = df_u$site[order(df_u$u95)])

p_u <- ggplot(df_u, aes(x = u95, y = site)) +
  geom_segment(aes(x = 0, xend = u95, y = site, yend = site),
               linewidth = 0.8, alpha = 0.6) +
  geom_point(size = 3) +
  labs(
    x = expression(u[0.95]~"(mm / 5 min)"),
    y = "Site"
  ) +
  theme_minimal() +
  # add the values of u95 on the right side of the points
  geom_text(aes(label = round(u95, 3), x = u95 + 0.05), hjust = 0, size = 4) +
  xlim(0, max(df_u$u95) + 0.1) +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text = element_text(size = 13),
    axis.title = element_text(size = 15)
  )

print(p_u)

ggsave(
  paste0(foldername, "threshold_u95_by_site.png"),
  plot = p_u,
  width = 8,
  height = 6
)
