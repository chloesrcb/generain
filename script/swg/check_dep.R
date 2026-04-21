
################################################################################
# Simulations over all episodes
################################################################################
adv_matrix <- as.matrix(adv_df_transfo[, c("vx_t", "vy_t")])
summary(adv_matrix)
speed <- sqrt(adv_matrix[, 1]^2 + adv_matrix[, 2]^2)
summary(speed)
speed[speed>150]
grid_omsev <- grid_coords_m

# convert advection from km/h to m/5min
adv_matrix <- adv_matrix * (1000 * 5 / 60)

# remove "invalid"
id_not_invalid <- which(adv_class$group != "invalid")
adv_matrix_filtered <- adv_matrix[id_not_invalid, , drop = FALSE]
list_episodes_filtered <- list_episodes[id_not_invalid]
s0_list_filtered <- s0_list[id_not_invalid]
u_list_filtered <- u_list[id_not_invalid]

M <- 100 # repetitions per episode
Nsim <- length(list_episodes_filtered)  # number of episodes to simulate
sims_by_ep <- vector("list", Nsim)
u_sim <- numeric(Nsim)
timesteps <- 0:11 #* 5 / 60 # every 5 minutes for 1 hour
set.seed(123)
for (j in seq_len(Nsim)) {
  ep_idx <- j
  s0_j  <- s0_list_filtered[ep_idx]
  u_j   <- u_list_filtered[ep_idx]
  adv_j <- adv_matrix_filtered[ep_idx, ]
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

N <- Nsim
length(sims_by_ep) * length(sims_by_ep[[1]])

# sims_by_ep[[i]][[m]] -> sims_all[[k]]
sims_all <- unlist(sims_by_ep, recursive = FALSE)

site_names <- colnames(sims_all[[1]]$X)


library(ggplot2)
library(dplyr)



# OBS cumul par épisode
Cobs <- sapply(seq_len(N), function(j) sum(as.matrix(list_episodes_filtered[[j]]), na.rm=TRUE))

# SIM cumul par épisode (N x M)
Csim <- matrix(NA_real_, nrow=N, ncol=M)
for (j in seq_len(N)) {
  for (m in seq_len(M)) {
    Csim[j,m] <- sum(sims_by_ep[[j]][[m]]$X, na.rm=TRUE)
  }
}

summary(Cobs)
summary(as.vector(Csim))

Cobs_pos <- Cobs
Csim_pos <- as.vector(Csim)
Csim_pos <- Csim_pos

# Violin sur log1p
library(ggplot2)
df <- rbind(
  data.frame(type="Observations", value=Cobs_pos),
  data.frame(type="Simulations", value=Csim_pos)
)
df$y <- log1p(df$value)

ggplot(df, aes(x=type, y=y)) +
  geom_violin(trim=FALSE, alpha=0.7, color=NA, scale = "width", bounds = c(0, Inf), bw = 0.5, fill = btfgreen) +
  geom_boxplot(width=0.12, outlier.shape=NA, fill="white") +
  theme_minimal() +
  labs(y="log(1 + Cumul)", x="") +
  btf_theme


q05 <- apply(Csim, 1, quantile, probs=0.05, na.rm=TRUE)
q95 <- apply(Csim, 1, quantile, probs=0.95, na.rm=TRUE)
covered90 <- (Cobs >= q05) & (Cobs <= q95)

mean(covered90, na.rm=TRUE)




q <- quantile(Cobs, probs = c(0.33, 0.66), na.rm=TRUE)

class_obs <- cut(
  Cobs,
  breaks = c(-Inf, q[1], q[2], Inf),
  labels = c("weak", "medium", "strong")
)


# function to discretize low values in sim to match obs
# low obs = 0, p, 2*p, etc
apply_discretization <- function(x, p) {
  x[x < 1e-2] <- 0
  x[x > 0 & x < p] <- p
  x[x > p & x < 2*p] <- 2*p
  # x[x > 2*p & x < 3*p] <- 3*p
  x
}


apply_measurement <- function(x, p) {
  x[x < 1e-2] <- 0
  x[x > 0 & x < p] <- p
  x
}


sims_all <- unlist(sims_by_ep, recursive = FALSE)

site_names <- c("cefe", "iem", "um", "cnrs")
site_names <- colnames(sims_all[[1]]$X)

for (site in site_names) {

  Xsim_site <- unlist(lapply(sims_all, function(sim) as.numeric(sim$X[, site])))

  Xobs_site <- unlist(lapply(list_episodes_filtered, function(ep) as.numeric(ep[, site])))
  Xobs_site <- Xobs_site[!is.na(Xobs_site)]

  Xsim_site <- Xsim_site[Xsim_site > 1e-2]

  # Xsim_meas <- apply_discretization(Xsim_site, p = 0.2152)
  Xsim_meas <- Xsim_site
  Xobs_meas <- Xobs_site

  Xsim_site <- Xsim_meas[Xsim_meas > 0]
  Xobs_site <- Xobs_meas[Xobs_meas > 0]

  config <- "above0"

  df <- rbind(
    data.frame(type = "Observations", value = Xobs_site),
    data.frame(type = "Simulations", value = Xsim_site)
  )

  bw_common <- bw.nrd0(df$value)

  df_log <- df %>% mutate(yplot = log1p(value))

  psite_log <- ggplot(df_log, aes(x = type, y = yplot)) +
    geom_violin(alpha = 0.7, trim = FALSE, color = NA,
                bw = 0.2, scale = "width", bounds= c(0, Inf), fill = btfgreen) +
    labs(x = "", y = "log(1 + Rain)") +
    theme_minimal() + btf_theme +
    geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA) +
    theme(legend.position = "none")

  folder_site <- paste0(im_folder, "swg/omsev/margins/", config, "/")
  if (!dir.exists(folder_site)) dir.create(folder_site, recursive = TRUE)

  ggsave(paste0(folder_site, "violin_log_", config, "_", site, ".png"),
         psite_log, width = 6, height = 4)

  psite_dens <- ggplot(df, aes(x = value, fill = type)) +
            geom_density(alpha = 0.4, bw = 0.2) +

            # boxplot en bas
            geom_boxplot(
              aes(y = -0.1, group = type),
              width = 0.15,
              alpha = 0.6,
              outlier.size = 0.8
            ) +

            labs(x = "Rainfall", y = "Density", fill = "Source") +
            theme_minimal() + btf_theme +
            xlim(0, 10)

  ggsave(paste0(folder_site, "density_box_", config, "_", site, ".png"),
         psite_dens, width = 8, height = 6)
}


# CUMULATIVE RAINFALL PER EPISODE AND SITE

sites_names <- colnames(rain)

# gestion of NA values
# in data if NA in obs, set sim to NA at the same episode and site
for (j in seq_len(Nsim)) {
  for (site in sites_names) {
    if (all(is.na(list_episodes_filtered[[j]][,site]))) {
      for (m in seq_along(sims_by_ep[[j]])) {
        sims_by_ep[[j]][[m]]$X[, site] <- NA
      }
    }
  }
}

# apply apply_measurement to all sims
for (j in seq_len(Nsim)) {
  for (m in seq_along(sims_by_ep[[j]])) {
    sims_by_ep[[j]][[m]]$X <- apply_discretization(sims_by_ep[[j]][[m]]$X, p = 0.2152)
  }
}

cumul_ep <- c()
for (j in seq_len(Nsim)) {
  ep_mat <- list_episodes_filtered[[j]]
  cumul_ep[j] <- sum(ep_mat, na.rm = TRUE)
}

df_cumul_obs <- data.frame(
  ep_id = seq_len(Nsim),
  cumul = cumul_ep
)

cumul_sim <- matrix(NA_real_, nrow=N, ncol=M)
mean_ep_sim <- numeric(Nsim)
for (j in seq_len(Nsim)) {
  cumul_j <- numeric(M)
  for (m in seq_len(M)) {
    cumul_sim[j,m] <- sum(sims_by_ep[[j]][[m]]$X, na.rm = TRUE)
    if(cumul_sim[j,m] < 1) {
      print(paste("Episode", j, "Rep", m))
      print(paste("Cumul sim:", cumul_sim[j,m]))   
     }
    cumul_j[m] <- cumul_sim[j,m]
  }
  mean_ep_sim[j] <- mean(cumul_j, na.rm = TRUE)
}

# [1] "Episode 237 Rep 27"
# [1] "Cumul sim: 0.865534598725329"
# ep_id <- 237
# rep <- 27
# sims_by_ep[[ep_id]][[rep]]$X

min_cum_sim <- min(mean_ep_sim, na.rm = TRUE)

df_cumul_sim <- as.data.frame(cumul_sim) %>%
  pivot_longer(cols = everything(), names_to = "rep", values_to = "cumul") %>%
  mutate(rep = as.integer(gsub("V", "", rep)))

df_cumul <- rbind(
  data.frame(type = "Observations", cumul = cumul_ep),
  data.frame(type = "Simulations", cumul = as.vector(cumul_sim))
)

# # # plot density
# ggplot(df_cumul_sim, aes(x = cumul)) +
#   geom_density(alpha = 0.4, bw = 20) +
#   labs(x = "Cumulative rainfall", y = "Density") +
#   btf_theme + 
#   xlim(0, 300)

# ggplot(df_cumul_obs, aes(x = cumul)) +
#   geom_density(alpha = 0.4, bw = 20) +
#   labs(x = "Cumulative rainfall", y = "Density") +
#   btf_theme + 
#   xlim(0, 300)

min_cum_sim <- min(df_cumul_sim$cumul, na.rm = TRUE)

df_cumul <- df_cumul %>% filter(is.finite(cumul))

bw_common <- bw.nrd(df_cumul$cumul)
ggplot(df_cumul, aes(x = cumul, fill = type)) +
  geom_density(alpha = 0.4, bw = 40) +
  labs(x = "Cumulative rainfall", y = "Density", fill = "Source") +
  geom_boxplot(
    aes(y = -0.0015, group = type),
    width = 0.002,
    alpha = 0.6,
    outlier.size = 0.1
  ) +
  btf_theme + 
  coord_cartesian(xlim=c(0,700))

# save plot
folder_site <- paste0(im_folder, "swg/omsev/cumuls/")
filename <- paste0(folder_site, "cumul_density_all.png")
if (!dir.exists(folder_site)) {
  dir.create(folder_site, recursive = TRUE)
}
ggsave(filename, width = 10, height = 6)  


# Cumulatives observed and simulated
cum_obs_mat <- t(vapply(
  seq_along(list_episodes_filtered),
  function(i) sum_over_time_by_site(list_episodes_filtered[[i]], sites_names),
  FUN.VALUE = setNames(numeric(length(sites_names)), sites_names)
))

cum_obs_long <- as.data.frame(cum_obs_mat) |>
  mutate(ep_id = seq_len(nrow(cum_obs_mat))) |>
  pivot_longer(-ep_id, names_to = "site", values_to = "cum_obs")

cum_sim_long <- purrr::map_dfr(seq_along(sims_by_ep), function(j) {
  purrr::map_dfr(seq_along(sims_by_ep[[j]]), function(m) {
    simX <- sims_by_ep[[j]][[m]]$X
    simX <- apply_measurement(simX, p = 0.2152)
    v <- sum_over_time_by_site(simX, sites_names)
    tibble(
      ep_id = j,
      rep   = m,
      site  = names(v),
      cum_sim = as.numeric(v)
    )
  })
})

selected_points_filtered <- selected_points[selected_points$speed_class != "invalid", ]

meta_ep <- tibble(
  ep_id = seq_len(length(list_episodes_filtered)),
  speed_class =   selected_points_filtered$speed_class
)

# Join OBS + SIM ep site
df_cum <- cum_obs_long |>
  left_join(meta_ep, by = "ep_id") |>
  left_join(cum_sim_long, by = c("ep_id", "site"))
table(df_cum$speed_class)

# Cumuls observed and simulated per site
df_cum_site <- df_cum |>
  group_by(ep_id, site) |>
  summarise(
    cum_obs = first(cum_obs),
    cum_sim = mean(cum_sim, na.rm = FALSE),
    .groups = "drop"
  )


df_cum <- cum_obs_long |>
  left_join(meta_ep, by = "ep_id") |>
  left_join(cum_sim_long, by = c("ep_id", "site"))
unique(df_cum$speed_class)
max(df_cum$cum_obs, na.rm=TRUE)
max(df_cum$cum_sim, na.rm=TRUE)
df_cum[is.na(df_cum$cum_obs), ]$cum_sim <- NA

head(sort(unique(df_cum$cum_obs)))
head(sort(unique(df_cum$cum_sim)))
summary(df_cum$cum_obs)
summary(df_cum$cum_sim)
obs_ep <- cum_obs_long |>
  group_by(ep_id) |>
  summarise(cumul_obs = sum(cum_obs, na.rm=TRUE), .groups="drop")

sim_eprep <- cum_sim_long |>
  group_by(ep_id, rep) |>
  summarise(cumul_sim = sum(cum_sim, na.rm=TRUE), .groups="drop")

df_all <- bind_rows(
  obs_ep |> transmute(ep_id, source="obs", cumul=cumul_obs),
  sim_eprep |> transmute(ep_id, source="sim", cumul=cumul_sim)
)

df_plot <- df_all |>
  left_join(meta_ep, by = "ep_id")
# density
ggplot(df_plot, aes(x = cumul, fill = source)) +
  geom_density(alpha = 0.4, bw = 40) +
  labs(x = "Episode cumul", y = "Density", fill = "Source") +
  btf_theme +
  geom_boxplot(
    aes(y = -0.0015, group = source),
    width = 0.002,
    alpha = 0.6,
    outlier.size = 0.1
  ) +
  coord_cartesian(xlim=c(0,700))

# save plot
folder_site <- paste0(im_folder, "swg/omsev/cumuls/")
filename <- paste0(folder_site, "cumul_density_episode_all_disc.png")
if (!dir.exists(folder_site)) {
  dir.create(folder_site, recursive = TRUE)
}
ggsave(filename, width = 10, height = 6)

df_cumul <- df_cumul %>% filter(is.finite(cumul))



cumul_episode <- function(
  ep_mat,
  sites,
  s0,
  radius,
  time_idx = NULL,
  transfo = FALSE
) {

  x <- as.matrix(ep_mat)
  if (transfo) {
     x <- apply_measurement(x, p = 0.2152)
  }

  is_time_by_site <- !is.null(colnames(x)) && all(colnames(x) %in% sites)
  if (is_time_by_site) {
    x <- x[, sites, drop = FALSE]
    nT <- nrow(x)
  } else {
    x <- x[sites, , drop = FALSE]
    nT <- ncol(x)
  }

  if (!is.null(time_idx)) {
    time_idx2 <- time_idx[time_idx >= 1 & time_idx <= nT]
    if (length(time_idx2) == 0) return(NA_real_)

    if (is_time_by_site) {
      x <- x[time_idx2, , drop = FALSE]
    } else {
      x <- x[, time_idx2, drop = FALSE]
    }
  }

  sum(x, na.rm = TRUE)
}




# ============================================================================
# CUMULS OBSERVED VS SIMULATED WITH NA-MATCHING
# ============================================================================

library(dplyr)
library(tidyr)
library(ggplot2)
library(purrr)
library(tibble)

sites_names <- colnames(rain)

Nep <- length(list_episodes_filtered)
M   <- length(sims_by_ep[[1]])

for (j in seq_along(list_episodes_filtered)) {
  obs_mat <- as.matrix(list_episodes_filtered[[j]][, sites_names, drop = FALSE])

  for (m in seq_along(sims_by_ep[[j]])) {
    sim_mat <- as.matrix(sims_by_ep[[j]][[m]]$X[, sites_names, drop = FALSE])

    sim_mat[is.na(obs_mat)] <- NA

    sims_by_ep[[j]][[m]]$X[, sites_names] <- sim_mat
  }
}

for (j in seq_along(sims_by_ep)) {
  for (m in seq_along(sims_by_ep[[j]])) {
    sims_by_ep[[j]][[m]]$X <- apply_discretization(
      sims_by_ep[[j]][[m]]$X,
      p = 0.2152
    )
  }
}

cumul_obs_raw  <- numeric(Nep)
cumul_obs_norm <- numeric(Nep)
n_avail_obs    <- integer(Nep)

cumul_sim_raw  <- matrix(NA_real_, nrow = Nep, ncol = M)
cumul_sim_norm <- matrix(NA_real_, nrow = Nep, ncol = M)

for (j in seq_along(list_episodes_filtered)) {
  obs_mat <- as.matrix(list_episodes_filtered[[j]][, sites_names, drop = FALSE])

  n_avail_obs[j] <- sum(!is.na(obs_mat))
  cumul_obs_raw[j] <- sum(obs_mat, na.rm = TRUE)

  if (n_avail_obs[j] > 0) {
    cumul_obs_norm[j] <- cumul_obs_raw[j] / n_avail_obs[j]
  } else {
    cumul_obs_norm[j] <- NA_real_
  }

  for (m in seq_along(sims_by_ep[[j]])) {
    sim_mat <- as.matrix(sims_by_ep[[j]][[m]]$X[, sites_names, drop = FALSE])

    cumul_sim_raw[j, m] <- sum(sim_mat, na.rm = TRUE)

    if (n_avail_obs[j] > 0) {
      cumul_sim_norm[j, m] <- cumul_sim_raw[j, m] / n_avail_obs[j]
    } else {
      cumul_sim_norm[j, m] <- NA_real_
    }
  }
}

df_cumul_obs_raw <- tibble(
  ep_id = seq_len(Nep),
  cumul = cumul_obs_raw,
  source = "Observations"
)

df_cumul_obs_norm <- tibble(
  ep_id = seq_len(Nep),
  cumul = cumul_obs_norm,
  source = "Observations"
)

df_cumul_sim_raw <- as.data.frame(cumul_sim_raw) %>%
  mutate(ep_id = seq_len(n())) %>%
  pivot_longer(
    cols = -ep_id,
    names_to = "rep",
    values_to = "cumul"
  ) %>%
  mutate(
    rep = as.integer(gsub("V", "", rep)),
    source = "Simulations"
  )

df_cumul_sim_norm <- as.data.frame(cumul_sim_norm) %>%
  mutate(ep_id = seq_len(n())) %>%
  pivot_longer(
    cols = -ep_id,
    names_to = "rep",
    values_to = "cumul"
  ) %>%
  mutate(
    rep = as.integer(gsub("V", "", rep)),
    source = "Simulations"
  )

df_plot_raw <- bind_rows(
  df_cumul_obs_raw,
  df_cumul_sim_raw %>% select(ep_id, cumul, source)
) %>%
  filter(is.finite(cumul))

df_plot_norm <- bind_rows(
  df_cumul_obs_norm,
  df_cumul_sim_norm %>% select(ep_id, cumul, source)
) %>%
  filter(is.finite(cumul))

cat("Raw cumulative rainfall summaries:\n")
print(df_plot_raw %>% group_by(source) %>% summarise(
  n = n(),
  min = min(cumul, na.rm = TRUE),
  q25 = quantile(cumul, 0.25, na.rm = TRUE),
  median = median(cumul, na.rm = TRUE),
  mean = mean(cumul, na.rm = TRUE),
  q75 = quantile(cumul, 0.75, na.rm = TRUE),
  max = max(cumul, na.rm = TRUE)
))

cat("\nNormalized cumulative rainfall summaries:\n")
print(df_plot_norm %>% group_by(source) %>% summarise(
  n = n(),
  min = min(cumul, na.rm = TRUE),
  q25 = quantile(cumul, 0.25, na.rm = TRUE),
  median = median(cumul, na.rm = TRUE),
  mean = mean(cumul, na.rm = TRUE),
  q75 = quantile(cumul, 0.75, na.rm = TRUE),
  max = max(cumul, na.rm = TRUE)
))

p_raw <- ggplot(df_plot_raw, aes(x = cumul, fill = source)) +
  geom_density(alpha = 0.4, bw = 40) +
  labs(
    x = "Episode cumulative rainfall",
    y = "Density",
    fill = "Source"
  ) +
  geom_boxplot(
    aes(y = -0.0015, group = source),
    width = 0.002,
    alpha = 0.6,
    outlier.size = 0.1
  ) +
  btf_theme +
  coord_cartesian(xlim = c(0, 700))

print(p_raw)

p_norm <- ggplot(df_plot_norm, aes(x = cumul, fill = source)) +
  geom_density(alpha = 0.4, bw = 0.8) +
  labs(
    x = "Normalized episode cumulative rainfall",
    y = "Density",
    fill = "Source"
  ) +
  geom_boxplot(
    aes(y = -0.05, group = source),
    width = 0.1,
    alpha = 0.5,
    outlier.size = 0.1
  ) +
  btf_theme +
  coord_cartesian(xlim = c(0, 10))

print(p_norm)


folder_site <- paste0(im_folder, "swg/omsev/cumuls/")
if (!dir.exists(folder_site)) {
  dir.create(folder_site, recursive = TRUE)
}

filename_raw <- paste0(folder_site, "cumul_density_episode_all_disc_masked.png")
ggsave(filename_raw, p_raw, width = 10, height = 6)

filename_norm <- paste0(folder_site, "cumul_density_episode_all_disc_masked_normalized.png")
ggsave(filename_norm, p_norm, width = 10, height = 6)

cat("Saved:\n")
cat(filename_raw, "\n")
cat(filename_norm, "\n")


mean_ep_sim_raw <- apply(cumul_sim_raw, 1, mean, na.rm = TRUE)
mean_ep_sim_norm <- apply(cumul_sim_norm, 1, mean, na.rm = TRUE)

df_compare_means <- tibble(
  ep_id = seq_len(Nep),
  obs_raw = cumul_obs_raw,
  sim_mean_raw = mean_ep_sim_raw,
  obs_norm = cumul_obs_norm,
  sim_mean_norm = mean_ep_sim_norm,
  n_avail_obs = n_avail_obs
)

head(df_compare_means)
summary(df_compare_means)

# ----------------------------------------------------------------------------
# 10) Optional scatterplots obs vs mean simulated
# ----------------------------------------------------------------------------
p_scatter_raw <- ggplot(df_compare_means, aes(x = obs_raw, y = sim_mean_raw)) +
  geom_point(alpha = 0.6) +
  geom_abline(intercept = 0, slope = 1, linetype = 2, color = "red") +
  labs(
    x = "Observed episode cumulative rainfall",
    y = "Mean simulated episode cumulative rainfall"
  ) +
  btf_theme +
  coord_cartesian(xlim = c(0, 700), ylim = c(0, 700))

print(p_scatter_raw)

filename_scatter_raw <- paste0(folder_site, "cumul_scatter_obs_vs_sim_mean_raw.png")
ggsave(filename_scatter_raw, p_scatter_raw, width = 7, height = 6)

p_scatter_norm <- ggplot(df_compare_means, aes(x = obs_norm, y = sim_mean_norm)) +
  geom_point(alpha = 0.6) +
  geom_abline(intercept = 0, slope = 1, linetype = 2, color = "red") +
  labs(
    x = "Observed normalized cumulative rainfall",
    y = "Mean simulated normalized cumulative rainfall"
  ) +
  btf_theme

print(p_scatter_norm)

filename_scatter_norm <- paste0(folder_site, "cumul_scatter_obs_vs_sim_mean_normalized.png")
ggsave(filename_scatter_norm, p_scatter_norm, width = 7, height = 6)

cat(filename_scatter_raw, "\n")
cat(filename_scatter_norm, "\n")


df_plot_raw %>%
  group_by(source) %>%
  summarise(
    q50 = quantile(cumul, 0.50, na.rm = TRUE),
    q75 = quantile(cumul, 0.75, na.rm = TRUE),
    q90 = quantile(cumul, 0.90, na.rm = TRUE),
    q95 = quantile(cumul, 0.95, na.rm = TRUE),
    q99 = quantile(cumul, 0.99, na.rm = TRUE),
    max = max(cumul, na.rm = TRUE)
  )


summary(df_compare_means$sim_mean_raw - df_compare_means$obs_raw)
cor(df_compare_means$obs_raw, df_compare_means$sim_mean_raw, use = "complete.obs")

df_compare_means %>%
  mutate(ratio = sim_mean_raw / obs_raw) %>%
  ggplot(aes(x = ratio)) +
  geom_histogram(bins = 50) +
  btf_theme

ggplot(df_compare_means, aes(x = obs_raw, y = sim_mean_raw)) +
  geom_point(alpha = 0.6) +
  geom_abline(slope = 1, intercept = 0, linetype = 2, color = "red") +
  btf_theme
