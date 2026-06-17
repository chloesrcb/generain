
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

site_names <- c("crbm")
site_names <- colnames(sims_all[[1]]$X)

for (site in site_names) {

  Xsim_site <- unlist(lapply(sims_all, function(sim) as.numeric(sim$X[, site])))

  Xobs_site <- unlist(lapply(list_episodes_filtered, function(ep) as.numeric(ep[, site])))
  Xobs_site <- Xobs_site[!is.na(Xobs_site)]

  Xsim_site <- Xsim_site[Xsim_site > 1e-2]

  Xsim_meas <- apply_discretization(Xsim_site, p = 0.2152)
  Xsim_meas <- Xsim_site
  Xobs_meas <- Xobs_site

  Xsim_site <- Xsim_meas[Xsim_meas > 0]
  Xobs_site <- Xobs_meas[Xobs_meas > 0]

  config <- "above0_disc"

  df <- rbind(
    data.frame(type = "Observations", value = Xobs_site),
    data.frame(type = "Simulations", value = Xsim_site)
  )

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
  geom_density(alpha = 0.4, bw = 100) +
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
folder_site <- paste0(im_folder, "swg/omsev/2025/cumuls/")
filename <- paste0(folder_site, "cumul_density_all.png")
if (!dir.exists(folder_site)) {
  dir.create(folder_site, recursive = TRUE)
}
ggsave(filename, width = 10, height = 6)  


# idem with log(1 + cumul)

df_cumul <- df_cumul %>% filter(is.finite(cumul)) %>%
  mutate(log_cumul = log1p(cumul))

ggplot(df_cumul, aes(x = log_cumul, fill = type)) +
  geom_density(alpha = 0.4, bw = 0.5) +
  labs(x = "log(1 + Cumulative rainfall)", y = "Density", fill = "Source") +
  geom_boxplot(
    aes(y = -0.0015, group = type),
    width = 0.02,
    alpha = 0.6,
    outlier.size = 0.1
  ) +
  btf_theme + 
  coord_cartesian(xlim=c(0,10))
# save plot
folder_site <- paste0(im_folder, "swg/omsev/2025/cumuls/")
filename <- paste0(folder_site, "cumul_log_density_all.png")
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
  geom_density(alpha = 0.4, bw = 5) +
  labs(x = "Episode cumul", y = "Density", fill = "Source") +
  btf_theme +
  geom_boxplot(
    aes(y = -0.0015, group = source),
    width = 0.002,
    alpha = 0.6,
    outlier.size = 0.01
  ) +
  coord_cartesian(xlim=c(0,500))

# save plot
folder_site <- paste0(im_folder, "swg/omsev/2025/cumuls/")
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




# CUMULS OBSERVED VS SIMULATED WITH NA-MATCHING
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


folder_site <- paste0(im_folder, "swg/omsev/2025/cumuls/")
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





library(dplyr)
library(tidyr)
library(ggplot2)
library(tibble)

################################################################################
# Cumulative rainfall: observations vs simulations
################################################################################

sites_names <- colnames(rain)

Nep <- length(list_episodes_filtered)
M   <- length(sims_by_ep[[1]])

p_tip <- 0.2152

# Pseudo-discretisation des simulations continues
apply_discretization <- function(x, p = 0.2152) {
  x[!is.finite(x)] <- NA_real_
  x[x < 1e-4] <- 0
  x[x > 0 & x < p] <- p
  x[x > p & x < 2 * p] <- 2 * p
  x
}

################################################################################
# 1. Apply observation NA mask to simulations + discretize simulations
################################################################################

for (j in seq_len(Nep)) {
  obs_mat <- as.matrix(list_episodes_filtered[[j]][, sites_names, drop = FALSE])

  for (m in seq_len(M)) {
    sim_mat <- as.matrix(sims_by_ep[[j]][[m]]$X[, sites_names, drop = FALSE])

    # same missingness pattern as observations
    sim_mat[is.na(obs_mat)] <- NA

    # pseudo measurement / tipping-bucket discretisation
    sim_mat <- apply_discretization(sim_mat, p = p_tip)

    sims_by_ep[[j]][[m]]$X[, sites_names] <- sim_mat
  }
}

################################################################################
# 2. Compute episode cumulative rainfall
################################################################################

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

################################################################################
# 3. Build plotting dataframe
################################################################################

df_cumul_plot <- bind_rows(
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
  mutate(source = factor(source, levels = c("Observations", "Simulations")))

################################################################################
# 4. Density + boxplot below
################################################################################

bw_cumul <- 40

p_cumul_density <- ggplot(df_cumul_plot, aes(x = cumul, fill = source, colour = source)) +
  geom_density(
    alpha = 0.35,
    bw = bw_cumul,
    linewidth = 0.4,
    # bounds = c(0, Inf)
  ) +
  geom_boxplot(
    aes(y = -0.0008, group = source),
    width = 0.0007,
    alpha = 0.6,
    outlier.size = 0.15,
    colour = "black",
    position = position_dodge(width = 0.0009)
  ) +
  labs(
    x = "Episode cumulative rainfall",
    y = "Density",
    fill = "Source",
    colour = "Source"
  ) +
  btf_theme +
  coord_cartesian(
    xlim = c(0, 500),
    ylim = c(-0.0013, NA),
    clip = "off",
  )

p_cumul_density

# save plot
folder_site <- paste0(im_folder, "swg/omsev/2025/cumuls/")
filename <- paste0(folder_site, "cumul_density_episode_all.png")
ggsave(filename, p_cumul_density, width = 10, height = 6)
