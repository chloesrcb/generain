
# remove objects from the environment
rm(list = ls())

# clear console
cat("\014")

# Load libraries and set theme
source("./script/load_libraries.R")

# library(generain)
library(ggplot2)
library(reshape2)
library(animation)
library(RandomFields)
library(RandomFieldsUtils)
library(dplyr)
library(tidyr)
library(lubridate)
library(purrr)
library(viridis)
library(magick)


# Get all files in the folder "R"
functions_folder <- "./R"
files <- list.files(functions_folder, full.names = TRUE)
# load all functions in files
invisible(lapply(files, function(f) source(f, echo = FALSE)))
library(latex2exp)

################################################################################
# FOCUS ON A SINGLE EPISODE
################################################################################


# get rain data
filename_rain <- paste0(data_folder, "omsev/omsev_5min/rain_mtp_5min_2019_2024_cleaned.csv")
rain_omsev <- read.csv(filename_rain)
head(rain_omsev)

# egpd fit
filename_egpd <- paste0(data_folder, "../thesis/resources/images/EGPD/OMSEV/2019_2024/egpd_results.csv")
egpd_params <- read.csv(filename_egpd)

# put dates as rownames
rownames(rain_omsev) <- rain_omsev$dates
rain <- rain_omsev[-1] # remove dates column

# get location of each rain gauge
filename_loc <- paste0(data_folder, "omsev/loc_rain_gauges.csv")
location_gauges <- read.csv(filename_loc)
location_gauges$Station <- c("iem", "mse", "poly", "um", "cefe", "cnrs",
                             "crbm", "archiw", "archie", "um35", "chu1",
                             "chu2", "chu3", "chu4", "chu5", "chu6", "chu7",
                             "cines", "brives", "hydro")


rain <- rain[, !(colnames(rain) %in% c("cines", "hydro", "brives"))]
location_gauges <- location_gauges[location_gauges$Station != "cines" &
                                   location_gauges$Station != "hydro" &
                                   location_gauges$Station != "brives", ]
dist_mat <- get_dist_mat(location_gauges)
df_dist <- reshape_distances(dist_mat)

sites_names <- colnames(rain)

sites_coords <- location_gauges[, c("Longitude", "Latitude")]

rownames(sites_coords) <- location_gauges$Station
sites_coords_sf <- st_as_sf(sites_coords, coords = c("Longitude", "Latitude"),
                            crs = 4326)
sites_coords_sf <- st_transform(sites_coords_sf, crs = 2154)
coords_m <- st_coordinates(sites_coords_sf)
grid_omsev <- as.data.frame(coords_m / 1000)
colnames(grid_omsev) <- c("Longitude", "Latitude")
rownames(grid_omsev) <- rownames(sites_coords)


etas_estimates <- c(0.6732099, 5.9844737)

params_vario <- list(
  beta1 = 1.5129632,
  beta2 = 4.4833122,
  alpha1 = 0.2138693,
  alpha2 = 0.7012997
)

convert_params <- function(beta1, beta2, alpha1, alpha2, c_x = 1, c_t = 1) {
  beta1_new <- beta1 / (c_x^alpha1)
  beta2_new <- beta2 / (c_t^alpha2)
  list(beta1 = beta1_new, beta2 = beta2_new)
}

# convert params from km/h to m/5min
c_x_m <- 1000    # for m
c_t_5min <- 12   # 1 hour = 12 * 5min

# convert params and ci to m/5min
params_m5min <- convert_params(result$par[1], result$par[2],
                               result$par[3], result$par[4],
                               c_x = c_x_m, c_t = c_t_5min)
params_vario <- list(
  beta1 = params_m5min$beta1,
  beta2 = params_m5min$beta2,
  alpha1 = result$par[3],
  alpha2 = result$par[4]
)

# etas_estimates <- c(0.868, 1.688)

# [1] 1.6854191 4.4309430 0.4847687 0.6759399 1 ,1



# u_s <- apply(rain, 2, function(x) quantile(x[x > 0], probs = q, na.rm = TRUE))

q <- 0.95
delta <- 12
dmin <- 1200  # m

# in rain remove when all data are NA
set_st_excess <- get_spatiotemp_excess(rain, quantile = q, remove_zeros = TRUE)


# Spatio-temporal neighborhood parameters
s0t0_set <- get_s0t0_pairs(grid_omsev, rain,
                            min_spatial_dist = dmin,
                            episode_size = delta,
                            set_st_excess = set_st_excess,
                            n_max_episodes = 10000,
                            latlon = FALSE)

selected_points <- s0t0_set
selected_points <- selected_points %>%
  mutate(t0_date = as.POSIXct(t0_date, format = "%Y-%m-%d %H:%M:%S", tz = "UTC"))

n_episodes <- length(selected_points$s0)
t0_list <- selected_points$t0
s0_list <- selected_points$s0
u_list <- selected_points$u_s0
list_episodes_points <- get_extreme_episodes(selected_points, rain,
                                     episode_size = delta, unif = FALSE)
list_episodes <- list_episodes_points$episodes

adv_filename <- paste(data_folder, "/omsev/adv_estim/combined_comephore_omsev/episode_advection_q",
                          q * 100, "_delta", delta, "_dmin", dmin,
                          ".csv", sep = "")
adv_df <- read.csv(adv_filename, sep = ",")

adv_df_transfo <- adv_df
adv_df_transfo$vnorm <- sqrt(adv_df$vx^2 + adv_df$vy^2)
adv_df_transfo$vnorm_t <- etas_estimates[1] * adv_df_transfo$vnorm^etas_estimates[2]

adv_df_transfo$vx <- ifelse(
  adv_df_transfo$vnorm > 0,
  adv_df_transfo$vx / adv_df_transfo$vnorm * adv_df_transfo$vnorm_t,
  0
)

adv_df_transfo$vy <- ifelse(
  adv_df_transfo$vnorm > 0,
  adv_df_transfo$vy / adv_df_transfo$vnorm * adv_df_transfo$vnorm_t,
  0
)

plot(adv_df$vx, adv_df$vy, pch = 19,
     xlab = "vx", ylab = "vy")

plot(adv_df_transfo$vx, adv_df_transfo$vy, pch = 19,
     xlab = "vx", ylab = "vy", xlim = c(-10,10), ylim = c(-10,10))

calm_threshold <- 0.1
# Compute speed
speed <- sqrt(adv_df$vx^2 + adv_df$vy^2)
is_calm <- speed < calm_threshold
summary(speed)
direction <- atan2(adv_df$vy, adv_df$vx) * (180 / pi)
# make sure direction is in [0, 360]
angle_deg <- (direction + 360) %% 360

speed_class <- ifelse(
  is_calm,
  "calm",
  as.character(cut(
    speed,
    breaks = c(0, 1, 2, Inf),
    labels = c("weak", "moderate", "strong"),
    include.lowest = TRUE
  ))
)

speed_class <- factor(speed_class, levels = c("calm","weak","moderate","strong"))

direction_class <- ifelse(
  is_calm,
  "none",
  as.character(cut(
    angle_deg,
    breaks = c(0, 90, 180, 270, 360),
    labels = c("E", "S", "W", "N"),
    include.lowest = TRUE
  ))
)
direction_class <- factor(direction_class, levels = c("none","E","S","W","N"))


adv_class <- data.frame(speed_class)
table(adv_class)
adv_class$group <- ifelse(
  speed_class == "calm",
  "calm",
  paste(speed_class)
)
table(adv_class$group)

adv_df <- adv_df_transfo

# convert adv to m/5min from km/h
adv_df$vx_final <- adv_df$vx * (1000 / 60) * 5
adv_df$vy_final <- adv_df$vy * (1000 / 60) * 5

episode_idx <- 279
real_episode <- list_episodes[[episode_idx]]
# remove columns with all NA
real_episode <- real_episode[, colSums(is.na(real_episode)) < nrow(real_episode)]
s0_real <- s0_list[episode_idx]
t0_real <- t0_list[episode_idx]
u_s0_real <- u_list[episode_idx]
adv_real <- as.numeric(adv_df[episode_idx, c("vx_final", "vy_final")])

coord_s0_real <- as.numeric(grid_omsev[rownames(grid_omsev) == s0_real, ])

ep_matrix <- list_episodes[[episode_idx]]
ep_values <- unlist(ep_matrix, use.names = FALSE)
ep_values <- ep_values[!is.na(ep_values)]

episode <- list_episodes[[episode_idx]]

res <- compute_p0_all_episodes(list_episodes)

p0_episode_obs <- res$p0_by_episode
p0_mean_by_site <- res$p0_mean
p0_target <- p0_mean_by_site

p0_values <- sapply(colnames(rain), function(s) {
  mean(rain[, s] == 0, na.rm = TRUE)
})
kappa_vect <- egpd_params[match(colnames(rain), egpd_params$Site), "kappa"]
xi_vect    <- egpd_params[match(colnames(rain), egpd_params$Site), "xi"]
sigma_vect <- egpd_params[match(colnames(rain), egpd_params$Site), "sigma"]
params_margins <- list(
  xi    = xi_vect,
  sigma = sigma_vect,
  kappa = kappa_vect,
  p0    = p0_values
)
times <- 0:(delta - 1)


selected_points$adv_group <- adv_class$group
head(selected_points)


X_s0_t0 <- numeric(length(list_episodes))
for (i in seq_along(list_episodes)) {
  ep <- list_episodes[[i]]
  s0 <- s0_list[i]
  X_s0_t0[i] <- ep[1, s0]
}

selected_points$X_s0_t0 <- X_s0_t0


q_init <- quantile(selected_points$X_s0_t0, probs = c(0.25, 0.75))

selected_points$init_class <- cut(
  selected_points$X_s0_t0,
  breaks = c(-Inf, q_init[1], q_init[2], Inf),
  labels = c("init_low", "init_mid", "init_high"),
  include.lowest = TRUE
)

head(selected_points)

# concat adv_class and init_class
# adv_init_class <- data.frame(
#   adv_group = adv_class$group,
#   init_class = selected_points$init_class
# )
# adv_init_class$group <- paste(adv_init_class$adv_group,
#                               adv_init_class$init_class,
#                               sep = "_")

# table(adv_init_class$group)

# selected_points$adv_init_group <- adv_init_class$group
# head(selected_points)

# Clustering of episodes with similar advection in list_episodes

adv_matrix <- as.matrix(adv_df[, c("vx_final", "vy_final")])
all_group_names <- unique(selected_points$adv_group)
group_adv <- "moderate"  # choose one group to simulate
#get list_episodes and s0_list for this group
indices_group <- which(selected_points$adv_group == group_adv)
list_episodes_group <- list_episodes[indices_group]
length(list_episodes_group)
s0_list_group <- s0_list[indices_group]
u_list_group <- u_list[indices_group]
adv_group <- adv_matrix[indices_group, ]

Nsim <- 1000
s0_sim <- integer(Nsim)
sims_group <- vector("list", Nsim)
for (i in seq_len(Nsim)) {
  idx <- sample(seq_along(list_episodes_group), 1)
  episode_i <- list_episodes_group[[idx]]
  s0_i <- s0_list_group[idx]
  s0_sim[i] <- s0_i
  # X_s0_t0 <- selected_points$X_s0_t0[indices_group[idx]]
  # cat("Simulating ", i, "/", Nsim,
  #     " (episode idx:", indices_group[idx],
  #     " s0:", s0_i,
  #     " X_s0_t0:", round(X_s0_t0,2), ")\n")
  
  sims_group[[i]] <- simulate_many_episodes(
    N = 1,
    u = 1000,
    u_emp = u_list_group[idx],
    params_vario = params_vario,
    params_margins = params_margins,
    coords = grid_omsev,
    times = times,
    adv = adv_group,
    t0 = 0,
    s0 = s0_i
  )[[1]]
}

df_sims <- do.call(rbind, lapply(seq_len(Nsim), function(i) {
  s0_i <- s0_sim[i]
  data.frame(
    time = seq_len(ncol(sims_group[[i]])),
    rainfall = sims_group[[i]][s0_i, ],
    sim_id = i
  )
}))


df_sim_summary <- df_sims %>%
  group_by(time) %>%
  summarise(
    sim_mean = mean(rainfall),
    sim_median = median(rainfall),
    sim_sd  = sd(rainfall),
    sim_low  = quantile(rainfall, 0.025),
    sim_high = quantile(rainfall, 0.975)
  )


df_real_group <- do.call(rbind, lapply(seq_along(list_episodes_group), function(idx) {
  ep <- list_episodes_group[[idx]]
  s0_i <- s0_list_group[idx]
  data.frame(
    time = seq_len(nrow(ep)),
    rainfall = ep[, s0_i],
    episode_id = idx
  )
}))


df_real_summary <- df_real_group %>%
  group_by(time) %>%
  summarise(
    real_mean = mean(rainfall),
    real_median = median(rainfall),
    real_low  = quantile(rainfall, 0.025),
    real_high = quantile(rainfall, 0.975)
  )

ggplot() +
  geom_ribbon(
    data = df_sim_summary,
    aes(x = time, ymin = sim_low, ymax = sim_high, fill = "Simulations"),
    alpha = 0.25
  ) +
  geom_line(
    data = df_sim_summary,
    aes(x = time, y = sim_median, color = "Simulations"),
    linewidth = 1
  ) +
  geom_line(
    data = df_real_summary,
    aes(x = time, y = real_median, color = "Observed episodes"),
    linewidth = 1
  ) +  geom_ribbon(
    data = df_real_summary,
    aes(x = time, ymin = real_low, ymax = real_high, fill = "Observed episodes"),
    alpha = 0.25
  ) +
  ylim(0, 11)+
  scale_color_manual(values = c("Simulations"="#6b6b6b","Observed episodes"="#c73535")) +
  scale_fill_manual(values = c("Simulations"="#bdbdbd","Observed episodes"="#e28b8b")) +
  labs(x="Time step", y="Rainfall (mm/5min)", color="", fill="") +
  theme_minimal() +
  btf_theme





# save this plot
ggsave(paste0(im_folder, "swg/omsev/sim_vs_real/CI_95sim_95obs_adv_group_", group_adv, ".png"),
       width = 8, height = 6, dpi = 300)

ggplot() +
  geom_ribbon(
    data = df_sim_summary,
    aes(x = time, ymin = sim_low, ymax = sim_high, fill = "Simulations"),
    alpha = 0.25
  ) +
  geom_line(
    data = df_sim_summary,
    aes(x = time, y = sim_median, color = "Median simulations"),
    size = 1
  ) +
  geom_line(
    data = df_real_group,
    aes(x = time, y = rainfall, group = episode_id, color = "Observed episodes"),
    size = 0.5,
    alpha = 0.3
  ) +
  ylim(0,11)+
  scale_color_manual(values = c("Median simulations"="#6b6b6b","Observed episodes"="#c73535")) +
  scale_fill_manual(values = c("Simulations"="#bdbdbd","Observed episodes"="#e28b8b")) +
  labs(x="Time step", y="Rainfall (mm/5min)", color="", fill="") +
  theme_minimal() +
  btf_theme

names(df_sim_summary)
names(df_real_summary)


# save this plot
ggsave(paste0(im_folder, "swg/omsev/sim_vs_real/adv_group_", group_adv, ".png"),
       width = 8, height = 6, dpi = 300)

X_s0_t0 <- c()
for (i in seq_along(s0_list_group)) {
  s0_i <- s0_list_group[i]
  episode_i <- list_episodes_group[[i]]
  X_s0_t0 <- c(X_s0_t0, episode_i[1, s0_i])
}
X_s0_t0_mean <- mean(X_s0_t0)
adv_mean <- colMeans(adv_group)
sims_from_mean <- vector("list", Nsim)
for (i in seq_len(Nsim)) {
  s0_random <- s0_list_group[sample(seq_along(s0_list_group), 1)]
  sims_from_mean[[i]] <- simulate_many_episodes(
    N = 1,
    u = 1000,
    u_emp = X_s0_t0_mean,   # on prend la moyenne
    params_vario = params_vario,
    params_margins = params_margins,
    coords = grid_omsev,
    times = times,
    adv = adv_mean,
    t0 = 0,
    s0 = s0_random
  )[[1]]
}


# --- Dataframe pour IC95% ---
df_sims <- do.call(rbind, lapply(seq_len(Nsim), function(i) {
  s0_i <- s0_list_group[sample(seq_along(s0_list_group), 1)]
  data.frame(
    time = seq_len(ncol(sims_from_mean[[i]])),
    rainfall = sims_from_mean[[i]][s0_i, ],
    sim_id = i
  )
}))

df_sim_summary <- df_sims %>%
  group_by(time) %>%
  summarise(
    sim_mean = mean(rainfall),
    sim_median = median(rainfall),
    sim_sd  = sd(rainfall),
    sim_low  = quantile(rainfall, 0.025),
    sim_high = quantile(rainfall, 0.975)
  )

# --- Tracé avec épisodes réels ---
df_real_group <- do.call(rbind, lapply(seq_along(list_episodes_group), function(idx) {
  ep <- list_episodes_group[[idx]]
  s0_i <- s0_list_group[idx]
  data.frame(
    time = seq_len(nrow(ep)),
    rainfall = ep[, s0_i],
    episode_id = idx
  )
}))

df_real_summary <- df_real_group %>%
  group_by(time) %>%
  summarise(
    real_mean = mean(rainfall),
    real_median = median(rainfall),
    real_low  = quantile(rainfall, 0.025),
    real_high = quantile(rainfall, 0.975)
  )
# 0.25 to 0.75 quantiles = 50 % CI
# 0.025 to 0.975 quantiles = 95 % CI
# 0.05 to 0.95 quantiles = 90 % CI
# 99 % CI = 0.005 to 0.995 quantiles

ggplot() +
  geom_ribbon(
    data = df_sim_summary,
    aes(x = time, ymin = sim_low, ymax = sim_high, fill = "Simulations"),
    alpha = 0.25
  ) +
  geom_line(
    data = df_sim_summary,
    aes(x = time, y = sim_median, color = "Simulations"),
    size = 1
  ) +
  geom_line(
    data = df_real_summary,
    aes(x = time, y = real_median, color = "Observed episodes"),
    size = 1
  ) +  geom_ribbon(
    data = df_real_summary,
    aes(x = time, ymin = real_low, ymax = real_high, fill = "Observed episodes"),
    alpha = 0.25
  ) +
  scale_color_manual(values = c("Simulations"="#6b6b6b","Observed episodes"="#c73535")) +
  scale_fill_manual(values = c("Simulations"="#bdbdbd","Observed episodes"="#e28b8b")) +
  labs(x="Time step", y="Rainfall (mm/5min)", color="", fill="") +
  theme_minimal() +
  btf_theme




# save this plot
ggsave(paste0(im_folder, "swg/omsev/sim_vs_real/CI_obs95_CI_sim95_adv_group_mean_initmean_", group_adv, ".png"),
       width = 12, height = 6, dpi = 300)


ggplot() +
  geom_ribbon(
    data = df_sim_summary,
    aes(x = time, ymin = sim_low, ymax = sim_high, fill = "Simulations"),
    alpha = 0.25
  ) +
  geom_line(
    data = df_sim_summary,
    aes(x = time, y = sim_median, color = "Median simulations"),
    size = 1
  ) +
  geom_line(
    data = df_real_group,
    aes(x = time, y = rainfall, group = episode_id, color = "Observed episodes"),
    size = 0.5,
    alpha = 0.3
  ) +
  scale_color_manual(values = c("Simulations"="#6b6b6b","Observed episodes"="#c73535")) +
  scale_fill_manual(values = c("Simulations"="#bdbdbd","Observed episodes"="#e28b8b")) +
  labs(x="Time step", y="Rainfall (mm/5min)", color="", fill="") +
  theme_minimal() +
  btf_theme

# save this plot
ggsave(paste0(im_folder, "swg/omsev/sim_vs_real/adv_group_mean_initmean_episodes_", group_adv, ".png"),
       width = 12, height = 6, dpi = 300)

# compute empirical chi vs theoretical chi
ep_id <- 100
s0 <- s0_list[ep_id]
adv_ep <- as.numeric(adv_df_transfo[ep_id, c("vx", "vy")])
times <- 0:11 # in 5 min
times_hour <- times / 12  # in hours
params <- c(params_vario$beta1, params_vario$beta2, params_vario$alpha1, params_vario$alpha2, adv_ep)
s0_coords <- as.numeric(sites_coords[rownames(sites_coords) == s0, ])
df_lags <- get_conditional_lag_vectors(sites_coords, s0_coords, t0=0, tau_vect=times_hour, latlon=TRUE)
chi_th <- theoretical_chi(params, df_lags, latlon=FALSE, distance = "euclidean")

# plot theoretical chi vs hnorm by tau
ggplot(chi_th, aes(x = hnorm, y = chi, color = factor(tau))) +
  geom_point() +
  geom_line() +
  labs(x = "Normalized spatial lag", y = "Theoretical chi", color = "Temporal lag (tau)") +
  theme_minimal()

# compute empirical chi
chi_st_emp <- chispatemp_empirical(rain_omsev, df_lags, quantile =0.95, remove_zeros = TRUE)
head(chi_st_emp)
head(chi_th)

# look at rain_omsev from 24 to 26 november 2021
# rain_episode <- rain_omsev %>%
#   filter(dates >= "2021-11-24 09:00:00" & dates <= "2021-11-26 01:00:00")
# rain_iem <- rain_episode$cefe
# plot(rain_iem, type="l",
#      xlab="Time step (5 min)", ylab="Total rainfall (mm/5min)",
#      main="Total rainfall over all stations from 24 to 26 Nov 2021")

df_chi <- chi_st_emp %>%
  left_join(chi_th %>% select(s1, s2, tau, chi), by = c("s1", "s2", "tau"))
head(df_chi)

ggplot(df_chi, aes(x = chi, y = chiemp2)) +
  geom_point(alpha = 0.6) +        # points
  geom_abline(slope = 1, intercept = 0, color = "red", linetype = "dashed") + # y=x
  labs(x = "Chi théorique", y = "Chi empirique") +
  theme_minimal()

ggplot(df_chi, aes(x = chi, y = chiemp2, color = factor(tau))) +
  geom_point(alpha = 0.6) +
  geom_abline(slope = 1, intercept = 0, color = "black", linetype = "dashed") +
  labs(x = "Chi théorique", y = "Chi empirique", color = "Lag temporel") +
  theme_minimal()

df_chi$h_bin
df_chi <- df_chi %>%
  mutate(h_bin = cut(hnorm,
                     breaks = seq(0, max(hnorm), length.out = 12),
                     include.lowest = TRUE))

ggplot(df_chi, aes(x = chi, y = chiemp2, color = factor(tau))) +
  geom_point(alpha = 0.7) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red", alpha = 0.5) +
  facet_wrap(~h_bin) +
  labs(x = "Theorical chi", y = "Empirical chi") +
  theme_minimal()



# save this plot
ggsave(paste0(im_folder, "swg/omsev/chi_emp_vs_chi_th_adv_group_", group_adv, ".png"),
       width = 12, height = 6, dpi = 300)


library(scoringRules)

# Nombre de pas de temps
n_time <- max(df_real_group$time)

# Pour chaque pas de temps, récupérer toutes les simulations et l'observation (moyenne si plusieurs épisodes)
crps_group <- sapply(seq_len(n_time), function(t) {
  sims_t <- df_all_sims_group %>%
    filter(time == t) %>%
    pull(rainfall)
  
  obs_t <- df_real_group %>%
    filter(time == t) %>%
    pull(rainfall)
  
  # Si plusieurs épisodes réels, prendre la moyenne
  obs_t_mean <- mean(obs_t)
  
  crps_sample(y = obs_t_mean, dat = sims_t)
})

# Résumé du CRPS
mean_crps_group <- mean(crps_group)
median_crps_group <- median(crps_group)

mean_crps_group
median_crps_group


library(scoringRules)
library(dplyr)

n_sims <- length(unique(df_all_sims_group$sim_id))
n_time <- max(df_all_sims_group$time)

crps_within_sims <- matrix(NA, nrow = n_time, ncol = n_sims)

for (i in 1:n_sims) {
  for (t in 1:n_time) {
    sim_i <- df_all_sims_group %>% filter(sim_id == i, time == t) %>% pull(rainfall)
    sims_other <- df_all_sims_group %>% filter(sim_id != i, time == t) %>% pull(rainfall)
    crps_within_sims[t, i] <- crps_sample(y = sim_i, dat = sims_other)
  }
}

# Moyenne sur les simulations et le temps
mean_crps_within <- mean(crps_within_sims)
median_crps_within <- median(crps_within_sims)

mean_crps_within
median_crps_within

n_time <- max(df_real_group$time)

crps_vs_obs <- sapply(seq_len(n_time), function(t) {
  sims_t <- df_all_sims_group %>%
    filter(time == t) %>%
    pull(rainfall)
  
  obs_t <- df_real_group %>%
    filter(time == t) %>%
    pull(rainfall)
  
  # Moyenne des épisodes réels si plusieurs
  obs_t_mean <- mean(obs_t)
  
  crps_sample(y = obs_t_mean, dat = sims_t)
})

# Résumé
mean_crps_vs_obs <- mean(crps_vs_obs)
median_crps_vs_obs <- median(crps_vs_obs)

mean_crps_vs_obs
median_crps_vs_obs


library(ggplot2)
library(dplyr)

# Résumé des simulations
df_summ_group <- df_all_sims_group %>%
  group_by(time) %>%
  summarize(
    median = median(rainfall),
    q025 = quantile(rainfall, 0.025),
    q975 = quantile(rainfall, 0.975)
  )

# Observations moyennes par pas de temps
df_obs_group <- df_real_group %>%
  group_by(time) %>%
  summarize(observed = mean(rainfall))

# Calcul CRPS vs observations
library(scoringRules)
crps_group <- sapply(seq_len(nrow(df_summ_group)), function(t) {
  sims_t <- df_all_sims_group %>% filter(time == t) %>% pull(rainfall)
  obs_t <- df_obs_group$observed[t]
  crps_sample(y = obs_t, dat = sims_t)
})

# Graphique
ggplot() +
  geom_ribbon(data = df_summ_group, aes(x = time, ymin = q025, ymax = q975),
              fill = "#979595", alpha = 0.3) +
  geom_line(data = df_summ_group, aes(x = time, y = median),
            color = "#979595", size = 1) +
  geom_line(data = df_obs_group, aes(x = time, y = observed),
            color = "#c73535", size = 1) +
  geom_line(aes(x = seq_len(nrow(df_summ_group)), y = crps_group*max(df_summ_group$median)/max(crps_group)),
            color = "blue", linetype = "dashed") +
  scale_y_continuous(
    name = "Rainfall (mm/5min)",
    sec.axis = sec_axis(~ . * max(crps_group) / max(df_summ_group$median),
                        name = "CRPS")
  ) +
  labs(
    x = "Time step",
    y = "Rainfall (mm/5min)",
    title = "Simulations vs Observations avec CRPS",
    subtitle = "Ribben: 95% CI simulations, rouge: observations, ligne bleue: CRPS"
  ) +
  theme_minimal()
