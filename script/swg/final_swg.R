rm(list = ls())

# clear console
cat("\014")

# Load libraries and set theme
source("./script/load_libraries.R")

# library(generain)
library(animation)
library(tidyr)
library(lubridate)
library(purrr)
library(viridis)
library(magick)
library(sf)
library(units)
library(ggspatial)
library(dplyr)
library(tidyr)
library(ggplot2)
library(forcats)

# Get all files in the folder "R"
functions_folder <- "./R"
files <- list.files(functions_folder, full.names = TRUE)
# load all functions in files
invisible(lapply(files, function(f) source(f, echo = FALSE)))
library(latex2exp)

################################################################################
# DATA AND COORDS
################################################################################

# get rain data
filename_rain <- paste0(data_folder, "omsev/omsev_5min/rain_mtp_5min_2019_2025.csv")
rain_omsev <- read.csv(filename_rain)
head(rain_omsev)

# egpd fit
filename_egpd <- paste0(data_folder, "../thesis/resources/images/EGPD/OMSEV/2019_2025/egpd_results.csv")
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
grid_coords_km <- sites_coords
grid_coords_m <- sites_coords
grid_coords_m$x_m <- (coords_m[, "X"] - min(coords_m[, "X"]))
grid_coords_m$y_m <- (coords_m[, "Y"] - min(coords_m[, "Y"]))
grid_coords_km$x_km <- (coords_m[, "X"] - min(coords_m[, "X"])) / 1000
grid_coords_km$y_km <- (coords_m[, "Y"] - min(coords_m[, "Y"]))  / 1000

# get distance matrix
grid_coords_m <- grid_coords_m[, c("x_m", "y_m")]
grid_coords_km <- grid_coords_km[, c("x_km", "y_km")]
colnames(grid_coords_m) <- c("Longitude", "Latitude")
colnames(grid_coords_km) <- c("Longitude", "Latitude")


dist_mat <- get_dist_mat(grid_coords_m, latlon = FALSE)

# Spatial chi
df_dist <- reshape_distances(dist_mat)
df_dist_km <- df_dist
df_dist_km$value <- df_dist$value / 1000

rownames(grid_coords_m) <- rownames(sites_coords)


################################################################################
## Marginal parameters
################################################################################

# keep rain rows with at least 50% non NA
# rain_complete <- rain[rowSums(is.na(rain)) <= (ncol(rain) / 2), ]
# head(rain_complete)
# rain <- rain_complete

# to compute p0 remove all cumul over 1h above 10 mm
window_time <- 12 # 1 hour = 12 * 5min
# rain_cumul_1h <- zoo::rollapply(rain, width = window_time, FUN = sum, align = "right", fill = NA)
# rain_cumul_1h[rain_cumul_1h > 22] <- NA

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

# for all sites names in params_margins$p0, put p0 to 0.95
# for (i in seq_along(params_margins$p0)) {
#   site <- names(params_margins$p0)[i]
#   if (site %in% colnames(rain)) {
#     params_margins$p0[i] <- 0.95
#   }
# }

#################################################################################
# GET EPISODES
#################################################################################
# rain_complete <- rain[rowSums(is.na(rain)) <= (ncol(rain) / 2), ]
# head(rain_complete)
# rain <- rain_complete
q <- 0.95
delta <- 12
dmin <- 1200  # m
times <- 0:(delta - 1)

# in rain remove when all data are NA
set_st_excess <- get_spatiotemp_excess(rain, quantile = q, remove_zeros = TRUE)

# Spatio-temporal neighborhood parameters
s0t0_set <- get_s0t0_pairs(grid_coords_m, rain,
                            min_spatial_dist = dmin,
                            episode_size = delta,
                            set_st_excess = set_st_excess,
                            n_max_episodes = 10000,
                            latlon = FALSE)
selected_points <- s0t0_set
selected_points <- selected_points %>%
  mutate(t0_date = as.POSIXct(t0_date, format = "%Y-%m-%d %H:%M:%S", tz = "UTC"))

n_episodes <- length(selected_points$s0)
table(selected_points$s0) # number of episodes per s0

# do a plot for the number of episodes per site
ggplot(selected_points, aes(x = s0)) +
  geom_bar(fill = btfgreen, alpha = 0.7) +
  labs(
    x = "Site (s0)",
    y = "Number of episodes"
  ) +
  theme_minimal() 

# save plot
filename_plot <- paste0(
  im_folder,
  "swg/omsev/2025/number_episodes_per_site_q",
  q * 100,
  "dmin",
  dmin,
  "delta",
  delta,
  ".png"
)
ggsave(
  filename = filename_plot,
  plot = last_plot(),
  width = 7,
  height = 5,
  dpi = 300
)

# number of episodes par t0 date
n_s0_per_t0 <- table(selected_points$t0_date)

t0_list <- selected_points$t0
s0_list <- selected_points$s0
u_list <- selected_points$u_s0

# get extreme episodes
list_episodes_points <- get_extreme_episodes(selected_points, rain,
                                     episode_size = delta, unif = FALSE)
list_episodes <- list_episodes_points$episodes
tau_vect <- 0:10
tmax <- max(tau_vect)
df_coords <- as.data.frame(grid_coords_m)
# Compute the lags and excesses for each conditional point
list_results <- mclapply(1:length(s0_list), function(i) {
  s0 <- s0_list[i]
  row_s0 <- which(rownames(df_coords) == s0)
  s0_coords <- df_coords[row_s0, ]
  episode <- list_episodes[[i]]
  ind_t0_ep <- 0 # index of t0 in the episode
  lags <- get_conditional_lag_vectors(df_coords, s0_coords, ind_t0_ep,
                                tau_vect, latlon = FALSE)
  u <- u_list[i]
  excesses <- empirical_excesses_rpar(episode, threshold = u,
                                      df_lags = lags, t0 = ind_t0_ep)
  # tau is in 5 minutes
  lags$tau <- lags$tau * 5 / 60 # convert to hours
  # lags$hnorm <- lags$hnorm / 1000 # convert to km
  list(lags = lags, excesses = excesses)
}, mc.cores = detectCores() - 1)

list_lags <- lapply(list_results, `[[`, "lags")
list_excesses <- lapply(list_results, `[[`, "excesses")
df_lags <- list_lags[[10]] # km/h
df_excesses <- list_excesses[[13]]
sum(df_excesses$kij)

length(list_episodes)

###################################################################################
# VARIOS PARAMETERS FROM KM/H TO M/5MIN
###################################################################################
# From results
# params_est <- c(1.2953654, 4.2212009, 0.2495586, 0.6657619, 3.8962660, 2.2208320)
params_est <-  c(1.101,3.716,0.103,0.676, 3.896, 2.221) #2025
etas_estimates <- params_est[5:6]

params_kmh <- list(
  beta1 = params_est[1],
  beta2 = params_est[2],
  alpha1 = params_est[3],
  alpha2 = params_est[4]
)

# convert params from km/h to m/5min
c_x_m <- 1000    # for m
c_t_5min <- 12   # 1 hour = 12 * 5min

# convert params and ci to m/5min
params_m5min_beta <- convert_params(params_kmh$beta1, params_kmh$beta2, params_kmh$alpha1, params_kmh$alpha2,
                               c_x = c_x_m, c_t = c_t_5min)

params_m5min <- list(
  beta1 = params_m5min_beta$beta1,
  beta2 = params_m5min_beta$beta2,
  alpha1 =  params_est[3],
  alpha2 = params_est[4]
)

beta1 <- params_m5min$beta1
beta2 <- params_m5min$beta2
alpha1 <- params_m5min$alpha1
alpha2 <- params_m5min$alpha2

##################################################################################
# TRANSFORM ADVECTION SPEEDS
##################################################################################
adv_filename <- paste(data_folder, "/omsev/adv_estim/2025/combined_comephore_omsev/episode_advection_q",
                          q * 100, "_delta", delta, "_dmin", dmin,
                          ".csv", sep = "")
adv_df_raw <- read.csv(adv_filename, sep = ",")
head(adv_df_raw)
setDT(selected_points)
setDT(adv_df_raw)

# is there na in adv_df_raw?
anyNA(adv_df_raw)

adv_df_raw[, t0_omsev := as.POSIXct(t0_omsev, format="%Y-%m-%d %H:%M:%S", tz="UTC")]
selected_points[, t0_date := as.POSIXct(t0_date, format="%Y-%m-%d %H:%M:%S", tz="UTC")]

adv_df_t0 <- adv_df_raw[, .(
  vx_final = vx_final[1],
  vy_final = vy_final[1]
), by = t0_omsev]

# is there na in adv_df_t0?
anyNA(adv_df_t0)
setkey(adv_df_t0, t0_omsev)

# get advections for selected points
selected_episodes <- adv_df_t0[selected_points, on = .(t0_omsev = t0_date)]
setnames(selected_episodes, c("vx_final","vy_final"), c("adv_x","adv_y"))
anyNA(selected_episodes)
which(is.na(selected_episodes$adv_x))
selected_episodes[is.na(selected_episodes$adv_x), ]
# Advection
V_episodes <- data.frame(
  v_x = selected_episodes$adv_x,
  v_y = selected_episodes$adv_y
)
adv_df <- V_episodes
colnames(adv_df) <- c("vx_final", "vy_final")
adv_df_transfo <- adv_df

# Transform advections with etas estimates
adv_df_transfo$vnorm <- sqrt(adv_df$vx_final^2 + adv_df$vy_final^2)
adv_df_transfo$vnorm_t <- etas_estimates[1] * adv_df_transfo$vnorm^etas_estimates[2]

adv_df_transfo$vx_t <- ifelse(
  adv_df_transfo$vnorm > 0,
  adv_df_transfo$vx_final / adv_df_transfo$vnorm * adv_df_transfo$vnorm_t,
  0
)

adv_df_transfo$vy_t <- ifelse(
  adv_df_transfo$vnorm > 0,
  adv_df_transfo$vy_final / adv_df_transfo$vnorm * adv_df_transfo$vnorm_t,
  0
)

plot(adv_df$vx_final, adv_df$vy_final, pch = 19,
     xlab = "vx", ylab = "vy")

plot(adv_df_transfo$vx_t, adv_df_transfo$vy_t, pch = 19,
     xlab = "vx", ylab = "vy")

# if too big advections after transformation, cap them
# max_adv <- 100  # km/h
# adv_df_transfo$vx_t <- pmax(pmin(adv_df_transfo$vx_t, max_adv), -max_adv)
# adv_df_transfo$vy_t <- pmax(pmin(adv_df_transfo$vy_t, max_adv), -max_adv)

# Compute speed
speed_raw <- sqrt(adv_df$vx_final^2 + adv_df$vy_final^2)
summary(speed_raw)
speed_transfo <- sqrt(adv_df_transfo$vx_t^2 + adv_df_transfo$vy_t^2)
summary(speed_transfo)

hist(speed_raw, breaks = 20, main = "Histogram of raw advection speeds", xlab = "Speed (km/h)")
hist(speed_transfo, breaks = 20, main = "Histogram of transformed advection speeds", xlab = "Speed (km/h)")
# Compute direction
direction_raw <- atan2(adv_df$vy_final, adv_df$vx_final) * (180 / pi)
# make sure direction is in [0, 360]
direction_transfo <- atan2(adv_df_transfo$vy_t, adv_df_transfo$vx_t) * (180 / pi)

summary(direction_raw)
summary(direction_transfo)
# same ok
direction <- direction_transfo
# make sure direction is in [0, 360]
angle_deg <- (direction + 360) %% 360
is_calm <- speed_transfo <= 1e-3
speed_class <- as.character(cut(
  speed_transfo,
  breaks = c(0, 1e-2, 0.5, 5, 170, Inf),
  labels = c("null", "still","weak","significant", "invalid"),
  include.lowest = TRUE
  ))
speed_class <- factor(speed_class, levels = c("null","still","weak","significant", "invalid"))
table(speed_class)

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

adv_class <- data.frame(
  speed = speed_raw,
  speed_transfo = speed_transfo,
  direction = direction,
  angle_deg = angle_deg,
  direction_class = direction_class,
  speed_class = speed_class,
  group = paste0(speed_class, "_", direction_class)
)
# if speed_class is invalid, set group to invalid
adv_class$group <- ifelse(adv_class$speed_class == "invalid", "invalid", adv_class$group)
head(adv_class)
table(adv_class$group)
# number of null advections
n_null_adv <- sum(adv_class$speed_class == "null")
# Wind rose plot
wind_rose_data <- adv_class %>%
  filter(direction_class != "none") %>%
  filter(speed_class != "invalid") %>%
  group_by(direction_class, speed_class) %>%
  summarise(count = n(), .groups = "drop")

# remove null advections from wind rose data
wind_rose_data <- wind_rose_data %>%
  filter(speed_class != "null")
ggplot(wind_rose_data, aes(x = direction_class, y = count, fill = speed_class)) +
  geom_col(position = "stack", width = 0.7) +
  scale_fill_manual(values = c("still" = "#fee5d9", "weak" = "#fcae91", "significant" = "#e6550d")) +
  labs(
  title = paste0("Number of null advections without direction: ", n_null_adv),
  x = "Direction",
  y = "Episodes count (advective only)",
  fill = "Speed class"
  ) +
  theme_minimal() +
  coord_polar(theta = "x", start = pi/4, direction = 1)+
  btf_theme

# save plot
foldername_plot <- paste0(im_folder,"swg/omsev/2025/")
filename_plot <- paste0(
  foldername_plot, "adv_wind_rose_95q1200dmin12delta_NSEW.png")

ggsave(
  filename = filename_plot,
  plot = last_plot(),
  width = 7,
  height = 7,
  dpi = 300
)

table(adv_class$group)

selected_points$speed_class <- adv_class$speed_class
selected_points$adv_group <- adv_class$group

# if needed, get X_s0_t0 for each episode
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

################################################################################
# Simulation for similar episodes advection
################################################################################
# for all sites names in params_margins$p0, put p0 to 0.95
# for (i in seq_along(params_margins$p0)) {
#   site <- names(params_margins$p0)[i]
#   if (site %in% colnames(rain)) {
#     params_margins$p0[i] <- 0.99
#   }
# }
grid_omsev <- grid_coords_m
adv_matrix <- as.matrix(adv_df_transfo[, c("vx_t", "vy_t")])
all_group_names <- unique(selected_points$speed_class)
# all_significant_groups <- all_group_names[grep("significant", all_group_names)]
group_adv <- all_group_names[4]  # choose one group to simulate
# get list_episodes and s0_list for this group
indices_group <- which(selected_points$speed_class == group_adv)
list_episodes_group <- list_episodes[indices_group]

# cum_1h <- sapply(list_episodes_group, function(ep) sum(ep, na.rm = TRUE))
# summary(cum_1h)
# thr_1h <- 550/24  # ~22.9
# list_episodes_keep <- list_episodes_group
# length(list_episodes_keep)
# indices_group_keep <- indices_group

# p0_wet <- sapply(seq_along(sites_names), function(j) {
#   x <- unlist(lapply(list_episodes_keep, function(ep) ep[, j]))
#   mean(x == 0, na.rm = TRUE)
# })
# names(p0_wet) <- sites_names
# summary(p0_wet)

# params_margins_wet <- list(
#   xi    = xi_vect,
#   sigma = sigma_vect,
#   kappa = kappa_vect,
#   p0    = p0_wet
# )

list_episodes_keep <- list_episodes

length(list_episodes_keep)
s0_list_keep <- s0_list
u_list_keep <- u_list
adv_group <- adv_matrix
Nsim <- 100
s0_sim <- integer(Nsim)
sims_group <- vector("list", Nsim)
adv_sim <- matrix(0, nrow = Nsim, ncol = 2)
idx <- sample(seq_along(list_episodes_keep), 1)
for (i in seq_len(Nsim)) {
  # idx <- sample(seq_along(list_episodes_keep), 1)
  # idx <- 20
  s0_i <- s0_list_keep[idx]
  s0_sim[i] <- s0_i
  u_i <- u_list_keep[idx]
  adv_i <- adv_group[idx, ]
  adv_sim[i, ] <- adv_i
  sims_group[[i]] <- simulate_many_episodes(
    N = 1,
    u_emp = u_i,
    params_vario = params_m5min,
    params_margins = params_margins,
    coords = grid_omsev, # km
    times = times,  # in hours
    adv = adv_i,
    t0 = 0,
    s0 = s0_i
  )[[1]]
}


# for one simulation, get cumul by sites
sim1 <- sims_group[[1]]$X
head(sim1)
cumul_sim1 <- colSums(sim1)

# for all simulations, get cumul by sites by simulation
cumul_sims <- matrix(0, nrow=Nsim, ncol=ncol(sim1))
sites_names <- colnames(sim1)
colnames(cumul_sims) <- sites_names
p <- 0.22

for (i in seq_len(Nsim)) {
  sim_i <- sims_group[[i]]$X
  
  # threshold <- 0.2*p
  sim_i[sim_i < 1e-3] <- 0
  # sim_i_discrete <- round(sim_i / p) * p
  cumul_sim_i <- colSums(sim_i)
  cumul_sims[i, ] <- cumul_sim_i
}

head(cumul_sims)
colSums(cumul_sims)

# keep only episode with enough activated sites ie 70% non na
list_episodes_keep  <- list_episodes_group[sapply(list_episodes_group, function(ep) {
  mean(!is.na(ep)) >= 0.5
})]
length(list_episodes_keep)
ep1 <- list_episodes[[1]]
cumul_eps <- matrix(0, nrow=length(list_episodes_keep), ncol=ncol(ep1))
colnames(cumul_eps) <- sites_names
for (i in seq_len(length(list_episodes_keep))) {
  ep_i <- list_episodes_keep[[i]]
    ep_i[ep_i < p*2] <- 0

  cumul_ep_i <- colSums(ep_i, na.rm = TRUE)
  cumul_eps[i, ] <- cumul_ep_i
}
head(cumul_eps)
mean_cumul_eps <- colMeans(cumul_eps)
med_cumul_eps <- 

par(mfrow=c(1,2))
boxplot(
  cumul_sims,
  xlab = "Site",
  ylab = "Cumulative rainfall",
  main = "Simulated cumul (boxplots) with observed mean per site",
  outline = FALSE
)

# points(
#   x = seq_along(cumul_eps),
#   y = cumul_eps,
#   pch = 19,
#   col = "red",
#   cex = 1.2
# )

boxplot(
  cumul_eps,
  xlab = "Site",
  ylab = "Cumulative rainfall",
  main = "Simulated cumul (boxplots) with observed mean per site",
  outline = FALSE
)

# do cumul by simulation and by episode, then compare them with a plot
cum_sim_mean <- rowSums(cumul_sims)

mean_cumul_eps_total <- sum(mean_cumul_eps)

cum_obs <- rowSums(cumul_eps)

boxplot(
  rowSums(cumul_sims),
  xlab = "Simulations",
  ylab = "Total cumulative rainfall",
  main = "Total cumulative rainfall by simulation",
  ylim = c(0, 500)
)



boxplot(cum_obs,
        xlab = "Observed episodes",
        ylab = "Total cumulative rainfall",
        main = "Total cumulative rainfall by observed episode",
        ylim = c(0, 500)
)

par(mfrow=c(1,1))
boxplot(list(Observed = cum_obs, Simulations = cum_sim_mean),
        ylab = "Total cumulative rainfall",
        ylim = c(0, 500))




# compare cumul by sites for sim and real episodes
cumul_comparison <- data.frame(
  site = sites_names,
  cumul_sim_mean = colMeans(cumul_sims),
  cumul_sim_sd = apply(cumul_sims, 2, sd),
  cumul_ep_mean = colMeans(cumul_eps),
  cumul_ep_sd = apply(cumul_eps, 2, sd)
)
n_obs_site <- colSums(!is.na(cumul_eps))

cumul_comparison <- data.frame(
  site = sites_names,
  n_obs = n_obs_site,
  cumul_sim_mean = colMeans(cumul_sims),
  cumul_sim_sd = apply(cumul_sims, 2, sd),
  cumul_ep_mean = colSums(cumul_eps) / n_obs_site,
  cumul_ep_sd = apply(cumul_eps, 2, sd, na.rm = TRUE)
)


cumul_long <- cumul_comparison %>%
  pivot_longer(cols = c(cumul_sim_mean, cumul_ep_mean), names_to = "type", values_to = "mean") %>%
  pivot_longer(cols = c(cumul_sim_sd, cumul_ep_sd), names_to = "type_sd", values_to = "sd") %>%
  filter((type == "cumul_sim_mean" & type_sd == "cumul_sim_sd") |
           (type == "cumul_ep_mean" & type_sd == "cumul_ep_sd")) %>%
  mutate(type = ifelse(type == "cumul_sim_mean", "Simulations", "Observations")) %>%
  select(site, type, mean, sd)
ggplot(cumul_long, aes(x = site, y = mean, fill = type)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7, alpha = 0.85) +
  geom_errorbar(
    aes(ymin = pmax(0, mean - sd), ymax = mean + sd),
    position = position_dodge(width = 0.8),
    width = 0.2
  ) +
  labs(
    x = "Site", y = "Cumulative rainfall (mm)",
    fill = NULL
  ) +
  theme_minimal() +
  btf_theme +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


# cumul by episode 
cumul_by_episode <- data.frame(
  episode_id = seq_len(nrow(cumul_eps)),
  cumul_ep = rowSums(cumul_eps)
)

cumul_by_episode_sim <- data.frame(
  sim_id = seq_len(nrow(cumul_sims)),
  cumul_sim = rowSums(cumul_sims)
)

cumul_comparison_episode <- data.frame(
  type = c(rep("Observed episodes", nrow(cumul_by_episode)), rep("Simulations", nrow(cumul_by_episode_sim))),
  cumul = c(cumul_by_episode$cumul_ep, cumul_by_episode_sim$cumul_sim)
)
min(cumul_comparison_episode$cumul[cumul_comparison_episode$type == "Observed episodes"])
min(cumul_comparison_episode$cumul)
# remove cumul < 1
cumul_comparison_episode <- cumul_comparison_episode[cumul_comparison_episode$cumul >= 12, ]
ggplot(cumul_comparison_episode, aes(x = type, y = cumul, fill = type)) +
  geom_boxplot(alpha = 0.7) +
  labs(x = "Type", y = "Total cumulative rainfall (mm)", title = "Total cumulative rainfall: Simulations vs Episodes") +
  theme_minimal() +
  theme(legend.position = "none") +
  btf_theme

# violin plot
ggplot(cumul_comparison_episode, aes(x = type, y = cumul, fill = type)) +
  geom_violin(alpha = 0.7, scale = "width") +
  labs(x = "Type", y = "Total cumulative rainfall (mm)", title = "Total cumulative rainfall: Simulations vs Episodes") +
  theme_minimal() +
  theme(legend.position = "none") +
  btf_theme







# check that we have exceedances at s0, t0 for each episode
u_sim <- sapply(seq_len(Nsim), function(i) u_list_group[which(s0_list_group == s0_sim[i])[1]])

all(sapply(seq_len(Nsim), function(i) {
  sim_i <- sims_group[[i]]
  s0_i <- s0_sim[i]
  u_i  <- u_sim[i]
  sim_i$X[1, s0_i] > u_i
}))


df_sims <- do.call(rbind, lapply(seq_len(Nsim), function(i) {
  s0_i <- s0_sim[i]
  sim_i <- sims_group[[i]]
  data.frame(
    time = seq_len(nrow(sim_i$X)),
    rainfall = sim_i$X[, s0_i],
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


ggplot() +
  geom_ribbon(
    data = df_sim_summary,
    aes(x = time, ymin = sim_low, ymax = sim_high, fill = "Simulations"),
    alpha = 0.30
  ) +

  geom_line(
    data = df_sim_summary,
    aes(x = time, y = sim_median, color = "Simulations"),
    linewidth = 1.1
  ) +
  geom_line(
    data = df_real_group,
    aes(x = time, y = rainfall, group = episode_id, color = "Observed episodes"),
    linewidth = 0.5,
    alpha = 0.25
  ) +

  ylim(0, 11) +

  scale_color_manual(
    values = c(
      "Simulations" = "#5c5c5c",
      "Observed episodes" = "#c73535"
    )
  ) +
  scale_fill_manual(
    values = c(
      "Simulations" = "#bdbdbd",
      "Observed episodes" = "#e28b8b"
    )
  ) +

  labs(
    x = "Time step",
    y = "Rainfall (mm / 5 min)",
    color = "",
    fill  = ""
  ) +

  theme_minimal() +
  btf_theme


# save plot
# foldername_plot <- paste0(
#   im_folder,
#   "swg/omsev/all/new/"
# )
# if (!dir.exists(foldername_plot)) {
#   dir.create(foldername_plot, recursive = TRUE)
# }
# filename <- paste0(
#   foldername_plot,
#   "simulated_vs_observed_rainfall_advgroup_",
#   group_adv,
#   "_n", Nsim,
#   ".png"
# )

# ggsave(
#   filename = filename,
#   plot = last_plot(),
#   width = 7,
#   height = 5,
#   dpi = 300
# )

# plot just simulated rainfall at s0
ggplot(df_sims, aes(x = time, y = rainfall, group = sim_id)) +
  geom_line(alpha = 0.7, color = "gray") +
  labs(x="Time step (5min)", y="Rainfall (mm/5min)") +
  theme_minimal() +
  btf_theme

# save plot
filename_simonly <- paste0(
  foldername_plot,
  "simulated_rainfall_only_advgroup_",
  group_adv,
  "_n", Nsim,
  ".png"
)
ggsave(
  filename = filename_simonly,
  plot = last_plot(),
  width = 8,
  height = 6,
  dpi = 300
)

#################################################################################################################
# on simulated episodes
sims_all <- sims_group
list_results_sim <- lapply(seq_along(sims_all), function(i) {

  Z_ep <- sims_all[[i]]$Z
  s0 <- s0_sim[i]
  row_s0 <- which(rownames(df_coords) == s0)
  s0_coords <- df_coords[row_s0, ]
  ind_t0_ep <- 0

  u_0 <- sims_all[[i]]$u_latent[s0]

  lags <- get_conditional_lag_vectors(
    df_coords, s0_coords, ind_t0_ep,
    tau_vect, latlon = FALSE
  )

  excesses <- empirical_excesses_rpar(
    data_rain   = Z_ep,
    threshold = u_0,
    df_lags   = lags,
    t0        = ind_t0_ep
  )

  lags$tau <- lags$tau * 5/60 # convert to hours
  excesses$tau <- excesses$tau * 5/60 # convert to hours
  lags$hnorm <- lags$hnorm / 1000 # convert to km
  if(sum(excesses$kij) == 0) {
    cat("No exceedances for simulated episode", i, "\n")
  }
  list(lags = lags, excesses = excesses)
})


list_lags_sim <- lapply(list_results_sim, `[[`, "lags")
list_excesses_sim <- lapply(list_results_sim, `[[`, "excesses")
lags_i <- list_lags_sim[[10]]
summary(lags_i$hnorm)
summary(lags_i$tau)
adv_i <- adv_sim[10, ]
adv_x <- adv_i[1]
adv_y <- adv_i[2]
sqrt(adv_x^2 + adv_y^2) * max(lags_i$tau)
excesses_i <- list_excesses_sim[[10]]
sum(excesses_i$kij)
head(excesses_i)
list_sumkij <- sapply(list_excesses_sim, function(df) sum(df$kij))
summary(list_sumkij)
table(list_sumkij > 0)


adv_df_sim <- as.data.frame(adv_sim)
colnames(adv_df_sim) <- c("vx", "vy")

dist_mat <- get_dist_mat(location_gauges) / 1000
df_dist <- reshape_distances(dist_mat)

h_used <- unlist(lapply(list_lags_sim, function(L) L$hnorm))
h_used <- h_used[is.finite(h_used)]
num_intervals <- 12
h_breaks <- quantile(h_used, probs = seq(0, 1, length.out = num_intervals + 1), na.rm = TRUE)
h_breaks <- unique(as.numeric(h_breaks))

h_breaks <- sort(h_breaks)


h_breaks[length(h_breaks)] <- 1.6


params_estimates <- as.numeric(c(params_kmh, etas_estimates))
params_estimates[1] <- params_estimates[1] / 2
params_estimates[3] <- params_estimates[3] / 2
p_sim <- plot_th_emp_chi(
  list_lags = list_lags_sim,
  list_excesses = list_excesses_sim,
  list_adv = adv_df_sim,
  params_estimates = params_estimates,
  tau_min = 0,
  h_breaks = h_breaks,
  latlon = FALSE,
  adv_transform = FALSE
)


# save plot
foldername_plot <- paste0(
  im_folder,
  "swg/omsev/all/new/"
)

if (!dir.exists(foldername_plot)) {
  dir.create(foldername_plot, recursive = TRUE)
}

filename_plot <- paste0(
  foldername_plot,
  "chi_vs_h_bytau_simulated_advgroup_",
  group_adv,
  "_n", Nsim,
  ".png"
)

# plot by tau
res_cmp_sim <- p_sim$res_cmp
res <- p_sim$res
# chi theoretical vs chi empirical with y=x line
ggplot(res_cmp_sim, aes(x = chi_theo_bar, y = chi_emp_bar, size = n_pairs)) +
  geom_point(alpha = 0.7, color=btfgreen) +
  geom_abline(slope = 1, intercept = 0, color = "red", linetype = "dashed") +
  labs(
    x = "Chi theoretical",
    y = "Chi empirical",
    size = "Number of pairs"
  ) +
  theme_minimal()

ggsave(
  filename = filename_plot,
  width = 7,
  height = 7,
  dpi = 300
)

library(boot)
bootstrap_ci <- function(data, n_boot = 1000, conf_level = 0.95) {
  boot_stat <- function(data, indices) {
    d <- data[indices]
    return(mean(d, na.rm = TRUE))
  }
  
  boot_res <- boot(data, statistic = boot_stat, R = n_boot)
  ci <- boot.ci(boot_res, type = "perc", conf = conf_level)
  return(c(ci$percent[4], ci$percent[5]))
}

res_boot <- res %>%
  dplyr::filter(.data$tau >= 0) %>%
  dplyr::group_by(.data$tau, .data$hbin) %>%
  dplyr::summarise(
    chi_emp_bar  = mean(.data$chi_emp, na.rm = TRUE),
    chi_theo_bar = mean(.data$chi_theo, na.rm = TRUE),
    ci_low  = bootstrap_ci(.data$chi_emp)[1],
    ci_high = bootstrap_ci(.data$chi_emp)[2],
    n_pairs = dplyr::n(),
    .groups = "drop"
  )


plot_all <- ggplot(res_boot, aes(
  x = chi_theo_bar,
  y = chi_emp_bar
)) +
  geom_errorbar(
    aes(ymin = ci_low, ymax = ci_high),
    width = 0,
    alpha = 0.4,
    color = btfgreen
  ) +
  geom_point(alpha = 0.7, color = btfgreen) +
  geom_abline(slope = 1, intercept = 0,
              color = "red", linetype = "dashed") +
  labs(
    title = "",
    x = "Theoretical Chi",
    y = "Empirical Chi"
  ) +
  theme_minimal() +
  btf_theme

plot_all

filename_plot <- paste0(
  foldername_plot,
  "chi_vs_h_bytau_simulated_advgroup_CI_",
  group_adv,
  "_n", Nsim,
  ".png"
)

ggsave(
  filename = filename_plot,
  plot = plot_all,
  width = 7,
  height = 7,
  dpi = 300
)
