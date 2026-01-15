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
library(sf)
library(dplyr)
library(ggplot2)
library(units)



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

stopifnot(nrow(grid_coords_m) == ncol(rain))
rownames(grid_coords_m) <- colnames(rain)

################################################################################
## Marginal parameters
################################################################################

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

#######################################################################################################
# GET EPISODES
#######################################################################################################

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
# number of episodes par t0 date
n_s0_per_t0 <- table(selected_points$t0_date)
# see the number of episodes per date in a figure

df_mult <- as.data.frame(n_s0_per_t0) %>%
    dplyr::filter(Freq > 0) %>%
    dplyr::rename(n_s0 = Freq)

ggplot(df_mult, aes(x = factor(n_s0))) +
  geom_bar(fill = btfgreen, alpha = 0.7) +
  labs(
    x = TeX("Number of episodes for the same $t_0$"),
    y = TeX("Number of $t_0$"),
    title = ""
  ) +
  theme_minimal() +
  btf_theme

# save plot
foldername_plot <- paste0(
  im_folder,
  "swg/omsev/"
)

filename_plot <- paste0(
  foldername_plot,
  "n_episodes_per_t0_",
  q * 100, "_dmin", dmin,
    "_delta", delta, ".png"
)

ggsave(
  filename = filename_plot,
  plot = last_plot(),
  width = 7,
  height = 5,
  dpi = 300
)

t0_list <- selected_points$t0
s0_list <- selected_points$s0
u_list <- selected_points$u_s0
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


#######################################################################################################
# VARIOS PARAMETERS FROM KM/H TO M/5MIN
#######################################################################################################
params_est <- c(1.090, 4.628, 0.225, 0.713, 1.621, 5.219)  # km/h
# params_est <- c(1.5, 4.4, 0.2, 0.7, 1.621, 5.219)  # km/h
params_est <- omsev_params
# [1] 1.0871823 4.6266826 0.2238579 0.7127469 1.6210000 5.2190000
etas_estimates <- params_est[5:6]

params_kmh <- list(
  beta1 = params_est[1],
  beta2 = params_est[2],
  alpha1 = params_est[3],
  alpha2 = params_est[4]
)

beta1 <- params_kmh$beta1
beta2 <- params_kmh$beta2
alpha1 <- params_kmh$alpha1
alpha2 <- params_kmh$alpha2


##################################################################################
# TRANSFORM ADVECTION SPEEDS
##################################################################################
adv_filename <- paste(data_folder, "/omsev/adv_estim/combined_comephore_omsev/episode_advection_q",
                          q * 100, "_delta", delta, "_dmin", dmin,
                          ".csv", sep = "")
adv_df_raw <- read.csv(adv_filename, sep = ",")
head(adv_df_raw)
setDT(selected_points)
setDT(adv_df_raw)

adv_df_raw[, t0_omsev := as.POSIXct(t0_omsev, format="%Y-%m-%d %H:%M:%S", tz="UTC")]
selected_points[, t0_date := as.POSIXct(t0_date, format="%Y-%m-%d %H:%M:%S", tz="UTC")]

adv_df_t0 <- adv_df_raw[, .(
  vx_final = vx_final[1],
  vy_final = vy_final[1]
), by = t0_omsev]

setkey(adv_df_t0, t0_omsev)

selected_episodes <- adv_df_t0[selected_points,
  on   = .(t0_omsev = t0_date),
  roll = 5*60
]

setnames(selected_episodes, c("vx_final","vy_final"), c("adv_x","adv_y"))

stopifnot(nrow(selected_episodes) == nrow(selected_points))


V_episodes <- data.frame(
  v_x = selected_episodes$adv_x,
  v_y = selected_episodes$adv_y
)


adv_df <- V_episodes

colnames(adv_df) <- c("vx_final", "vy_final")
adv_df_transfo <- adv_df
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

# Compute speed
speed_raw <- sqrt(adv_df$vx_final^2 + adv_df$vy_final^2)
summary(speed_raw)
speed_transfo <- sqrt(adv_df_transfo$vx_t^2 + adv_df_transfo$vy_t^2)
summary(speed_transfo)


direction <- atan2(adv_df$vy_final, adv_df$vx_final) * (180 / pi)
# make sure direction is in [0, 360]
direction <- atan2(adv_df$vy_final, adv_df$vx_final) * (180 / pi)
# make sure direction is in [0, 360]
angle_deg <- (direction + 360) %% 360
is_calm <- speed_raw == 0
speed_class <- as.character(cut(
  speed_raw,
  breaks = c(0, 0.5, 1, Inf),
  labels = c("still","weak","significant"),
  include.lowest = TRUE
  ))

speed_class <- factor(speed_class, levels = c("still","weak","significant"))
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


# number of null advections
n_null_adv <- sum(adv_class$direction_class == "none")
# Wind rose plot
wind_rose_data <- adv_class %>%
  filter(direction_class != "none") %>%
  group_by(direction_class, speed_class) %>%
  summarise(count = n(), .groups = "drop")

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
  coord_polar(theta = "x", start = pi/4, direction = 1)

# save plot
foldername_plot <- paste0(
  im_folder,
  "swg/omsev/"
)

filename_plot <- paste0(
  foldername_plot,
  "adv_wind_rose_95q1200dmin12delta_NSEW.png"
)

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
grid_omsev <- grid_coords_km
adv_matrix <- as.matrix(adv_df_transfo[, c("vx_t", "vy_t")])
all_group_names <- unique(selected_points$adv_group)
all_significant_groups <- all_group_names[grep("significant", all_group_names)]
group_adv <- all_significant_groups[3]  # choose one group to simulate
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
adv_sim <- matrix(0, nrow = Nsim, ncol = 2)
for (i in seq_len(Nsim)) {
  idx <- sample(seq_along(list_episodes_group), 1)
  s0_i <- s0_list_group[idx]
  s0_sim[i] <- s0_i
  u_i <- u_list_group[idx]
  adv_i <- adv_group[idx, ]
  adv_sim[i, ] <- adv_i
  sims_group[[i]] <- simulate_many_episodes(
    N = 1,
    u = 1000,
    u_emp = u_i,
    params_vario = params_kmh,
    params_margins = params_margins,
    coords = grid_omsev, # km
    times = times * 5 / 60,  # in hours
    adv = adv_i,
    t0 = 0,
    s0 = s0_i
  )[[1]]
}


# check that we have exceedances at s0, t0 for each episode
u_sim <- sapply(sims_group, function(sim) sim$u_latent)
all(sapply(seq_len(Nsim), function(i) {
  sims_group[[i]]$Z[s0_sim[i], 1] >= u_sim[i]
}))


df_sims <- do.call(rbind, lapply(seq_len(Nsim), function(i) {
  s0_i <- s0_sim[i]
  data.frame(
    time = seq_len(ncol(sims_group[[i]]$X)),
    rainfall = sims_group[[i]]$X[s0_i, ],
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
foldername_plot <- paste0(
  im_folder,
  "swg/omsev/all/new/"
)
if (!dir.exists(foldername_plot)) {
  dir.create(foldername_plot, recursive = TRUE)
}
filename <- paste0(
  foldername_plot,
  "simulated_vs_observed_rainfall_advgroup_",
  group_adv,
  "_n", Nsim,
  ".png"
)

ggsave(
  filename = filename,
  plot = last_plot(),
  width = 7,
  height = 5,
  dpi = 300
)

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
list_results_sim <- lapply(seq_along(sims_all), function(i) {

  Z_ep <- sims_all[[i]]$Z
  s0 <- s0_sim[i]
  row_s0 <- which(rownames(df_coords) == s0)
  s0_coords <- df_coords[row_s0, ]
  ind_t0_ep <- 0

  u_0 <- sims_all[[i]]$u_latent

  lags <- get_conditional_lag_vectors(
    df_coords, s0_coords, ind_t0_ep,
    tau_vect, latlon = FALSE
  )

  excesses <- empirical_excesses_rpar(
    data_rain   = t(Z_ep),
    threshold = u_0,
    df_lags   = lags,
    t0        = ind_t0_ep
  )

  lags$tau <- lags$tau * 5/60 # convert to hours
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
if (length(h_breaks) < 3) stop("Not enough unique breaks. Reduce num_intervals.")

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


res <- p_sim$res
res <- res %>%
  mutate(
    adv_norm = sqrt(adv_x^2 + adv_y^2),
    # projection de h sur la direction du vent (km)
    h_par = ifelse(adv_norm > 0, (hx*adv_x + hy*adv_y)/adv_norm, NA_real_),
    # composante transverse
    h_perp = ifelse(adv_norm > 0, sqrt(pmax(h^2 - h_par^2, 0)), NA_real_),
    flow = case_when(
      is.na(h_par) ~ "calm",
      h_par >= 0 ~ "downwind",
      TRUE ~ "upwind"
    )
  )


###########################################################################################
# Check marginals
###########################################################################################
grid_omsev <- grid_coords_km
adv_matrix <- as.matrix(adv_df_transfo[, c("vx_t", "vy_t")])
Nsim <- 1000
s0_sim <- integer(Nsim)
sims_all <- vector("list", Nsim)
adv_sim <- matrix(0, nrow = Nsim, ncol = 2)
for (i in seq_len(Nsim)) {
  idx <- sample(seq_along(list_episodes), 1)
  s0_i <- s0_list[idx]
  s0_sim[i] <- s0_i
  u_i <- u_list[idx]
  adv_i <- adv_matrix[idx, ]
  adv_sim[i, ] <- adv_i
  sims_all[[i]] <- simulate_many_episodes(
    N = 1,
    u = 1000,
    u_emp = u_i,
    params_vario = params_kmh,
    params_margins = params_margins,
    coords = grid_omsev, # km
    times = times * 5 / 60,  # in hours
    adv = adv_i,
    t0 = 0,
    s0 = s0_i
  )[[1]]
}

Xsim_all <- do.call(cbind, lapply(sims_all, function(sim) sim$X))
sites <- rownames(Xsim_all)
stopifnot(!is.null(sites))

Hgpd <- function(x, sigma, xi) {
  x <- pmax(x, 0)
  if (abs(xi) < 1e-10) return(1 - exp(-x / sigma))
  1 - (1 + xi * x / sigma)^(-1/xi)
}

Qgpd <- function(p, sigma, xi) {
  p <- pmin(pmax(p, 1e-12), 1 - 1e-12)
  if (abs(xi) < 1e-10) return(-sigma * log(1 - p))
  sigma/xi * ((1 - p)^(-xi) - 1)
}

r_egpd_pos <- function(n, kappa, sigma, xi) {
  u <- runif(n)
  Qgpd(u^(1 / kappa), sigma, xi)
}

r_hurdle_egpd <- function(n, p0, kappa, sigma, xi, include_zeros = FALSE) {
  if (!include_zeros) return(r_egpd_pos(n, kappa, sigma, xi))
  z0 <- runif(n) < p0
  x <- numeric(n)
  npos <- sum(!z0)
  if (npos > 0) x[!z0] <- r_egpd_pos(npos, kappa, sigma, xi)
  x
}


library(ggplot2)
library(dplyr)
library(tidyr)

# Split-violin geom (version sans ggname)
GeomSplitViolin <- ggproto(
  "GeomSplitViolin", GeomViolin,
  draw_group = function(self, data, panel_scales, coord, draw_quantiles = NULL) {
    data <- transform(data,
                      xminv = x - violinwidth * (x - xmin),
                      xmaxv = x + violinwidth * (xmax - x))
    grp <- data[1, "group"]
    newdata <- data
    newdata$x <- if (grp %% 2 == 1) data$xminv else data$xmaxv
    newdata <- newdata[order(newdata$y), ]
    newdata <- rbind(newdata[1, ], newdata, newdata[nrow(newdata), ], newdata[1, ])
    newdata[c(1, nrow(newdata)-1, nrow(newdata)), "x"] <- newdata[1, "x"]
    GeomPolygon$draw_panel(newdata, panel_scales, coord)
  }
)

geom_split_violin <- function(mapping = NULL, data = NULL, stat = "ydensity",
                              position = "identity", ..., trim = TRUE, scale = "area",
                              na.rm = FALSE, show.legend = NA, inherit.aes = TRUE) {
  layer(
    data = data, mapping = mapping, stat = stat, geom = GeomSplitViolin,
    position = position, show.legend = show.legend, inherit.aes = inherit.aes,
    params = list(trim = trim, scale = scale, na.rm = na.rm, ...)
  )
}

site <- rownames(sims_all[[1]]$X)[4]
lt <- ncol(sims_all[[1]]$X)

Xsim_site <- do.call(rbind, lapply(sims_all, function(sim) as.numeric(sim$X[site, ])))

Xobs_site <- do.call(rbind, lapply(list_episodes, function(ep) as.numeric(ep[, site])))

stopifnot(ncol(Xsim_site) == lt, ncol(Xobs_site) == lt)

df <- data.frame(
  value = c(Xobs_site, Xsim_site),
  type = rep(c("Observed episodes", "Simulated episodes"))
) %>%
  filter(is.finite(value) & value > 0.2) %>%
  mutate(yplot = value)

ggplot(df, aes(x = type, y = yplot, fill = type)) +
  geom_violin(alpha = 0.7, trim = FALSE, color = NA) +
  labs(
       x = "", y = "Rain") +
  theme_minimal() + btf_theme +
  theme(legend.position = "none")

df_log <- data.frame(
  value = c(Xobs_site, Xsim_site),
  type = rep(c("Observed episodes", "Simulated episodes"))
) %>%
  filter(is.finite(value) & value > 0.1) %>%
  mutate(yplot = log1p(value))

ggplot(df_log, aes(x = type, y = yplot, fill = type)) +
  geom_violin(alpha = 0.7, trim = FALSE, color = NA) +
  labs(x = "", y = "log(1 + rain)") +
  theme_minimal() + btf_theme +
  theme(legend.position = "none")


sites <- rownames(sims_all[[1]]$X)
lt <- ncol(sims_all[[1]]$X)

Xsim_all <- lapply(sims_all, function(sim) sim$X)

sites <- rownames(sims_all[[1]]$X)
sites_plot <- sites_names[13:17]  # choose sites to plot

df_all <- do.call(rbind, lapply(sites_plot, function(site) {

  Xobs_site <- unlist(lapply(list_episodes, function(ep) ep[, site]))
  Xobs_site <- Xobs_site[is.finite(Xobs_site)]

  Xsim_site <- unlist(lapply(sims_all, function(sim) sim$X[site, ]))
  Xsim_site <- Xsim_site[is.finite(Xsim_site)]

  data.frame(
    site  = site,
    value = c(Xobs_site, Xsim_site),
    type  = rep(c("Observed", "Simulated"))
  )
}))


df_all <- df_all %>%
  filter(is.finite(value) & value > 0.22) %>%
  mutate(yplot = value)

# Plot
ggplot(df_all, aes(x = type, y = yplot, fill = type)) +
  geom_violin(alpha = 1, trim = FALSE, color = NA) +
  facet_wrap(~site, scales = "free_y", ncol = 3) +
  labs(
    x = "", y = "Rainfall (mm/5min)"
  ) +
  theme_minimal() + btf_theme +
  theme(legend.position = "none")

ggsave(
  filename = paste0(
    foldername_plot,
    "marginal_distributions_observed_vs_simulated_episodes_",
    "_n", Nsim, "_3",
    ".png"
  ),
  plot = last_plot(),
  width = 10,
  height = 8,
  dpi = 300
)


################################################################################
# On grid OMSEV are
################################################################################



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
  s0_index <- which(rownames(coords) == s0_pixel_id)

  sim <- sim_rpareto_coords(
    beta1 = params_vario$beta1,
    beta2 = params_vario$beta2,
    alpha1 = params_vario$alpha1,
    alpha2 = params_vario$alpha2,
    adv    = adv,
    coords = coords,
    t      = times,
    t0_index     = t0 + 1,
    s0_index     = s0_index,
    threshold = u
  )

  Z <- sim$Z
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


params_margins_common <- list(
  p0    = mean(params_margins$p0),
  xi    = mean(params_margins$xi),
  sigma = mean(params_margins$sigma),
  kappa = mean(params_margins$kappa)
)




library(sf)

# Convertir sites en sf et définir bounding box
sites_sf <- st_as_sf(sites_coords, coords = c("Longitude", "Latitude"), crs = 2154)
bbox <- st_bbox(sites_sf)
bbox_sf <- st_as_sfc(st_bbox(c(xmin = 3.848, ymin = 43.628,
                               xmax = 3.870, ymax = 43.639), crs = 4326))


# Transformer en projection métrique (Lambert 93)
bbox_m <- st_transform(bbox_sf, crs = 2154)
bbox <- st_bbox(bbox_m)

# Créer grille régulière
x_seq <- seq(bbox$xmin, bbox$xmax, by = 100)  # 100 m
y_seq <- seq(bbox$ymin, bbox$ymax, by = 100)

grid_sf <- expand.grid(x = x_seq, y = y_seq) %>%
  st_as_sf(coords = c("x", "y"), crs = 2154)
grid_coords <- st_coordinates(grid_sf)
grid_df <- as.data.frame(grid_coords)
rownames(grid_df) <- paste0("pixel_", seq_len(nrow(grid_df)))
colnames(grid_df) <- c("Longitude", "Latitude")


grid_sf <- st_as_sf(grid_df, coords = c("Longitude", "Latitude"), crs = 2154)
# add pixels id
grid_sf$pixel_id <- rownames(grid_df)
grid_latlon <- st_transform(grid_sf, crs = 4326)
grid_latlon_coords <- st_coordinates(grid_latlon)

grid_df <- as.data.frame(grid_latlon_coords)
rownames(grid_df) <- grid_sf$pixel_id
colnames(grid_df) <- c("Longitude", "Latitude")

# plot on a map pixel and rain gauges
library(ggplot2)
library(ggspatial)  

sites_sf <- st_as_sf(sites_coords, coords = c("Longitude", "Latitude"), crs = 4326)
# add site names
sites_sf$site_name <- rownames(sites_coords)
ggplot() +
  geom_sf(data = bbox_sf, fill = "lightblue", alpha = 0.3) +
  geom_sf(data = grid_latlon, color = "grey80", size = 0.5) +
  geom_sf(data = sites_sf, color = "red", size = 2) +
  geom_sf_text(
    data = sites_sf,
    aes(label = site_name),  # adapte au nom réel de ta colonne
    color = "red",
    size = 3,
    nudge_y = 0.00005
  ) +
  theme_minimal() +
  btf_theme

# Transformer en lat/lon pour ggplot
grid_latlon_poly <- st_transform(grid_l93_poly, crs = 4326)

res <- 100  # 100 m

grid_l93_poly <- st_make_grid(bbox_m, cellsize = res, what = "polygons") |>
  st_sf() |>
  dplyr::mutate(pixel_id = paste0("pixel_", dplyr::row_number()))
grid_l93_pts <- st_centroid(grid_l93_poly)

coords_l93 <- st_coordinates(grid_l93_pts)
coords_df <- as.data.frame(coords_l93)
rownames(coords_df) <- grid_l93_poly$pixel_id
colnames(coords_df) <- c("Longitude", "Latitude")  # (en L93 ici)

grid_latlon_poly <- st_transform(grid_l93_poly, crs = 4326)

# Plot
ggplot() +
  geom_sf(data = grid_latlon_poly, fill = NA, color = "grey80", size = 0.3) +
  geom_sf(data = sites_sf, color = "red", size = 2) +
  theme_minimal() +
  theme(
    axis.title = element_blank(),
    axis.ticks = element_blank(),
    panel.grid = element_blank()
  ) 

# ggplot() +
#   geom_sf(data = grid_latlon_poly, fill = NA, color = "grey80", linewidth = 0.2) +
#   geom_sf_text(
#     data = grid_latlon_poly,
#     aes(label = pixel_id),
#     size = 2,
#     color = "blue"
#   ) +
#   geom_sf(data = sites_sf, color = "red", size = 2) +
#   geom_sf_text(
#     data = sites_sf,
#     aes(label = site_name),  # adapte au nom réel de ta colonne
#     color = "red",
#     size = 3,
#     nudge_y = 0.00005
#   ) +
#   coord_sf(expand = FALSE) +
#   theme_minimal()


# save this plot
ggsave(paste0(im_folder, "swg/omsev/grid_pixels_rain_gauges.png"),
       width = 8, height = 6, dpi = 300)


# choose random conditioning pixel
pixel_ids <- rownames(grid_df)
s0_pixel_id <- sample(pixel_ids, 1)

u_emp <- mean(u_list)
adv <- c(-1, 2)
sim_episode <- sim_episode_grid(
  params_vario = params_kmh,
  params_margins_common = params_margins_common,
  coords = grid_df,
  times = times * 5 / 60,
  adv = adv,
  t0 = 0,
  s0_pixel_id = s0_pixel_id,
  u = 1000,
  u_emp = u_emp
)

df_sim_episode <- data.frame(
  time = seq_len(ncol(sim_episode)),
  rainfall = sim_episode[s0_pixel_id, ]
)

ggplot(df_sim_episode, aes(x = time, y = rainfall)) +
  geom_line(color = btfgreen, size = 1) +
  labs(x = "Time step (5min)", y = "Rainfall (mm/5min)") +
  theme_minimal() +
  btf_theme




grid_latlon_poly <- st_transform(grid_latlon_poly, 4326)
sites_sf <- st_transform(sites_sf, 4326)

stopifnot(all(grid_latlon_poly$pixel_id %in% rownames(sim_episode)))

fill_limits <- range(sim_episode, na.rm = TRUE)

episode_idx <- 1  # adjust if needed
dir_frames <- file.path(im_folder, paste0("swg/omsev/frames_episode_", episode_idx))
dir.create(dir_frames, recursive = TRUE, showWarnings = FALSE)


# Polygones en L93 + centroïdes
grid_l93_poly <- st_transform(grid_latlon_poly, 2154)
grid_l93_cent <- st_centroid(grid_l93_poly)

cent_xy <- st_coordinates(grid_l93_cent)
cent_df <- data.frame(
  pixel_id = grid_l93_poly$pixel_id,
  x = cent_xy[,1],
  y = cent_xy[,2]
)

# Vitesse advection (km/h)
vx <- adv[1]
vy <- adv[2]

# longueur de flèche (visuelle) : déplacement sur dt en heures, puis facteur de zoom
dt_arrow_hours <- 5/60
scale_extra <- 1   # augmente si encore trop court
dx_m <- vx * dt_arrow_hours * 1000 * scale_extra
dy_m <- vy * dt_arrow_hours * 1000 * scale_extra


for (tt in seq_len(nT)) {

  df_t <- data.frame(
    pixel_id = rownames(sim_episode),
    rain = sim_episode[, tt]
  )

  map_t <- grid_latlon_poly %>%
    left_join(df_t, by = "pixel_id")

  # ---- barycentre de pluie (en L93) ----
  tmp <- df_t %>%
    left_join(cent_df, by = "pixel_id") %>%
    mutate(w = pmax(rain, 0))  # ou pmax(rain - p, 0) si tu veux enlever le bruit

  # Si pluie quasi nulle partout, pas de barycentre -> pas de flèche
  if (sum(tmp$w, na.rm = TRUE) > 0) {

    bx <- sum(tmp$x * tmp$w, na.rm = TRUE) / sum(tmp$w, na.rm = TRUE)
    by <- sum(tmp$y * tmp$w, na.rm = TRUE) / sum(tmp$w, na.rm = TRUE)

    arrow_line_l93 <- st_sfc(
      st_linestring(rbind(
        c(bx, by),
        c(bx + dx_m, by + dy_m)
      )),
      crs = 2154
    ) |> st_sf()

    arrow_line_4326 <- st_transform(arrow_line_l93, 4326)

    bary_pt_4326 <- st_transform(st_sfc(st_point(c(bx, by)), crs = 2154) |> st_sf(), 4326)
  } else {
    arrow_line_4326 <- NULL
    bary_pt_4326 <- NULL
  }

  # ---- plot ----
  p <- ggplot() +
    geom_sf(data = map_t, aes(fill = rain), color = NA) +
    geom_sf(data = s0_poly, fill = NA, color = "white", linewidth = 1.2) +
    geom_sf(data = sites_sf, color = "white", size = 3) +
    geom_sf(data = sites_sf, color = "grey", size = 2) +
    scale_fill_viridis_c(
      option = "C",
      limits = fill_limits,
      na.value = "transparent",
      name = "Rainfall\n(mm/5min)"
    ) +
    labs(title = paste("Simulated rainfall field — t =", tt)) +
    theme_minimal() +
    theme(
      axis.title = element_blank(),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      panel.grid = element_blank(),
      plot.title = element_text(size = 12)
    ) +
    coord_sf(expand = FALSE)

  # ---- ajouter barycentre + flèche si dispo et advection non nulle ----
  if (!is.null(arrow_line_4326) && sqrt(vx^2 + vy^2) > 0.1) {

    p <- p +
      # barycentre (petit point)
      geom_sf(data = bary_pt_4326, color = "white", size = 3) +
      geom_sf(data = bary_pt_4326, color = "red", size = 2) +
      # flèche (halo + rouge)
      geom_sf(
        data = arrow_line_4326,
        arrow = arrow(length = unit(0.5, "cm")),
        color = "white",
        linewidth = 2.8
      ) +
      geom_sf(
        data = arrow_line_4326,
        arrow = arrow(length = unit(0.45, "cm")),
        color = "red",
        linewidth = 1.4
      )
  }

  out_png <- file.path(dir_frames, sprintf("frame_%03d.png", tt))
  ggsave(out_png, plot = p, width = 7, height = 6, dpi = 150)
}


library(magick)

files <- list.files(dir_frames, pattern = "frame_\\d+\\.png$", full.names = TRUE)
files <- sort(files)

img <- image_read(files)
gif <- image_animate(img, fps = 1)  # 4 images/sec (ajuste à ton goût)

out_gif <- file.path(im_folder, paste0("swg/omsev/simulated_episode_", episode_idx, ".gif"))
image_write(gif, path = out_gif)
out_gif





