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

grid_omsev_km <- grid_omsev

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
tau_vect <- 0:10
tmax <- max(tau_vect)
df_coords <- as.data.frame(grid_omsev)
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
df_lags <- list_lags[[10]]
df_excesses <- list_excesses[[13]]
sum(df_excesses$kij)


#######################################################################################################
# VARIOS PARAMETERS FROM KM/H TO M/5MIN
#######################################################################################################
params_est <- c(1.507 / 2, 4.466 /2 , 0.206, 0.703, 1.092, 5.666)  # km/h
etas_estimates <- params_est[5:6]

params_vario_kmh <- list(
  beta1 = params_est[1],
  beta2 = params_est[2],
  alpha1 = params_est[3],
  alpha2 = params_est[4]
)
params_kmh <- c(
  params_vario_kmh$beta1,
  params_vario_kmh$beta2,
  params_vario_kmh$alpha1,
  params_vario_kmh$alpha2
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
params_m5min <- convert_params(params_vario_kmh$beta1, params_vario_kmh$beta2,
                                params_vario_kmh$alpha1, params_vario_kmh$alpha2,
                               c_x = c_x_m, c_t = c_t_5min)
params_vario_m5min <- list(
  beta1 = params_m5min$beta1,
  beta2 = params_m5min$beta2,
  alpha1 = params_vario_kmh$alpha1,
  alpha2 = params_vario_kmh$alpha2
)

##################################################################################
# TRANSFORM ADVECTION SPEEDS
##################################################################################
adv_filename <- paste(data_folder, "/omsev/adv_estim/combined_comephore_omsev/episode_advection_q",
                          q * 100, "_delta", delta, "_dmin", dmin,
                          ".csv", sep = "")
adv_df <- read.csv(adv_filename, sep = ",")

# count number of advection speed between 0 and 0.1 km/h and with mean_dy_comephore == 0
n_between_0_and_01 <- sum((sqrt(adv_df$vx_final^2 + adv_df$vy_final^2) > 0) &
                           (sqrt(adv_df$vx_final^2 + adv_df$vy_final^2) < 0.1))
cat("Number of episodes with advection between 0 and 0.1 km/h:", n_between_0_and_01, "\n")

index_01_com0 <- which((sqrt(adv_df$vx_final^2 + adv_df$vy_final^2) > 0) &
                      (sqrt(adv_df$vx_final^2 + adv_df$vy_final^2) < 0.1) &
                      (adv_df$mean_dy_kmh_comephore == 0))
adv_df[index_01_com0, c("vx_final", "vy_final")] <- 0

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
     xlab = "vx", ylab = "vy", xlim = c(-10,10), ylim = c(-10,10))

# Compute speed
speed_raw <- sqrt(adv_df$vx_final^2 + adv_df$vy_final^2)
summary(speed_raw)
speed_transfo <- sqrt(adv_df_transfo$vx_t^2 + adv_df_transfo$vy_t^2)
summary(speed_transfo)

direction <- atan2(adv_df$vy_final, adv_df$vx_final) * (180 / pi)
# make sure direction is in [0, 360]
angle_deg <- (direction + 360) %% 360
is_calm <- speed_raw == 0
speed_class <-as.character(cut(
    speed_raw,
    breaks = c(0, 0.5, 1, Inf),
    labels = c("still","weak","significant"),
    include.lowest = TRUE
  ))

speed_class <- factor(speed_class, levels = c("still","weak","significant"))
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
  speed_class = speed_class,
  direction = direction,
  direction_class = direction_class,
  group = paste0(speed_class, "_", direction_class)
)

table(adv_class$group)
adv_df <- adv_df_transfo # km/h 

selected_points$adv_group <- adv_class$group

adv_df_m5min <- adv_df_transfo
adv_df_m5min$vx_final <- adv_df_transfo$vx_final * 1000 / 12  # m/5min
adv_df_m5min$vy_final <- adv_df_transfo$vy_final * 1000 / 12  # m/5min
adv_df_m5min$vx_t <- adv_df_transfo$vx_t * 1000 / 12  # m/5min
adv_df_m5min$vy_t <- adv_df_transfo$vy_t * 1000 / 12  # m/5min


head(adv_class)
# plot of histogram of advection speeds
ggplot(adv_class, aes(x = direction_class)) +
  geom_bar(aes(fill = speed_class), position = "dodge", color = "black") +
  scale_fill_brewer(palette = "YlGnBu", name = "Speed class") +
  labs(title = "", x = "Direction class", y = "Episodes count") +
  theme_minimal() +
  btf_theme


# Plot advection classes
plot_dir_df <- adv_class
plot_dir_df$speed_class <- speed_class
# remove none direction class
plot_dir_df <- plot_dir_df[plot_dir_df$group != "still_none", ]
plot_dir_df$direction_class <- factor(
  plot_dir_df$direction_class,
  levels = c("N", "E", "S", "W")
)

# none direction number
none_count <- sum(adv_class$direction_class == "none")
cat("Number of episodes with 'none' direction class:", none_count, "\n")

ggplot(plot_dir_df, aes(x = direction_class, fill = speed_class)) +
  geom_bar(stat = "count", width = 1, color = "black") +
  coord_polar(start = -45 * (pi/180)) +
  scale_fill_brewer(palette = "YlGnBu", name = "Speed class") +
  labs(title = paste0("Number of episodes without direction: ", none_count), x = "Direction class", y = "Episodes count (advective only)") +
  theme_minimal() +
  btf_theme

min_spatial_dist <- dmin
file_name <- paste0(
  im_folder, "/swg/omsev/adv_classes_",
  q * 100, "q",
  min_spatial_dist[1], "dmin",
  delta[1], "delta_NSEW.png"
)

ggsave(filename = file_name)


ggplot(plot_dir_df, aes(x = speed_class)) +
  geom_bar(stat = "count", width = 0.5, fill = btfgreen, alpha = 0.6) +
  scale_x_discrete(drop = TRUE) +
  labs(title = "", x = "Speed class (km/h)", y = "Episodes count") +
  theme_minimal() +
  btf_theme

ggsave(filename = paste0(im_folder, "/optim/comephore/adv_classes_",
  q*100, "q", min_spatial_dist, "dmin", delta, "delta_speedonly.png"))


################################################################################
# initial excess by episode
################################################################################

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
adv_matrix <- as.matrix(adv_df_transfo[, c("vx_t", "vy_t")])
all_group_names <- unique(selected_points$adv_group)
group_adv <- all_group_names[12]  # choose one group to simulate
#get list_episodes and s0_list for this group
indices_group <- which(selected_points$adv_group == group_adv)
list_episodes_group <- list_episodes[indices_group]
length(list_episodes_group)
s0_list_group <- s0_list[indices_group]
u_list_group <- u_list[indices_group]
adv_group <- adv_matrix[indices_group, ]
Nsim <- 100
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
    params_vario = params_vario_kmh,
    params_margins = params_margins,
    coords = grid_omsev,
    times = times,
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


# save plot

foldername_plot <- paste0(
  im_folder,
  "swg/omsev/all/"
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

################################################################################
# Validation chi
################################################################################

dist_mat <- get_dist_mat(location_gauges) / 1000
df_dist <- reshape_distances(dist_mat)
n_hbins <- 10
h_all <- df_dist$value

h_breaks <- quantile(
  h_all,
  probs = seq(0, 1, length.out = n_hbins + 1),
  na.rm = TRUE
)

h_breaks <- unique(as.numeric(h_breaks))
h_breaks[length(h_breaks)] <- 1.6


# on simulated episodes
list_results_sim <- lapply(seq_along(sims_group), function(i) {

  Z_ep <- sims_group[[i]]$Z
  s0 <- s0_sim[i]
  row_s0 <- which(rownames(df_coords) == s0)
  s0_coords <- df_coords[row_s0, ]
  ind_t0_ep <- 0

  u_0 <- sims_group[[i]]$u_latent

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

  lags$tau <- lags$tau * 5 / 60  # convert to hours
  if(sum(excesses$kij) == 0) {
    cat("No exceedances for simulated episode", i, "\n")
  }
  list(lags = lags, excesses = excesses)
})



list_lags_sim <- lapply(list_results_sim, `[[`, "lags")
list_excesses_sim <- lapply(list_results_sim, `[[`, "excesses")


list_sumkij <- sapply(list_excesses_sim, function(df) sum(df$kij))
summary(list_sumkij)
table(list_sumkij > 0)


adv_df_sim <- as.data.frame(adv_sim)
colnames(adv_df_sim) <- c("vx", "vy")

params_estimates <- as.numeric(c(params_vario_kmh, etas_estimates))

p_sim <- plot_th_emp_chi(
  list_lags_sim,
  list_excesses_sim,
  wind_df_filtered = adv_df_sim,
  params_estimates = params_estimates,
  tau_min = 0,
  h_breaks = h_breaks,
  latlon = FALSE,
  adv_transform = FALSE
)

print(p_sim$plots$all)

# save plot
foldername_plot <- paste0(
  im_folder,
  "swg/omsev/simu"
)

if (!dir.exists(foldername_plot)) {
  dir.create(foldername_plot, recursive = TRUE)
}
ggsave(
  filename = paste0(
    foldername_plot,
    "chi_vs_h_bytau_simulated_advgroup_",
    group_adv,
    "_n", Nsim,
    ".png"
  ),
  plot = last_plot(),
  width = 7,
  height = 7,
  dpi = 300
)

# plot by tau
res_cmp_sim <- p_sim$res_cmp

library(ggplot2)

ggplot(res_cmp_sim, aes(x = chi_theo_bar, y = chi_emp_bar, color = factor(tau))) +
  geom_point(alpha = 0.7, size = 2) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "black") +
  scale_color_viridis_d(name = expression(tau)) +
  labs(
    x = expression(chi[theo]),
    y = expression(chi[emp]),
    title = expression("Empirical vs theoretical r-extremogram"),
    subtitle = "Colour = temporal lag τ"
  ) +
  theme_minimal()


# get h_mid from hbin of form  [0,0.228]
res_cmp_sim$h_mid <- sapply(res_cmp_sim$hbin, function(hb) {
  hb_clean <- gsub("\\[|\\]|\\(|\\)", "", hb)
  limits <- as.numeric(strsplit(hb_clean, ",")[[1]])
  mean(limits)
})
