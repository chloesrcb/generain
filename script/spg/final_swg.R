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
# location_gauges$Station <- c("iem", "mse", "poly", "um", "cefe", "cnrs",
#                              "crbm", "archiw", "archie", "um35", "chu1",
#                              "chu2", "chu3", "chu4", "chu5", "chu6", "chu7",
#                              "cines", "brives", "hydro")

# rain <- rain[, !(colnames(rain) %in% c("cines", "hydro", "brives"))]
# location_gauges <- location_gauges[location_gauges$Station != "cines" &
#                                    location_gauges$Station != "hydro" &
#                                    location_gauges$Station != "brives", ]
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
# params_est <-  c(1.101,3.716,0.103,0.676, 3.896, 2.221) #2025
params_est <-  c(1.076,3.911,0.116,0.654, 3.896, 2.221)
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
adv_filename <- paste(data_folder, "/omsev/adv_estim/2025_all/combined_comephore_omsev/episode_advection_q",
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


# get episodes in december 2025
episodes_dec2025 <- selected_points[selected_points$t0_date >= as.POSIXct("2025-12-01 00:00:00", tz = "UTC") &
                                      selected_points$t0_date < as.POSIXct("2026-01-01 00:00:00", tz = "UTC"), ]

# get corresponding advections
sp_dec2025 <- selected_episodes[selected_points$t0_date >= as.POSIXct("2025-12-22 00:00:00", tz = "UTC") &
                                      selected_points$t0_date < as.POSIXct("2025-12-23 00:00:00", tz = "UTC"), ]

# get an episode with advection speed > 10 km/h in selected_episodes
selected_episodes$speed <- sqrt(selected_episodes$adv_x^2 + selected_episodes$adv_y^2)

sp_above10 <- which(selected_episodes$speed > 1)

episode_idx <- 191
sp_adv <- selected_episodes[episode_idx, ]
