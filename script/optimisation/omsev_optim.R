rm(list = ls())

# clear console
cat("\014")

library(data.table)
library(latex2exp)
library(lubridate)
library(fuzzyjoin)
library(grid)

muse <- FALSE

if (muse) {
  # Get the muse folder
  folder_muse <- "/home/serrec/work_rainstsimu/rpareto/"
  setwd(folder_muse)
  # Load libraries and set theme
  source("load_libraries.R")
  im_folder <- "./images"
  source("config_omsev.R")
  data_folder <- "./data/"
  ncores <- 27
} else {
  # Load libraries and set theme
  source("./script/load_libraries.R")
  source("./script/optimisation/config_omsev.R")
  ncores <- detectCores() - 1
}

# Get all files in the folder "R"
functions_folder <- "./R"
files <- list.files(functions_folder, full.names = TRUE)
# # load all functions in files
invisible(lapply(files, function(f) source(f, echo = FALSE)))

# get rain data from omsev
# filename_omsev <- paste0(data_folder,
#                          "omsev/omsev_5min/rain_mtp_5min_2019_2024_cleaned.csv")
# # # write.csv(rain, filename_omsev, row.names = FALSE)
# rain <- read.csv(filename_omsev)
# head(rain)
# tail(rain)

# rain_omsev <- rain
# colnames(rain_omsev)
filename_loc <- paste0(data_folder,
                           "omsev/loc_rain_gauges.csv")
# get location of each rain gauge
location_gauges <- read.csv(filename_loc)
location_gauges$Station <- c("iem", "mse", "poly", "um", "cefe", "cnrs",
                             "crbm", "archiw", "archie", "um35", "chu1",
                             "chu2", "chu3", "chu4", "chu5", "chu6", "chu7",
                             "cines", "brives", "hydro")

# get rain data
filename_rain <- paste0(data_folder, "omsev/omsev_5min/rain_mtp_5min_2019_2025.csv")
rain_omsev <- read.csv(filename_rain)
head(rain_omsev)
tail(rain_omsev)
# put dates as rownames
rownames(rain_omsev) <- rain_omsev$dates
# remove dates after january 2025
# rain_omsev <- rain_omsev[rownames(rain_omsev) < "2025-01-31 00:00:00", ]
rain <- rain_omsev[-1] # remove dates column
# rain$mse
# Get distances matrix
dist_mat <- get_dist_mat(location_gauges)
df_dist <- reshape_distances(dist_mat)

# remove cines, hydro, brives
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
grid_coords_m <- as.data.frame(coords_m)
colnames(grid_coords_m) <- c("Longitude", "Latitude")
rownames(grid_coords_m) <- rownames(sites_coords)
grid_coords_km <- as.data.frame(coords_m / 1000)
colnames(grid_coords_km) <- c("Longitude", "Latitude")
rownames(grid_coords_km) <- rownames(sites_coords)


################################################################################
# VARIOGRAM --------------------------------------------------------------------
################################################################################

# get first date with at least 3 active sites
# rain$active_sites <- rowSums(rain > 0, na.rm = TRUE)
# rain_filtered <- rain[rain$active_sites >= 3, ]
# head(rain_filtered)
# first_date <- rownames(rain_filtered)[1]
# rain <- rain[rownames(rain) >= first_date, ]
# remove active_sites column
# rain <- rain[, -ncol(rain)]
# rain <- rain_complete
# in rain remove when all data are NA
set_st_excess <- get_spatiotemp_excess(rain, quantile = q, remove_zeros = TRUE)

# verify that the excess is above the threshold
# get list of sites and times
list_s <- set_st_excess$list_s
list_t <- set_st_excess$list_t
list_u <- set_st_excess$list_u

min_spatial_dist <- 1200
delta <- 12
# Spatio-temporal neighborhood parameters
s0t0_set <- get_s0t0_pairs(grid_coords_m, rain,
                            min_spatial_dist = min_spatial_dist,
                            episode_size = delta,
                            set_st_excess = set_st_excess,
                            n_max_episodes = 10000,
                            latlon = FALSE)

selected_points <- s0t0_set
length(selected_points$s0)
selected_points$t0
same_temporality <- function(t0_i, t0_j, delta, step_min = 5) {
  start_i <- t0_i
  end_i   <- t0_i + (delta - 1) * step_min * 60
  start_j <- t0_j
  end_j   <- t0_j + (delta - 1) * step_min * 60
  
  !(end_i < start_j || end_j < start_i)
}


t0_vec <- selected_points$t0_date
t0_vec <- as.POSIXct(t0_vec, format = "%Y-%m-%d %H:%M:%S", tz = "UTC")
n_ep <- length(t0_vec)

temporal_overlap <- matrix(FALSE, n_ep, n_ep)

for (i in 1:n_ep) {
  for (j in 1:n_ep) {
    temporal_overlap[i, j] <- same_temporality(t0_vec[i], t0_vec[j], delta)
  }
}

has_temporal_neighbor <- apply(temporal_overlap, 1, function(x) sum(x) > 1)
table(has_temporal_neighbor)


# dmin <- c(500, 600, 700, 800, 900, 1000, 1100, 1200, 1300, 1400, 1500, 1600)
# n_episodes <- sapply(dmin, function(d) {
#   s0t0_set <- get_s0t0_pairs(grid_coords_m, rain,
#                               min_spatial_dist = d,
#                               episode_size = delta,
#                               set_st_excess = set_st_excess,
#                               n_max_episodes = 10000,
#                               latlon = FALSE)
#   length(s0t0_set$s0)
# })

# # plot n_episodes by dmin
# df_n_episodes <- data.frame(dmin = dmin, n_episodes = n_episodes)
# ggplot(df_n_episodes, aes(x = dmin, y = n_episodes)) +
#   geom_line(color=btfgreen) +
#   geom_point(color=btfgreen) +
#   ylim(0, 800) +
#   geom_vline(xintercept = 1000, linetype = "dashed", color = "red") +
#   labs(x = "Minimum spatial distance (m)", y = "Number of episodes") +
#   btf_theme

# filename <- paste0(im_folder, "/optim/omsev/choice_config/n_episodes_dmin1000_delta", delta, ".png")

# if (!dir.exists(paste0(im_folder, "/optim/omsev/choice_config/"))) {
#   dir.create(paste0(im_folder, "/optim/omsev/choice_config/"), recursive = TRUE)
# }

# ggsave(filename, width = 7, height = 5)



# deltas <- seq(2, 24, by = 2)  # in steps of 5 minutes
# n_episodes <- sapply(deltas, function(d) {
#    s0t0_set <- get_s0t0_pairs(grid_coords_m, rain,
#                               min_spatial_dist = min_spatial_dist,
#                               episode_size = d,
#                               set_st_excess = set_st_excess,
#                               n_max_episodes = 10000,
#                               latlon = FALSE)
#    length(s0t0_set$s0)
# })

# # plot n_episodes by dmin
# df_n_episodes <- data.frame(delta = deltas * 5, n_episodes = n_episodes)
# ggplot(df_n_episodes, aes(x = delta, y = n_episodes)) +
#   geom_line(color=btfgreen) +
#   geom_point(color=btfgreen) +
#   ylim(0, 500) +
#   geom_vline(xintercept = delta * 5, linetype = "dashed", color = "red") +
#   labs(x = "Episode size (minutes)", y = "Number of episodes") +
#   btf_theme

# filename <- paste0(im_folder, "/optim/omsev/choice_config/n_episodes_delta_dmin", min_spatial_dist, ".png")

# if (!dir.exists(paste0(im_folder, "/optim/omsev/choice_config/"))) {
#   dir.create(paste0(im_folder, "/optim/omsev/choice_config/"), recursive = TRUE)
# }

# ggsave(filename, width = 7, height = 5)

selected_points <- selected_points %>%
  mutate(t0_date = as.POSIXct(t0_date, format = "%Y-%m-%d %H:%M:%S", tz = "UTC"))
n_episodes <- nrow(selected_points)
# For each (s0, t0, u_s0), check if rain[t0, s0] > u_s0
excess_check_s0t0 <- s0t0_set[, {
  rain_val <- rain[t0, s0]
  is_excess <- rain_val > u_s0
  list(rain_value = rain_val, is_excess = is_excess)
}, by = .(s0, t0, u_s0)]


# Threshold histogram
df_threshold <- data.frame(u_s0 = selected_points$u_s0)
breaks <- seq(floor(min(df_threshold$u_s0)), ceiling(max(df_threshold$u_s0)), by = 0.1)

n_episodes <- length(selected_points$s0)
t0_list <- selected_points$t0
s0_list <- selected_points$s0

list_episodes_points <- get_extreme_episodes(selected_points, rain,
                                     episode_size = delta, unif = FALSE)
list_episodes <- list_episodes_points$episodes

selected_points$t0_date <- as.POSIXct(selected_points$t0_date, format = "%Y-%m-%d %H:%M:%S", tz = "UTC")


# Round to the next hour
selected_points$t0_date_rounded <- ceiling_date(selected_points$t0_date, "hour")

datetimes <- unique(selected_points$t0_date)


datetimes_hour <- unique(selected_points$t0_date_rounded)


# save datetime list to csv
datetime_filename <- paste(data_folder, "/omsev/t0_episodes_q", q * 100,
                           "_delta", delta, "_dmin", min_spatial_dist,
                           ".csv", sep = "")
write.csv(data.frame(t0_date = datetimes_hour), datetime_filename, row.names = FALSE)


# save datetime list to csv
datetime_filename <- paste(data_folder, "/omsev/t0_5min_episodes_q", q * 100,
                           "_delta", delta, "_dmin", min_spatial_dist,
                           ".csv", sep = "")
write.csv(data.frame(t0_date = datetimes), datetime_filename, row.names = FALSE)


adv_filename <- paste(data_folder, "/omsev/adv_estim/2025/combined_comephore_omsev/episode_advection_q",
                          q * 100, "_delta", delta, "_dmin", min_spatial_dist,
                          ".csv", sep = "")
adv_df <- read.csv(adv_filename, sep = ",")

head(adv_df)
tail(adv_df)
# get index of advection speed below 0.1 km/h
adv_speed <- sqrt(adv_df$vx_final^2 + adv_df$vy_final^2)
summary(adv_speed)
idx_low_adv <- which(adv_speed > 0 & adv_speed < 0.1)

# get index of advection speed below 0.1 km/h and with mean_dx_kmh_comphore and mean_dy_kmh_comephore at 0
idx_low_speed_adv_comephore <- which(adv_speed > 0 & adv_speed < 0.1 & adv_df$mean_dx_kmh_comephore == 0 & adv_df$mean_dy_kmh_comephore == 0)
print(paste("Number of episodes with advection speed below 0.1 km/h and comephore advection at 0:", length(idx_low_speed_adv_comephore)))
total_n_episodes <- nrow(adv_df)
print(paste("Total number of episodes:", total_n_episodes))
# put these advections to 0
# adv_df$vx_final[idx_low_speed_adv_comephore] <- 0
# adv_df$vy_final[idx_low_speed_adv_comephore] <- 0

# get only matching episodes from selected_points
# convert adv_df$t0 to POSIXct
adv_df$t0_omsev <- as.POSIXct(adv_df$t0_omsev, format="%Y-%m-%d %H:%M:%S", tz="UTC")
selected_points$t0_date <- as.POSIXct(selected_points$t0_date, format="%Y-%m-%d %H:%M:%S", tz="UTC")
matching_indices <- sapply(selected_points$t0_date, function(t) {
  diffs <- abs(difftime(adv_df$t0_omsev, t, units = "secs"))
  idx <- which.min(diffs)
  if (diffs[idx] > 60*5) {  # tolérance = 5 minutes
    return(NA)
  } else {
    return(idx)
  }
})

# remove NA indices
matching_indices <- matching_indices[!is.na(matching_indices)]
# get only matching rows
adv_df <- adv_df[matching_indices,]
tail(adv_df)
rownames(adv_df) <- NULL  # reset row names

selected_episodes <- selected_points
selected_episodes$adv_x <- rep(NA, nrow(selected_episodes))
selected_episodes$adv_y <- rep(NA, nrow(selected_episodes))

# get adv values for each episode according to the t0_date
setDT(selected_points)
setDT(adv_df)

adv_df[, t0_omsev := as.POSIXct(t0_omsev, format="%Y-%m-%d %H:%M:%S", tz="UTC")]
selected_points[, t0_date := as.POSIXct(t0_date, format="%Y-%m-%d %H:%M:%S", tz="UTC")]

adv_df_u <- adv_df[, .(
  vx_final = vx_final[1],
  vy_final = vy_final[1]
), by = t0_omsev]

setkey(adv_df_u, t0_omsev)

selected_episodes <- adv_df_u[selected_points,
  on   = .(t0_omsev = t0_date)
]

setnames(selected_episodes, c("vx_final","vy_final"), c("adv_x","adv_y"))


V_episodes <- data.frame(
  v_x = selected_episodes$adv_x,
  v_y = selected_episodes$adv_y
)

colnames(V_episodes) <- c("vx", "vy")


# remove NA
V_episodes <- V_episodes[!is.na(V_episodes$vx) & !is.na(V_episodes$vy), ]

# max of advection
max_adv <- max(sqrt(V_episodes$vx^2 + V_episodes$vy^2))
# number of 0 adv
n_zero_adv <- sum(V_episodes$vx == 0 & V_episodes$vy == 0)
cat("Number of episodes with zero advection:", n_zero_adv, "\n")

n_between_0_and_01 <- sum((sqrt(V_episodes$vx^2 + V_episodes$vy^2) > 0) &
                           (sqrt(V_episodes$vx^2 + V_episodes$vy^2) < 0.1))
cat("Number of episodes with advection between 0 and 0.1 km/h:", n_between_0_and_01, "\n")

tau_vect <- 0:10

selected_episodes_nona <- selected_episodes[!is.na(selected_episodes$adv_x) & !is.na(selected_episodes$adv_y), ]
s0_list <- selected_episodes_nona$s0
list_episodes_points <- get_extreme_episodes(selected_episodes_nona, rain,
                              episode_size = delta, unif = FALSE)
list_episodes <- list_episodes_points$episodes
u0_list <- selected_episodes_nona$u_s0
tmax <- max(tau_vect)
df_coords <- as.data.frame(grid_coords_km)
# Compute the lags and excesses for each conditional point
list_results <- mclapply(1:length(s0_list), function(i) {
  s0 <- s0_list[i]
  row_s0 <- which(rownames(df_coords) == s0)
  s0_coords <- df_coords[row_s0, ]
  episode <- list_episodes[[i]]
  ind_t0_ep <- 0 # index of t0 in the episode
  lags <- get_conditional_lag_vectors(df_coords, s0_coords, ind_t0_ep,
                                tau_vect, latlon = FALSE)
  u <- u0_list[i]
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
df_lags$hnorm
df_excesses <- list_excesses[[15]]
sum(df_excesses$kij)


# -------------------------------------------------------------------------
# Detect temporal overlaps and restrict spatial domain only for those episodes
# -------------------------------------------------------------------------

selected_episodes$t0_date <- as.POSIXct(
  selected_episodes$t0_date,
  format = "%Y-%m-%d %H:%M:%S",
  tz = "UTC"
)

# same_temporality <- function(t0_i, t0_j, delta, step_min = 5) {
#   abs(as.numeric(difftime(t0_i, t0_j, units = "mins"))) < delta * step_min
# }

# t0_vec <- selected_episodes$t0_date
# n_ep <- length(t0_vec)

# temporal_overlap <- matrix(FALSE, n_ep, n_ep)

# for (i in 1:n_ep) {
#   for (j in 1:n_ep) {
#     if (i != j) {
#       temporal_overlap[i, j] <- same_temporality(t0_vec[i], t0_vec[j], delta)
#     }
#   }
# }

# has_temporal_overlap <- apply(temporal_overlap, 1, any)
# table(has_temporal_overlap)

# choose local radius for episodes with temporal overlap
# IMPORTANT: adapt depending on unit of hnorm
# d_local <- 0.8   # if hnorm is in km
# d_local <- 600 # if hnorm is in meters

# list_lags_overlap <- vector("list", length(list_lags))
# list_excesses_overlap <- vector("list", length(list_excesses))

# for (i in seq_along(list_lags)) {
#   if (has_temporal_overlap[i]) {
#     keep <- list_lags[[i]]$hnorm <= d_local
#   } else {
#     keep <- rep(TRUE, nrow(list_lags[[i]]))
#   }

#   list_lags_overlap[[i]] <- list_lags[[i]][keep, , drop = FALSE]
#   list_excesses_overlap[[i]] <- list_excesses[[i]][keep, , drop = FALSE]
# }

# selected_episodes_filtered <- selected_episodes
# list_episodes_filtered <- list_episodes
# list_lags_filtered <- list_lags
# list_excesses_filtered <- list_excesses

# get comephore estimates
# filename_com_res <- paste(data_folder,
#         "comephore/optim_results/euclidean/free_eta/results_com_classes_q95_delta30_dmin5.csv",
#         sep = "")
# com_results <- read.csv(filename_com_res)
com_results <- c(0.1579947, 0.9923270, 0.5800128, 0.6627969, 3.8962660, 2.2208320)

# check for na in adv and wind
hmax <- max(dist_mat) / 1000


# Basic diagnostics on advection (optional)
selected_episodes$adv_speed <- sqrt(selected_episodes$adv_x^2 +
                                    selected_episodes$adv_y^2)
selected_episodes$adv_direction <- atan2(selected_episodes$adv_y,
                                         selected_episodes$adv_x) * (180 / pi)
selected_episodes$adv_direction[selected_episodes$adv_direction < 0] <-
  selected_episodes$adv_direction[selected_episodes$adv_direction < 0] + 360
hist(selected_episodes$adv_speed, main = "Histogram of advection speed",
     xlab = "Speed (km/h)", breaks = 30)

# get advection with speed above 0.2 km/h
speed_threshold <- 0
adv_speed <- sqrt(selected_episodes$adv_x^2 +
                    selected_episodes$adv_y^2)

# max(adv_speed)
# # put all adv < speed_threshold to 0
# adv_speed[adv_speed < speed_threshold] <- 0
# selected_episodes$adv_speed <- adv_speed
# selected_episodes$adv_x[adv_speed == 0] <- 0
# selected_episodes$adv_y[adv_speed == 0] <- 0
# filtered <- FALSE
# if (filtered) {
#   selected_episodes_filtered <- selected_episodes[adv_speed >= speed_threshold, ]
#     list_episodes_filtered <- list_episodes[adv_speed > speed_threshold]
#     list_lags_filtered <- list_lags[adv_speed > speed_threshold]
#     list_excesses_filtered <- list_excesses[adv_speed > speed_threshold]
# } else {
  selected_episodes_filtered <- selected_episodes
    list_episodes_filtered <- list_episodes
    list_lags_filtered <- list_lags
    list_excesses_filtered <- list_excesses

# }

V_episodes_filtered <- data.frame(
    vx = selected_episodes_filtered$adv_x,
    vy = selected_episodes_filtered$adv_y
)
number_episode <- length(list_episodes_filtered)

print(paste("Number of episodes after filtering:", number_episode))

# look in rain omsev for which dates there are at least 3 sites activated
# rain_omsev$active_sites <- rowSums(rain_omsev[, -1] > 0, na.rm = TRUE)
# rain_omsev_filtered <- rain_omsev[rain_omsev$active_sites >= 3, ]
# head(rain_omsev_filtered)
# head(rain_omsev)
# remove episodes before 2019-10-01
# id_filter_date <- which(selected_episodes_filtered$t0_date < as.POSIXct("2019-11-01", tz = "UTC"))
# selected_episodes_filtered <- selected_episodes_filtered[-id_filter_date, ]
# list_episodes_filtered <- list_episodes_filtered[-id_filter_date]
# list_lags_filtered <- list_lags_filtered[-id_filter_date]
# list_excesses_filtered <- list_excesses_filtered[-id_filter_date]
# V_episodes_filtered <- V_episodes_filtered[-id_filter_date, ]
# length(list_episodes_filtered)
# OPTIMISATION ---------------------------------------------------------------
# [1] 1.1533284 4.1566903 0.2059604 0.6466280 avec les 414 episodes 

init_com_class <- com_results
eta1_class <- init_com_class[5]
eta2_class <- init_com_class[6]
head(V_episodes_filtered)
result <- optim(
  par = init_com_class[1:4],
  fn = neg_ll_composite_fixed_eta,
  list_lags = list_lags_filtered,
  list_episodes = list_episodes_filtered,
  list_excesses = list_excesses_filtered,
  hmax = hmax,
  wind_df = V_episodes_filtered,
  latlon = FALSE,
  distance = "euclidean",
  fixed_eta1 = eta1_class,
  fixed_eta2 = eta2_class,
  method = "L-BFGS-B",
  lower = c(1e-08, 1e-08, 1e-08, 1e-08),
  upper = c(10, 10, 1.999, 1.999),
  control = list(maxit = 20000, trace = 1)
)

if(result$convergence == 0) {
  cat("Optimization converged successfully\n")
} else {
  warning("Optimization did not converge")
}
# [1] 1.0603714 4.4023551 0.1923471 0.6743526 5.3556250 2.1357170 dmin=1000
result

convert_params <- function(beta1, beta2, alpha1, alpha2, c_x = 1, c_t = 1) {
  beta1_new <- beta1 / (c_x^alpha1)
  beta2_new <- beta2 / (c_t^alpha2)
  list(beta1 = beta1_new, beta2 = beta2_new)
}

par_m5min <- convert_params(result$par[1], result$par[2],
                               result$par[3], result$par[4],
                               c_x = 1000, c_t = 12)

params_estim_m5min <- c(par_m5min$beta1,
                        par_m5min$beta2,
                        result$par[3],
                        result$par[4],
                        eta1_class,
                        eta2_class)

param_estim_kmh <- c(result$par[1],
                     result$par[2],
                     result$par[3],
                     result$par[4],
                     eta1_class,
                     eta2_class)


if (result$convergence != 0) {
  warning("Optimization did not converge")
} else {
  cat("Optimization converged successfully\n")
}
result_df <- data.frame(beta1 = result$par[1],
                        beta2 = result$par[2],
                        alpha1 = result$par[3],
                        alpha2 = result$par[4])

eta_type <- "fixed_eta"
foldername_res <- file.path(
  paste0(data_folder, "omsev/optim_results"),
  distance_type,
  eta_type
)

if (!dir.exists(foldername_res)) {
  dir.create(foldername_res, recursive = TRUE)
}
# filename <- paste0(foldername_res, "/results_q",
#            q * 100, "_delta", delta, "_dmin", min_spatial_dist, ".csv")
# # write.csv(result_df, filename, row.names = FALSE)

# result_df <- read.csv(filename)

params_est <- as.numeric(result_df)

params_est <- c(params_est, eta1_class, eta2_class)

# Jackknife CI ---------------------------------------------------------------

# Initial parameters for jackknife (from full data optimization)
init_param_jk <- as.numeric(result_df)

#### MONTH-YEAR JACKKNIFE ------------------------------------------------------

selected_episodes_filtered$t0_date <- as.POSIXct(
  selected_episodes_filtered$t0_date,
  format = "%Y-%m-%d %H:%M:%S",
  tz = "UTC"
)
ncores <- detectCores() - 1
month_vec <- format(selected_episodes_filtered$t0_date, "%m")
year_vec  <- format(selected_episodes_filtered$t0_date, "%Y")
month_year_vec <- paste(year_vec, month_vec, sep = "-")

unique_month_years <- sort(unique(month_year_vec))
n_month_years <- length(unique_month_years)

jack_estimates_list <- parallel::mclapply(unique_month_years, function(month_year) {
  cat("Excluding month-year:", month_year, "\n")
  exclude_idx <- which(month_year_vec == month_year)
  
  jack_episodes <- list_episodes_filtered[-exclude_idx]
  jack_lags <- list_lags_filtered[-exclude_idx]
  jack_excesses <- list_excesses_filtered[-exclude_idx]
  jack_wind <- V_episodes_filtered[-exclude_idx, , drop = FALSE]
  
  res <- tryCatch({
    optim(par = init_param_jk[1:4], fn = neg_ll_composite_fixed_eta,
      list_lags = jack_lags, list_episodes = jack_episodes,
      list_excesses = jack_excesses, hmax = hmax,
      wind_df = jack_wind,
      latlon = FALSE,
      distance = "euclidean",
      fixed_eta1 = eta1_class,
      fixed_eta2 = eta2_class,
      method = "L-BFGS-B",
      lower = c(1e-08, 1e-08, 1e-08, 1e-08),
      upper = c(10, 10, 1.999, 1.999),
      control = list(maxit = 10000),
      hessian = FALSE)
  }, error = function(e) NULL)
  
  if (!is.null(res) && res$convergence == 0) {
    return(res$par)
  } else {
    print("Optimization failed or did not converge")
    return(rep(NA, length(init_param_jk)))
  }
}, mc.cores = ncores)


jack_estimates <- do.call(rbind, jack_estimates_list)
jack_estimates <- na.omit(jack_estimates)
n_eff <- nrow(jack_estimates)

foldername_jk <- paste0(data_folder, "omsev/optim_results/jackknife_estimates/")
if (!dir.exists(foldername_jk)) {
  dir.create(foldername_jk, recursive = TRUE)
}
filename <- paste0(foldername_jk, "all_results_jk_by_monthyear_n",
           n_eff, "_q", q*100, "_delta", delta, "_dmin", min_spatial_dist, ".csv")
write.csv(jack_estimates, filename, row.names = FALSE)


jack_estimates <- do.call(rbind, jack_estimates_list)
jack_estimates <- na.omit(jack_estimates)
n_eff <- nrow(jack_estimates)

log_jack_estimates <- log(jack_estimates)
log_init_param_jk  <- log(init_param_jk[1:4])

pseudo_values_log <- matrix(NA, nrow = n_eff, ncol = 4)
for (i in 1:n_eff) {
  pseudo_values_log[i, ] <- n_eff * log_init_param_jk - (n_eff - 1) * log_jack_estimates[i, ]
}

jack_mean_pseudo_log <- colMeans(pseudo_values_log)
jack_se_log          <- apply(pseudo_values_log, 2, sd) / sqrt(n_eff)

z <- qnorm(0.975)
lower_ci_log <- jack_mean_pseudo_log - z * jack_se_log
upper_ci_log <- jack_mean_pseudo_log + z * jack_se_log

estimate_jk_naturel <- exp(jack_mean_pseudo_log)
lower_ci <- exp(lower_ci_log)
upper_ci <- exp(upper_ci_log)

jackknife_monthyear_results <- data.frame(
  Parameter     = c("beta1", "beta2", "alpha1", "alpha2"),
  Estimate_full = result$par,
  Estimate_jk   = estimate_jk_naturel,
  StdError_log  = jack_se_log,
  CI_lower      = lower_ci,
  CI_upper      = upper_ci
)


filename <- paste0(foldername_jk, "results_jk_by_monthyear_n",
           n_eff, "_q", q*100, "_delta", delta, "_dmin", min_spatial_dist, ".csv")
write.csv(jackknife_monthyear_results, filename, row.names = FALSE)
jackknife_monthyear_results <- read.csv(filename)
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

beta1_m5min_i <- jack_estimates[,1] / (1000 ^ jack_estimates[,3])
beta2_m5min_i <- jack_estimates[,2] / (12   ^ jack_estimates[,4])

b1L <- jackknife_monthyear_results$CI_lower[1]
b1U <- jackknife_monthyear_results$CI_upper[1]
a1L <- jackknife_monthyear_results$CI_lower[3]
a1U <- jackknife_monthyear_results$CI_upper[3]

b2L <- jackknife_monthyear_results$CI_lower[2]
b2U <- jackknife_monthyear_results$CI_upper[2]
a2L <- jackknife_monthyear_results$CI_lower[4]
a2U <- jackknife_monthyear_results$CI_upper[4]

# conversion correcte
beta1_m5_low  <- b1L / (1000 ^ a1U)
beta1_m5_high <- b1U / (1000 ^ a1L)
beta2_m5_low  <- b2L / (12   ^ a2U)
beta2_m5_high <- b2U / (12   ^ a2L)


# Pseudo-valeurs sur l’échelle m/5min
n_eff <- nrow(jack_estimates)
b1_full_m5 <- result$par[1] / (1000 ^ result$par[3])
b2_full_m5 <- result$par[2] / (12   ^ result$par[4])

b1_jk <- jack_mean_pseudo[1] / (1000 ^ jack_mean_pseudo[3])
b2_jk <- jack_mean_pseudo[2] / (12   ^ jack_mean_pseudo[4])

# Table of results
param_table_m5min <- data.frame(
  Parameter = c("beta1", "beta2", "alpha1", "alpha2"),
  Estimate_full_kmh = result$par,
  Estimate_jk_kmh = jack_mean_pseudo,
  CI_lower_kmh = jackknife_monthyear_results$CI_lower,
  CI_upper_kmh = jackknife_monthyear_results$CI_upper,
  Estimate_full_m5min = c(params_m5min$beta1, params_m5min$beta2,
                          result$par[3], result$par[4]),
  Estimate_jk_m5min = c(b1_jk, b2_jk,
                          jack_mean_pseudo[3], jack_mean_pseudo[4]),
  CI_lower_m5min = c(beta1_m5_low, beta2_m5_low,
                     jackknife_monthyear_results$CI_lower[3],
                     jackknife_monthyear_results$CI_lower[4]),
  CI_upper_m5min = c(beta1_m5_high, beta2_m5_high,
                     jackknife_monthyear_results$CI_upper[3],
                     jackknife_monthyear_results$CI_upper[4])
)
print(param_table_m5min)

number_of_episodes <- length(list_episodes_filtered)
cat("Number of episodes used in optimization:", number_of_episodes, "\n")
adv_threshold <- speed_threshold
cat("Advection speed threshold (km/h):", adv_threshold, "\n")

# number of advection put to 0
n_zero_adv <- sum(V_episodes_filtered$vx == 0 & V_episodes_filtered$vy == 0)
cat("Number of episodes with zero advection after filtering:", n_zero_adv, "\n")


omsev_params <- params_est

dist_mat <- get_dist_mat(location_gauges) / 1000
df_dist <- reshape_distances(dist_mat)
n_hbins <- 12
h_all <- df_dist$value

h_breaks_omsev <- quantile(
  h_all,
  probs = seq(0, 1, length.out = n_hbins + 1),
  na.rm = TRUE
)

h_breaks_omsev <- unique(as.numeric(h_breaks_omsev))
h_breaks_omsev[length(h_breaks_omsev)] <- 1.6

chi_omsev <- compute_group_chi(
  list_lags = list_lags_filtered,
  list_excesses = list_excesses_filtered,
  wind_df = V_episodes_filtered,
  params = omsev_params,
  tau_fixed = 0,
  h_breaks = h_breaks_omsev,
  adv_transform = TRUE
)

corr_omsev <- summary_correlation(chi_omsev$res_cmp)
message(sprintf("OMSEV chi Spearman correlation: %.3f", corr_omsev))

# plot chi_come emp vs theo
print(chi_omsev$plots$all)

res_om_cmp <- chi_omsev$res_cmp
ggplot(res_om_cmp, aes(x = chi_theo_bar, y = chi_emp_bar, size = n_pairs)) +
  geom_point(alpha = 0.7, color = btfgreen) +
  geom_abline(linetype = "dashed", color = "red") +
  theme_bw() +
  scale_size_continuous(range = c(1, 4)) +
  labs(x = TeX("Theoretical $\\chi$"), y = TeX("Empirical $\\chi$"),
      size = "Number of pairs") +
  btf_theme +
  theme(legend.position = "right",
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 14))

# save plot
foldername <- paste0(im_folder, "workflows/full_pipeline/2025/")
if (!dir.exists(foldername)) {
  dir.create(foldername, recursive = TRUE)
}
filename <- paste0(foldername, "omsev_chi_th_emp_dmin1200.png")
ggsave(filename,
       width = 7,
       height = 6)



# --- REMPLACE TON BLOC DE CALCUL DE L'IC PAR CELUI-CI ---

jack_estimates <- do.call(rbind, jack_estimates_list)
jack_estimates <- na.omit(jack_estimates)
n_eff <- nrow(jack_estimates)

# 1. Passage à l'échelle log pour stabiliser la variance et forcer la positivité
log_jack_estimates <- log(jack_estimates)
log_init_param_jk  <- log(init_param_jk[1:4])

# 2. Calcul des pseudo-valeurs sur l'échelle LOG
# Note : Il est plus rigoureux d'utiliser la vraie estimation totale (log_init_param_jk) plutôt que la moyenne du jackknife
pseudo_values_log <- matrix(NA, nrow = n_eff, ncol = 4)
for (i in 1:n_eff) {
  pseudo_values_log[i, ] <- n_eff * log_init_param_jk - (n_eff - 1) * log_jack_estimates[i, ]
}

# 3. Moyenne et Standard Error sur l'échelle LOG
jack_mean_pseudo_log <- colMeans(pseudo_values_log)
jack_se_log          <- apply(pseudo_values_log, 2, sd) / sqrt(n_eff)

# 4. Calcul de l'IC sur l'échelle LOG
z <- qnorm(0.975)
lower_ci_log <- jack_mean_pseudo_log - z * jack_se_log
upper_ci_log <- jack_mean_pseudo_log + z * jack_se_log

# 5. RETOUR À L'ÉCHELLE NATURELLE VIA EXP() -> Plus de valeurs négatives possibles !
estimate_jk_naturel <- exp(jack_mean_pseudo_log) # Estimation corrigée du biais
lower_ci <- exp(lower_ci_log)
upper_ci <- exp(upper_ci_log)

jackknife_monthyear_results <- data.frame(
  Parameter     = c("beta1", "beta2", "alpha1", "alpha2"),
  Estimate_full = result$par,
  Estimate_jk   = estimate_jk_naturel,
  StdError_log  = jack_se_log, # L'erreur standard est maintenant relative (multiplicative)
  CI_lower      = lower_ci,
  CI_upper      = upper_ci
)
