library(data.table)
library(latex2exp)
library(lubridate)
library(fuzzyjoin)
library(grid)

muse <- TRUE

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
#                          "omsev/omsev_5min/rain_mtp_5min_2019_2024.RData")

# load(filename_omsev)

filename_loc <- paste0(data_folder,
                           "omsev/loc_rain_gauges.csv")
# get location of each rain gauge
location_gauges <- read.csv(filename_loc)
location_gauges$Station <- c("iem", "mse", "poly", "um", "cefe", "cnrs",
                             "crbm", "archiw", "archie", "um35", "chu1",
                             "chu2", "chu3", "chu4", "chu5", "chu6", "chu7",
                             "cines", "brives", "hydro")

# get rain data
filename_rain <- paste0(data_folder, "omsev/omsev_5min/rain_mtp_5min_2019_2024_cleaned.csv")
rain_omsev <- read.csv(filename_rain)


rain_test <- rain_omsev
rain_test$nb_sites_non_NA <- apply(rain_test[ , -1], 1, function(x) sum(!is.na(x)))
date_2_sites <- rain_test$dates[which(rain_test$nb_sites_non_NA >= 2)[1]]
date_3_sites <- rain_test$dates[which(rain_test$nb_sites_non_NA >= 3)[1]]
date_4_sites <- rain_test$dates[which(rain_test$nb_sites_non_NA >= 4)[1]]
date_5_sites <- rain_test$dates[which(rain_test$nb_sites_non_NA >= 5)[1]]

# begin in 2020-01-01
rain_omsev <- rain_omsev[rain_omsev$dates >= date_5_sites, ]
head(rain_omsev)

# put dates as rownames
rownames(rain_omsev) <- rain_omsev$dates
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

################################################################################
# WLSE results -----------------------------------------------------------------
################################################################################
# foldername <- paste0(data_folder, "omsev/WLSE/")
# df_spa_result <- read.csv(paste0(foldername, "wlse_spa_summary.csv"))
# df_temp_result <- read.csv(paste0(foldername, "wlse_temp_summary.csv"))
# df_result_all <- merge(df_spa_result, df_temp_result, by = NULL)

# # select one row
# df_result <- df_result_all[df_result_all$q_spa == 0.95 &
#                              df_result_all$q_temp == 0.95, ]

# beta1 <- df_result$beta1
# beta2 <- df_result$beta2
# alpha1 <- df_result$alpha1
# alpha2 <- df_result$alpha2

# # estimates
# wlse_omsev <- c(beta1, beta2, alpha1, alpha2) # m / 5 min

################################################################################
# VARIOGRAM --------------------------------------------------------------------
################################################################################

################################################################################

# in rain remove when all data are NA
set_st_excess <- get_spatiotemp_excess(rain, quantile = q, remove_zeros = TRUE)

# verify that the excess is above the threshold
# get list of sites and times
list_s <- set_st_excess$list_s
list_t <- set_st_excess$list_t
list_u <- set_st_excess$list_u

# Spatio-temporal neighborhood parameters
s0t0_set <- get_s0t0_pairs(grid_coords_km, rain,
                            min_spatial_dist = min_spatial_dist,
                            episode_size = delta,
                            set_st_excess = set_st_excess,
                            n_max_episodes = 10000,
                            latlon = FALSE)

selected_points <- s0t0_set

selected_points <- selected_points %>%
  mutate(t0_date = as.POSIXct(t0_date, format = "%Y-%m-%d %H:%M:%S", tz = "UTC"))

# Make sure rain is in matrix form
# site names in columns, time in rows
# column names must match `s0t0_set$s0`
stopifnot(all(s0t0_set$s0 %in% colnames(rain)))

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
                                     episode_size = delta, unif = FALSE,
                                     beta = 0)
list_episodes <- list_episodes_points$episodes

selected_points$t0_date <- as.POSIXct(selected_points$t0_date, format = "%Y-%m-%d %H:%M:%S", tz = "UTC")

# Round to the next hour
selected_points$t0_date_rounded <- ceiling_date(selected_points$t0_date, "hour")

adv_filename <- paste(data_folder, "/omsev/adv_estim/bary_omsev/episode_advection_q",
                          q * 100, "_delta", delta, "_dmin", min_spatial_dist,
                          ".csv", sep = "")
adv_df <- read.csv(adv_filename, sep = ",")
head(adv_df)

# # count 0 adv
# n_zero_adv <- sum(adv_df$dx_comb_kmh == 0 & adv_df$dy_comb_kmh == 0)
# cat("Number of episodes with zero advection:", n_zero_adv, "\n")

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
adv_df <- adv_df[matching_indices, ]
rownames(adv_df) <- NULL  # reset row names

selected_episodes <- selected_points
selected_episodes$adv_x <- rep(NA, nrow(selected_episodes))
selected_episodes$adv_y <- rep(NA, nrow(selected_episodes))

# get adv values for each episode according to the t0_date

library(data.table)
setDT(selected_points)
setDT(adv_df)
setkey(adv_df, t0_omsev)
selected_episodes <- adv_df[selected_points, roll = 5*60, on = .(t0_omsev = t0_date)]
colnames(selected_episodes)[which(names(selected_episodes) == "mean_dx_kmh_omsev")] <- "adv_x"
colnames(selected_episodes)[which(names(selected_episodes) == "mean_dy_kmh_omsev")] <- "adv_y"
V_episodes <- data.frame(
  v_x = selected_episodes$adv_x,
  v_y = selected_episodes$adv_y
)

colnames(V_episodes) <- c("vx", "vy")

# remove NA
V_episodes <- V_episodes[!is.na(V_episodes$vx) & !is.na(V_episodes$vy), ]

# number of 0 adv
n_zero_adv <- sum(V_episodes$vx == 0 & V_episodes$vy == 0)
cat("Number of episodes with zero advection:", n_zero_adv, "\n")

# # put 0 for NA adv values
# V_episodes$vx[is.na(V_episodes$vx)] <- 0
# V_episodes$vy[is.na(V_episodes$vy)] <- 0

tau_vect <- 0:10
# thresholds_by_site <- apply(rain, 2, function(col) {
#   col <- col[!is.na(col) & col > 0]
#   if (length(col) < 30) return(NA)
#   quantile(col, probs = q)
# })

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
  col_s0 <- which(colnames(rain) == s0)
  s0_coords <- df_coords[col_s0, ]
  episode <- list_episodes[[i]]
  ind_t0_ep <- 0 # index of t0 in the episode
  lags <- get_conditional_lag_vectors(df_coords, s0_coords, ind_t0_ep,
                                tau_vect, latlon = FALSE)
  u <- u0_list[i]
  excesses <- empirical_excesses_rpar(episode, quantile = u, threshold = TRUE,
                                      df_lags = lags, t0 = ind_t0_ep)
  # tau is in 5 minutes
  lags$tau <- lags$tau * 5 / 60 # convert to hours
  list(lags = lags, excesses = excesses)
}, mc.cores = detectCores() - 1)

list_lags <- lapply(list_results, `[[`, "lags")
list_excesses <- lapply(list_results, `[[`, "excesses")
df_lags <- list_lags[[1]]
df_excesses <- list_excesses[[1]]
sum(df_excesses$kij)


# get comephore estimates
filename_com_res <- paste(data_folder, 
        "/comephore/optim_results/lalpha/free_eta/combined_optim_results.csv",
        sep = "")
com_results <- read.csv(filename_com_res)

# get estimates from comephore
params_com <- com_results[com_results$q == q*100 &
                           com_results$delta == 30 &
                           com_results$dmin == 5, ]
init_params_com <- c(params_com$beta1, params_com$beta2,
                     params_com$alpha1, params_com$alpha2,
                     params_com$eta1, params_com$eta2)
# check for na in adv and wind
V_episodes <- V_episodes[!is.na(V_episodes$vx) & !is.na(V_episodes$vy), ]
hmax <- max(dist_mat) / 1000 # convert to km


# OPTIMISATION ---------------------------------------------------------------
# Composite likelihood optimisation with fixed eta1 and eta2
# only beta1, beta2, alpha1, alpha2 are estimated
init_params_com <- c(params_com$beta1, params_com$beta2,
                     params_com$alpha1, params_com$alpha2,
                     params_com$eta1, params_com$eta2)

# try with 0 adv
V_episodes_noadv <- V_episodes
V_episodes_noadv$vx <- rep(0, nrow(V_episodes_noadv))
V_episodes_noadv$vy <- rep(0, nrow(V_episodes_noadv))

V_episodes_meanadv <- V_episodes
V_episodes_meanadv$vx <- rep(mean(V_episodes$vx, na.rm = TRUE), nrow(V_episodes_meanadv))
V_episodes_meanadv$vy <- rep(mean(V_episodes$vy, na.rm = TRUE), nrow(V_episodes_meanadv))


result <- optim(
  par = init_params_com[1:4],
  fn = neg_ll_composite_fixed_eta,
  list_lags = list_lags,
  list_episodes = list_episodes,
  list_excesses = list_excesses,
  hmax = hmax,
  wind_df = V_episodes,
  latlon = FALSE,
  distance = "lalpha",
  fixed_eta1 = init_params_com[5],
  fixed_eta2 = init_params_com[6],
  method = "L-BFGS-B",
  lower = c(1e-08, 1e-08, 1e-08, 1e-08),
  upper = c(10, 10, 1.999, 1.999),
  control = list(maxit = 20000, trace = 1)
)

result

grad(objfun, result$par)

objfun(result$par)
objfun(result$par + c(1e-4, 0, 0, 0))
objfun(result$par - c(1e-4, 0, 0, 0))


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

filename <- paste0(foldername_res, "/results_q",
           q * 100, "_delta", delta, "_dmin", min_spatial_dist, 
           starting_year, ".csv")
write.csv(result_df, filename, row.names = FALSE)

# Jackknife CI ---------------------------------------------------------------

# Initial parameters for jackknife (from full data optimization)
init_param_jk <- as.numeric(result_df)
selected_episodes$t0_date <- as.POSIXct(selected_episodes$t0_date,
                                        format="%Y-%m-%d %H:%M:%S", tz="UTC")

# Define season boundaries (using day-of-year)
get_season <- function(date) {
  yday <- as.integer(format(date, "%j"))
  if (yday >= 80 & yday <= 171) {
    return("Spring")
  } else if (yday >= 172 & yday <= 263) {
    return("Summer")
  } else if (yday >= 264 & yday <= 354) {
    return("Autumn")
  } else {
    return("Winter")
  }
}

season_vec <- sapply(selected_episodes$t0_date, get_season)
year_vec <- format(selected_episodes$t0_date, "%Y")
season_year_vec <- paste(season_vec, year_vec)

unique_season_years <- sort(unique(season_year_vec))
n_season_years <- length(unique_season_years)

jack_estimates_list <- parallel::mclapply(unique_season_years, function(season_year) {
  cat("Excluding season-year:", season_year, "\n")
  exclude_idx <- which(season_year_vec == season_year)
  jack_episodes <- episodes_opt[-exclude_idx]
  jack_lags <- lags_opt[-exclude_idx]
  jack_excesses <- excesses_opt[-exclude_idx]
  jack_wind <- V_episodes[-exclude_idx, , drop = FALSE]
  res <- tryCatch({
    optim(par = init_param_jk[1:4], fn = neg_ll_composite_fixed_eta,
      list_lags = jack_lags, list_episodes = jack_episodes,
      list_excesses = jack_excesses, hmax = hmax,
      wind_df = jack_wind,
      latlon = FALSE,
      distance = "lalpha",
      fixed_eta1 = params_com$eta1,
      fixed_eta2 = params_com$eta2,
      method = "L-BFGS-B",
      lower = c(1e-08, 1e-08, 1e-08, 1e-08),
      upper = c(10, 10, 1.999, 1.999),
      control = list(maxit = 10000),
      hessian = FALSE)
  }, error = function(e) NULL)
  if (!is.null(res)) {
    return(res$par)
  } else {
    return(rep(NA, length(init_param_jk)))
  }
}, mc.cores = ncores)

jack_estimates <- do.call(rbind, jack_estimates_list)
jack_estimates <- na.omit(jack_estimates)
n_eff <- nrow(jack_estimates)

filename <- paste0(data_folder, "omsev/optim_results/jackknife_estimates/all_results_jk_by_seasonyear_n", 
           n_eff, "_q", q*100, "_delta", delta, "_dmin", min_spatial_dist, ".csv")
write.csv(jack_estimates, filename, row.names = FALSE)

jack_mean <- colMeans(jack_estimates)
pseudo_values <- matrix(NA, nrow = n_eff, ncol = length(init_param_jk))
for (i in 1:n_eff) {
  pseudo_values[i, ] <- n_eff * result$par - (n_eff - 1) * jack_estimates[i, ]
}

jack_mean_pseudo <- colMeans(pseudo_values)
jack_se <- apply(pseudo_values, 2, sd) / sqrt(n_eff)

z <- qnorm(0.975)
lower_ci <- jack_mean_pseudo - z * jack_se
upper_ci <- jack_mean_pseudo + z * jack_se

jackknife_seasonyear_results <- data.frame(
  Parameter = c("beta1", "beta2", "alpha1", "alpha2", "eta1", "eta2"),
  Estimate_full = result$par,
  Estimate_jk   = jack_mean_pseudo,
  StdError = jack_se,
  CI_lower = lower_ci,
  CI_upper = upper_ci
)


filename <- paste0(data_folder, "omsev/optim_results/jackknife_estimates/results_jk_by_seasonyear_n", 
           n_eff, "_q", q*100, "_delta", delta, "_dmin", min_spatial_dist, ".csv")
write.csv(jackknife_seasonyear_results, filename, row.names = FALSE)
