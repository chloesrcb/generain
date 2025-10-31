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
grid_coords_km <- as.data.frame(coords_m / 1000)
colnames(grid_coords_km) <- c("Longitude", "Latitude")
rownames(grid_coords_km) <- rownames(sites_coords)
# grid_coords_m <- sites_coords
# grid_coords_m$x_m <- (coords_m[, "X"] - min(coords_m[, "X"]))
# grid_coords_m$y_m <- (coords_m[, "Y"] - min(coords_m[, "Y"]))
# grid_coords_km$x_km <- (coords_m[, "X"] - min(coords_m[, "X"])) / 1000
# grid_coords_km$y_km <- (coords_m[, "Y"] - min(coords_m[, "Y"]))  / 1000

# # get distance matrix
# grid_coords_m <- grid_coords_m[, c("x_m", "y_m")]
# grid_coords_km <- grid_coords_km[, c("x_km", "y_km")]
# colnames(grid_coords_m) <- c("Longitude", "Latitude")
# colnames(grid_coords_km) <- c("Longitude", "Latitude")
# dist_mat <- get_dist_mat(grid_coords_m, latlon = FALSE)

# # Spatial chi
# df_dist <- reshape_distances(dist_mat)
# df_dist_km <- df_dist
# df_dist_km$value <- df_dist$value / 1000

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
length(selected_points$s0)
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
                                     episode_size = delta, unif = FALSE)
list_episodes <- list_episodes_points$episodes

selected_points$t0_date <- as.POSIXct(selected_points$t0_date, format = "%Y-%m-%d %H:%M:%S", tz = "UTC")

# Round to the next hour
selected_points$t0_date_rounded <- ceiling_date(selected_points$t0_date, "hour")



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


adv_filename <- paste(data_folder, "/omsev/adv_estim/combined_comephore_omsev/episode_advection_q",
                          q * 100, "_delta", delta, "_dmin", min_spatial_dist,
                          ".csv", sep = "")
adv_df <- read.csv(adv_filename, sep = ",")

head(adv_df)

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
tail(adv_df)
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
colnames(selected_episodes)[which(names(selected_episodes) == "vx_final")] <- "adv_x"
colnames(selected_episodes)[which(names(selected_episodes) == "vy_final")] <- "adv_y"
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
df_excesses <- list_excesses[[10]]
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
init_params_com <- c(0.2593, 0.5494, 0.5181, 0.9604, 1.4247, 4.8454)


# check for na in adv and wind
hmax <- max(dist_mat) / 1000


# OPTIMISATION ---------------------------------------------------------------
# Composite likelihood optimisation with fixed eta1 and eta2

# # try with 0 adv
# V_episodes_noadv <- V_episodes
# V_episodes_noadv$vx <- rep(0, nrow(V_episodes_noadv))
# V_episodes_noadv$vy <- rep(0, nrow(V_episodes_noadv))

# V_episodes_meanadv <- V_episodes
# V_episodes_meanadv$vx <- rep(mean(V_episodes$vx, na.rm = TRUE), nrow(V_episodes_meanadv))
# V_episodes_meanadv$vy <- rep(mean(V_episodes$vy, na.rm = TRUE), nrow(V_episodes_meanadv))


result <- optim(
  par = init_params_com[1:4],
  fn = neg_ll_composite_fixed_eta,
  list_lags = list_lags,
  list_episodes = list_episodes,
  list_excesses = list_excesses,
  hmax = hmax,
  wind_df = V_episodes,
  latlon = FALSE,
  distance = "euclidean",
  fixed_eta1 = init_params_com[5],
  fixed_eta2 = init_params_com[6],
  method = "L-BFGS-B",
  lower = c(1e-08, 1e-08, 1e-08, 1e-08),
  upper = c(10, 10, 1.999, 1.999),
  control = list(maxit = 20000, trace = 1)
)

result


# result_full <- optim(
#   par = init_params_com,
#   fn = neg_ll_composite,
#   list_lags = list_lags,
#   list_episodes = list_episodes,
#   list_excesses = list_excesses,
#   hmax = hmax,
#   wind_df = V_episodes,
#   latlon = FALSE,
#   distance = "euclidean",
#   method = "L-BFGS-B",
#   lower = c(1e-08, 1e-08, 1e-08, 1e-08, -1e-08, -1e-08),
#   upper = c(10, 10, 1.999, 1.999, 10, 10),
#   control = list(maxit = 20000, trace = 1)
# )


# result_full

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
filename <- paste0(foldername_res, "/results_q",
           q * 100, "_delta", delta, "_dmin", min_spatial_dist, ".csv")
write.csv(result_df, filename, row.names = FALSE)

# Jackknife CI ---------------------------------------------------------------

# Initial parameters for jackknife (from full data optimization)
init_param_jk <- as.numeric(result_df)

#### MONTH-YEAR JACKKNIFE ------------------------------------------------------

selected_episodes$t0_date <- as.POSIXct(selected_episodes$t0_date,
                                        format="%Y-%m-%d %H:%M:%S", tz="UTC")

month_vec <- format(selected_episodes$t0_date, "%m")
year_vec  <- format(selected_episodes$t0_date, "%Y")
month_year_vec <- paste(year_vec, month_vec, sep = "-")

unique_month_years <- sort(unique(month_year_vec))
n_month_years <- length(unique_month_years)

jack_estimates_list <- parallel::mclapply(unique_month_years, function(month_year) {
  cat("Excluding month-year:", month_year, "\n")
  exclude_idx <- which(month_year_vec == month_year)
  
  jack_episodes <- list_episodes[-exclude_idx]
  jack_lags <- list_lags[-exclude_idx]
  jack_excesses <- list_excesses[-exclude_idx]
  jack_wind <- V_episodes[-exclude_idx, , drop = FALSE]
  
  res <- tryCatch({
    optim(par = init_param_jk[1:4], fn = neg_ll_composite_fixed_eta,
      list_lags = jack_lags, list_episodes = jack_episodes,
      list_excesses = jack_excesses, hmax = hmax,
      wind_df = jack_wind,
      latlon = FALSE,
      distance = "euclidean",
      fixed_eta1 = init_params_com[5],
      fixed_eta2 = init_params_com[6],
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

foldername_jk <- paste0(data_folder, "omsev/optim_results/jackknife_estimates/")
if (!dir.exists(foldername_jk)) {
  dir.create(foldername_jk, recursive = TRUE)
}
filename <- paste0(foldername_jk, "all_results_jk_by_monthyear_n", 
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

jackknife_monthyear_results <- data.frame(
  Parameter = c("beta1", "beta2", "alpha1", "alpha2"),
  Estimate_full = result$par,
  Estimate_jk   = jack_mean_pseudo,
  StdError = jack_se,
  CI_lower = lower_ci,
  CI_upper = upper_ci
)

filename <- paste0(foldername_jk, "results_jk_by_monthyear_n",
           n_eff, "_q", q*100, "_delta", delta, "_dmin", min_spatial_dist, ".csv")
write.csv(jackknife_monthyear_results, filename, row.names = FALSE)

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

# Après avoir obtenu jack_estimates (colonnes: beta1, beta2, alpha1, alpha2)
beta1_m5min_i <- jack_estimates[,1] / (1000 ^ jack_estimates[,3])
beta2_m5min_i <- jack_estimates[,2] / (12   ^ jack_estimates[,4])

# bornes croisées correctes
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
