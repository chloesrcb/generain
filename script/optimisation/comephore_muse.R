# remove objects from the environment
rm(list = ls())

# clear console
cat("\014")

muse <- TRUE

if (muse) {
  # Get the muse folder
  folder_muse <- "/home/serrec/work_rainstsimu/rpareto/"
  setwd(folder_muse)
  # Load libraries and set theme
  source("load_libraries.R")
  im_folder <- "./images"
  source("config_com.R")
  data_folder <- "./data/"
  ncores <- 27
} else {
  # Load libraries and set theme
  source("./script/load_libraries.R")
  source("./script/optimisation/config_com.R")
  ncores <- detectCores() - 1
}

# Get all files in the folder "R"
functions_folder <- "./R"
files <- list.files(functions_folder, full.names = TRUE)
invisible(lapply(files, function(f) source(f, echo = FALSE)))

library(latex2exp)
library(sf)
library(dplyr)
library(ggplot2)
library(gridExtra)

eta_type <- if (!is.na(fixed_eta1) && !is.na(fixed_eta2)) {
  "fixed_eta"
} else if (!is.na(fixed_eta1)) {
  "fixed_eta1"
} else if (!is.na(fixed_eta2)) {
  "fixed_eta2"
} else {
  "free_eta"
}

# Folder name to save the data
foldername_res <- file.path(
  paste0(data_folder, "comephore/optim_results"),
  distance_type,
  eta_type
)

if (!dir.exists(foldername_res)) {
  message("Folder created: ", foldername_res)
  dir.create(foldername_res, recursive = TRUE)
}

# LOAD DATA ####################################################################
filename_com <- paste0(data_folder, "comephore/comephore_2008_2024_5km.csv")
comephore_raw <- read.csv(filename_com, sep = ",")
filename_loc <- paste0(data_folder, "comephore/loc_px_zoom_5km.csv")
loc_px <- read.csv(filename_loc, sep = ",")

# remove pixel in loc_px that are not in comephore_raw
loc_px <- loc_px[loc_px$pixel_name %in% colnames(comephore_raw)[-1], ]
rownames(loc_px) <- NULL

df_comephore <- as.data.frame(comephore_raw)
# colnames(df_comephore)[1] <- "date"
df_comephore$date <- as.POSIXct(df_comephore$date,
                                format = "%Y-%m-%d %H:%M:%S", tz = "UTC")
rownames(df_comephore) <- format(df_comephore$date, "%Y-%m-%d %H:%M:%S")
comephore <- df_comephore[-1] # remove dates column

# DISTANCE AND COORDS ##########################################################
# Get number of sites
nsites <- nrow(loc_px)

# Get coords
sites_coords <- loc_px[, c("Longitude", "Latitude")]
rownames(sites_coords) <- loc_px$pixel_name
sites_coords_sf <- st_as_sf(sites_coords, coords = c("Longitude", "Latitude"),
                            crs = 4326)
sites_coords_sf <- st_transform(sites_coords_sf, crs = 2154)
coords_m <- st_coordinates(sites_coords_sf)
grid_coords_km <- as.data.frame(coords_m / 1000)
colnames(grid_coords_km) <- c("Longitude", "Latitude")
rownames(grid_coords_km) <- rownames(sites_coords)

# Spatial chi WLSE #############################################################
foldername <- paste0(data_folder, "/comephore/WLSE/")

# get csv
df_spa_result <- read.csv(paste0(foldername, "wlse_spa_summary.csv"))
df_temp_result <- read.csv(paste0(foldername, "wlse_temp_summary.csv"))
df_result_all <- merge(df_spa_result, df_temp_result, by = NULL)

df_result_all <- df_result_all[, c("q_spa", "q_temp",
                                   "beta1", "alpha1", "beta2", "alpha2")]

# round values to 4 decimal places
df_result_all$beta1 <- round(df_result_all$beta1, 4)
df_result_all$alpha1 <- round(df_result_all$alpha1, 4)
df_result_all$beta2 <- round(as.numeric(df_result_all$beta2), 4)
df_result_all$alpha2 <- round(as.numeric(df_result_all$alpha2), 4)

# choice of one row
df_result <- df_result_all[df_result_all$q_spa == q &
                           df_result_all$q_temp == q, ]

beta1 <- df_result$beta1
beta2 <- df_result$beta2
alpha1 <- df_result$alpha1
alpha2 <- df_result$alpha2
# CHOOSE EXTREME EPISODE FOR R-PARETO ##########################################

# get central site from sites_coords
set_st_excess <- get_spatiotemp_excess(data = comephore, quantile = q,
                                      remove_zeros = TRUE)
first_ts <- as.POSIXct(rownames(comephore)[1], tz = "UTC")
starting_year <- year(first_ts)
starting_year_suffix <- if (starting_year == 2008) "" else paste0("_from", starting_year)

list_s <- set_st_excess$list_s
list_t <- set_st_excess$list_t
list_u <- set_st_excess$list_u

# check that we have excess 
for (i in seq_along(list_s)) {
  s0 <- list_s[[i]]                     # station name as string
  t0 <- list_t[[i]][1]                  # numeric index
  u_s0 <- list_u[[i]][1]                # numeric threshold

  # Check if the value at (t0, s0) is above threshold
  if (comephore[t0, s0] <= u_s0) {
    message <- paste("Excess is not above threshold for s =", s0, "and t =", t0)
    print(message)
    stop()
  }
}

# Spatio-temporal neighborhood parameters
s0t0_set <- get_s0t0_pairs(grid_coords_km, comephore,
                            min_spatial_dist = min_spatial_dist,
                            episode_size = delta,
                            set_st_excess = set_st_excess,
                            n_max_episodes = 10000,
                            latlon = FALSE,
                            beta = 0)

selected_points <- s0t0_set

n_episodes <- length(selected_points$s0)
t0_list <- selected_points$t0
s0_list <- selected_points$s0
u_list <- selected_points$u_s0
list_episodes_points <- get_extreme_episodes(selected_points, comephore,
                                     episode_size = delta, unif = FALSE)

list_episodes <- list_episodes_points$episodes
s0_list <- selected_points$s0
t0_list <- selected_points$t0
u_list <- selected_points$u_s0
episode <- list_episodes[[1]]

library(parallel)

tau_vect <- 0:10
tmax <- max(tau_vect)
df_coords <- as.data.frame(grid_coords_km)

# Compute the lags and excesses for each conditional point
list_results <- parallel::mclapply(1:length(s0_list), function(i) {
  s0 <- s0_list[i]
  s0_coords <- df_coords[s0, ]
  u <- u_list[i]
  episode <- list_episodes[[i]]
  ind_t0_ep <- 0 # index of t0 in the episode
  lags <- get_conditional_lag_vectors(df_coords, s0_coords, ind_t0_ep,
                                tau_vect, latlon = FALSE)
  excesses <- empirical_excesses_rpar(episode, threshold = u, # !!!!!!!
                                  df_lags = lags, t0 = ind_t0_ep)
  list(lags = lags, excesses = excesses)
}, mc.cores = detectCores() - 1)

list_lags <- lapply(list_results, `[[`, "lags")
list_excesses <- lapply(list_results, `[[`, "excesses")

# ADD ADVECTION ESTIMATES ######################################################
adv_filename <- paste0(
  data_folder, "comephore/adv_estim/advection_results_q",
  q * 100, "_delta", delta, "_dmin", min_spatial_dist,
  starting_year_suffix, ".csv"
)
adv_df <- read.csv(adv_filename, sep = ",")
adv_df$t0 <- as.POSIXct(adv_df$t0, tz = "UTC")
selected_points$t0_date <- as.POSIXct(selected_points$t0_date, tz = "UTC")

matching_indices <- match(selected_points$t0_date, adv_df$t0)
matching_indices <- matching_indices[!is.na(matching_indices)]
adv_df <- adv_df[matching_indices, ]
rownames(adv_df) <- NULL

selected_episodes <- selected_points
selected_episodes$adv_x <- NA
selected_episodes$adv_y <- NA

for (i in 1:nrow(selected_episodes)) {
  t0_date <- selected_episodes$t0_date[i]
  adv_row <- adv_df[adv_df$t0 == t0_date, ]
  if (nrow(adv_row) > 0) {
    selected_episodes$adv_x[i] <- adv_row$mean_dx_kmh[1]
    selected_episodes$adv_y[i] <- adv_row$mean_dy_kmh[1]
  }
}


wind_df <- data.frame(
  vx = selected_episodes$adv_x, # * 3.6,  # m/s -> km/h
  vy = selected_episodes$adv_y # * 3.6
)
head(wind_df)

# OPTIMIZATION #################################################################
hmax <- 7
stopifnot(length(list_episodes) == nrow(wind_df),
          length(list_lags)     == nrow(wind_df),
          length(list_excesses) == nrow(wind_df))
lags <- list_lags[[1]]
excesses <- list_excesses[[1]]

print("Starting optimization")
init_param <- c(beta1, beta2, alpha1, alpha2, 1, 1)
result <- optim(par = init_param, fn = neg_ll_composite,
          list_lags = list_lags, list_episodes = list_episodes,
          list_excesses = list_excesses, hmax = NA,
          wind_df = wind_df,
          latlon = FALSE,
          distance = distance_type,
          method = "L-BFGS-B",
          lower = c(1e-08, 1e-08, 1e-08, 1e-08, 1e-08, 1e-08),
          upper = c(10, 10, 1.999, 1.999, 10, 10),
          control = list(maxit = 10000, trace = 1))

# tau = 0:10, time window adv = +/-2 hours
# hmax = NA
# 1] 0.3247256 0.6725507 0.3863810 0.7220611 1.7808056 5.0345370
# final  value 413678.487943 




# } else if (eta_type == "fixed_eta1") {
#   init_param <- c(beta1, beta2, alpha1, alpha2, 1)
#   result <- optim(par = init_param, fn = neg_ll_composite_fixed_eta,
#           list_lags = list_lags, list_episodes = list_episodes,
#           list_excesses = list_excesses, hmax = hmax,
#           wind_df = wind_df,
#           latlon = FALSE,
#           distance = distance_type,
#           fixed_eta1 = fixed_eta1,
#           method = "L-BFGS-B",
#           lower = c(1e-08, 1e-08, 1e-08, 1e-08, 1e-08),
#           upper = c(10, 10, 1.999, 1.999, 10),
#           control = list(maxit = 10000, trace = 1),
#           hessian = FALSE)
# } else if (eta_type == "fixed_eta2") {
#   init_param <- c(beta1, beta2, alpha1, alpha2, 1)
#   result <- optim(par = init_param, fn = neg_ll_composite_fixed_eta,
#           list_lags = list_lags, list_episodes = list_episodes,
#           list_excesses = list_excesses, hmax = hmax,
#           wind_df = wind_df,
#           latlon = FALSE,
#           distance = distance_type,
#           fixed_eta2 = fixed_eta2,
#           method = "L-BFGS-B",
#           lower = c(1e-08, 1e-08, 1e-08, 1e-08, 1e-08),
#           upper = c(10, 10, 1.999, 1.999, 10),
#           control = list(maxit = 10000, trace = 1),
#           hessian = FALSE)
#   } else if (eta_type == "fixed_eta") {
#   init_param <- c(beta1, beta2, alpha1, alpha2)
#   result <- optim(par = init_param, fn = neg_ll_composite_fixed_eta,
#           list_lags = list_lags, list_episodes = list_episodes,
#           list_excesses = list_excesses, hmax = hmax,
#           wind_df = wind_df,
#           latlon = FALSE,
#           distance = distance_type,
#           fixed_eta1 = fixed_eta1,
#           fixed_eta2 = fixed_eta2,
#           method = "L-BFGS-B",
#           lower = c(1e-08, 1e-08, 1e-08, 1e-08),
#           upper = c(10, 10, 1.999, 1.999),
#           control = list(maxit = 10000, trace = 1),
#           hessian = FALSE)
# }
if (result$convergence != 0) {
  stop("Optimization did not converge")
} else {
  print("Optimization converged")
}
# SAVE RESULTS #################################################################
result_df <- data.frame(
  beta1 = result$par[1],
  beta2 = result$par[2],
  alpha1 = result$par[3],
  alpha2 = result$par[4],
  eta1 = result$par[5],
  eta2 = result$par[6],
  nll = result$value
)

# Folder name to save the data
foldername_res <- file.path(
  paste0(data_folder, "comephore/optim_results"),
  distance_type,
  eta_type
)
filename <- paste0(foldername_res, "/results_q",
           q * 100, "_delta", delta, "_dmin", min_spatial_dist,
           starting_year_suffix, ".csv")

write.csv(result_df, filename, row.names = FALSE)
print(result_df)
cat("Results saved to", filename, "\n")

# # JACKKNIFE CI #################################################################
# Initial parameters for jackknife (from full data optimization)
init_param_jk <- as.numeric(result_df)[-1]
selected_episodes$t0_date <- as.POSIXct(selected_episodes$t0_date, format="%Y-%m-%d %H:%M:%S", tz="UTC")

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
  jack_episodes <- list_episodes[-exclude_idx]
  jack_lags <- list_lags[-exclude_idx]
  jack_excesses <- list_excesses[-exclude_idx]
  jack_wind <- wind_df[-exclude_idx, , drop = FALSE]
  res <- tryCatch({
    optim(par = init_param_jk, fn = neg_ll_composite,
      list_lags = jack_lags, list_episodes = jack_episodes,
      list_excesses = jack_excesses, hmax = hmax,
      wind_df = jack_wind,
      latlon = FALSE,
      distance = distance_type,
      method = "L-BFGS-B",
      lower = c(1e-08, 1e-08, 1e-08, 1e-08, 1e-08, 1e-08),
      upper = c(10, 10, 1.999, 1.999, 10, 10),
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

foldername <- file.path(
  paste0(data_folder, "comephore/optim_results"),
  distance_type,
  eta_type,
  "jackknife_estimates"
)
filename <- paste0(foldername, "/all_results_jk_by_seasonyear_n", 
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


filename <- paste0(foldername, "/results_jk_by_seasonyear_n", 
           n_eff, "_q", q*100, "_delta", delta, "_dmin", min_spatial_dist, ".csv")
write.csv(jackknife_seasonyear_results, filename, row.names = FALSE)
