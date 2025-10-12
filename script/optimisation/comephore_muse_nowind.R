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
# # load all functions in files
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
  paste0(data_folder, "comephore/optim_results/"),
  distance_type,
  eta_type,
  "no_wind"
)

if (!dir.exists(foldername_res)) {
  message("Folder created: ", foldername_res)
  dir.create(foldername_res, recursive = TRUE)
}

# LOAD DATA ####################################################################
filename_com <- paste0(data_folder, "comephore/comephore_2008_2024_5km.csv")
# filename_com <- paste0(data_folder, "comephore/zoom_5km.csv")
comephore_raw <- read.csv(filename_com, sep = ",")
# filename_loc <- paste0(data_folder, "comephore/coords_pixels_10km.csv")
filename_loc <- paste0(data_folder, "comephore/loc_px_zoom_5km.csv")
loc_px <- read.csv(filename_loc, sep = ",")

# remove pixel in loc_px that are not in comephore_raw
loc_px <- loc_px[loc_px$pixel_name %in% colnames(comephore_raw)[-1], ]
# reindex
rownames(loc_px) <- NULL
df_comephore <- as.data.frame(comephore_raw)
df_comephore$date <- as.POSIXct(df_comephore$date,
                                format = "%Y-%m-%d %H:%M:%S", tz = "UTC")

# Take only data after 2007
df_comephore <- df_comephore[df_comephore$date >= "2008-01-01", ]
rownames(df_comephore) <- format(as.POSIXct(df_comephore$date), "%Y-%m-%d %H:%M:%S")
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
df_dist_km$value <- round(df_dist$value / 1000, 1)


# Spatial chi WLSE #############################################################
h_vect <- sort(unique(df_dist_km$value))
h_vect <- h_vect[h_vect > 0]  # remove 0
hmax <- h_vect[10] # 10th value
hmax <- 10
tmax <- 10

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
comephore_subset <- comephore[rownames(comephore) >= "2008-01-01", ]
set_st_excess <- get_spatiotemp_excess(data = comephore_subset, quantile = q,
                                      remove_zeros = TRUE)
starting_year <- as.character(year(rownames(comephore_subset)[1]))
if (starting_year == "2008") {
  starting_year <- ""
} else {
  starting_year <- paste0("_from", starting_year)
}

# remove date column if exists
if ("date" %in% colnames(comephore_subset)) {
  comephore_subset <- comephore_subset[,
                                -c(which(colnames(comephore_subset) == "date"))]
}

list_s <- set_st_excess$list_s
list_t <- set_st_excess$list_t
list_u <- set_st_excess$list_u

# check that we have excess 
for (i in seq_along(list_s)) {
  s0 <- list_s[[i]]                     # station name as string
  t0 <- list_t[[i]][1]                  # numeric index
  u_s0 <- list_u[[i]][1]                # numeric threshold

  # Check if the value at (t0, s0) is above threshold
  if (comephore_subset[t0, s0] <= u_s0) {
    message <- paste("Excess is not above threshold for s =", s0, "and t =", t0)
    print(message)
    stop()
  }
}

# unique(list_u) # check unique excess values
# Spatio-temporal neighborhood parameters
s0t0_set <- get_s0t0_pairs(grid_coords_km, comephore_subset,
                            min_spatial_dist = min_spatial_dist,
                            episode_size = delta,
                            set_st_excess = set_st_excess,
                            n_max_episodes = 10000,
                            latlon = FALSE)

selected_points <- s0t0_set

n_episodes <- length(selected_points$s0)
t0_list <- selected_points$t0
s0_list <- selected_points$s0
u_list <- selected_points$u_s0
list_episodes_points <- get_extreme_episodes(selected_points, comephore_subset,
                                     episode_size = delta, unif = FALSE,
                                     beta = 0)

list_episodes <- list_episodes_points$episodes
s0_list <- selected_points$s0
t0_list <- selected_points$t0
u_list <- selected_points$u_s0

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
  excesses <- empirical_excesses_rpar(episode, quantile = u,
                                  threshold = TRUE, # !!!!!!!
                                  df_lags = lags, t0 = ind_t0_ep)
  lags$tau <- lags$tau
  list(lags = lags, excesses = excesses)
}, mc.cores = detectCores() - 1)

list_lags <- lapply(list_results, `[[`, "lags")
list_excesses <- lapply(list_results, `[[`, "excesses")

# ADD ADVECTION ESTIMATES ######################################################

selected_episodes <- selected_points
selected_episodes$adv_x <- rep(0, nrow(selected_episodes))
selected_episodes$adv_y <- rep(0, nrow(selected_episodes))

# get adv values for each episode
ind_NA_adv <- which(is.na(selected_episodes$adv_x) | is.na(selected_episodes$adv_y))
selected_episodes_nona <- selected_episodes[-ind_NA_adv, ]

wind_df <- selected_episodes[, c("adv_x", "adv_y")]
colnames(wind_df) <- c("vx", "vy")
length(wind_df$vx) # should be the same as number of episodes

# OPTIMIZATION #################################################################

ind_NA <- which(is.na(wind_df$vx))
ind_NA <- ind_NA_adv
wind_opt <- wind_df
if (any(ind_NA > 0)) {
  # remove these episodes
  # wind_opt <- wind_df[-ind_NA, ]
  episodes_opt <- list_episodes[-ind_NA]
  lags_opt <- list_lags[-ind_NA]
  excesses_opt <- list_excesses[-ind_NA]
} else {
  wind_opt <- wind_df
  episodes_opt <- list_episodes
  lags_opt <- list_lags
  excesses_opt <- list_excesses
}



filename <- paste0(foldername_res, "/results_q",
           q * 100, "_delta", delta, "_dmin", min_spatial_dist, 
           starting_year, ".csv")

# do_optim <- TRUE
# if (file.exists(filename)) {
#   result_df <- read.csv(filename, sep = ",")
#   print("File exists, loading results")
#   do_optim <- FALSE
# }

init_param <- c(beta1, beta2, alpha1, alpha2, 1, 1)

print("Starting optimization")
result <- optim(par = init_param, fn = neg_ll_composite,
        list_lags = lags_opt, list_episodes = episodes_opt,
        list_excesses = excesses_opt, hmax = hmax,
        wind_df = wind_opt,
        latlon = FALSE,
        distance = "lalpha",
        method = "L-BFGS-B",
        lower = c(1e-08, 1e-08, 1e-08, 1e-08, 1e-08, 1e-08),
        upper = c(10, 10, 1.999, 1.999, 10, 10),
        control = list(maxit = 10000, trace = 1),
        hessian = FALSE)
if (result$convergence != 0) {
  stop("Optimization did not converge")
}
print("Optimization converged")
result_df <- data.frame(beta1 = result$par[1],
                        beta2 = result$par[2],
                        alpha1 = result$par[3],
                        alpha2 = result$par[4],
                        eta1 = result$par[5],
                        eta2 = result$par[6],
                        nll = result$value)
write.csv(result_df, filename, row.names = FALSE)



# JACKKNIFE CI #################################################################
# Initial parameters for jackknife (from full data optimization)
init_param_jk <- as.numeric(result_df)
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
  jack_episodes <- episodes_opt[-exclude_idx]
  jack_lags <- lags_opt[-exclude_idx]
  jack_excesses <- excesses_opt[-exclude_idx]
  jack_wind <- wind_opt[-exclude_idx, , drop = FALSE]
  res <- tryCatch({
    optim(par = init_param_jk, fn = neg_ll_composite,
      list_lags = jack_lags, list_episodes = jack_episodes,
      list_excesses = jack_excesses, hmax = hmax,
      wind_df = jack_wind,
      latlon = FALSE,
      distance = "lalpha",
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

foldername <- paste0(data_folder, "comephore/optim_results/jackknife_estimates/no_wind/")
if (!dir.exists(foldername)) {
  dir.create(foldername, recursive = TRUE)
}

filename <- paste0(data_folder, "comephore/optim_results/jackknife_estimates/no_wind/all_results_jk_by_seasonyear_n", 
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

filename <- paste0(data_folder, "comephore/optim_results/jackknife_estimates/no_wind/results_jk_by_seasonyear_n", 
           n_eff, "_q", q*100, "_delta", delta, "_dmin", min_spatial_dist, ".csv")
write.csv(jackknife_seasonyear_results, filename, row.names = FALSE)
