# remove objects from the environment
rm(list = ls())

# clear console
cat("\014")

# Get the muse folder
folder_muse <- "/home/serrec/work_rainstsimu/rpareto/"
setwd(folder_muse)
# Load libraries and set theme
source("load_libraries.R")
im_folder <- "./images"
source("config.R")

# Get all files in the folder "R"
functions_folder <- "./R"
files <- list.files(functions_folder, full.names = TRUE)
# # load all functions in files
invisible(lapply(files, function(f) source(f, echo = FALSE)))
library(latex2exp)
data_folder <- "./data/"
library(sf)
library(dplyr)
library(ggplot2)
library(gridExtra)

# LOAD DATA ####################################################################
filename_com <- paste0(data_folder, "comephore/comephore_2008_2024_5km.csv")
# filename_com <- paste0(data_folder, "comephore/zoom_5km.csv")
comephore_raw <- read.csv(filename_com, sep = ",")
# filename_loc <- paste0(data_folder, "comephore/coords_pixels_10km.csv")
filename_loc <- paste0(data_folder, "comephore/loc_px_zoom_5km.csv")
loc_px <- read.csv(filename_loc, sep = ",")

# remove pixel in loc_px that are not in comephore_raw
loc_px <- loc_px[loc_px$pixel_name %in% colnames(comephore_raw)[-1], ]
nrow(loc_px) # number of pixels
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

hmax <- 7

q_no0_spa <- 0.97
chispa_df <- spatial_chi_alldist(df_dist_km, data_rain = comephore,
                          quantile = q_no0_spa, hmax = hmax, zeros = FALSE)

# etachispa_df <- data.frame(chi = eta(chispa_df$chi),
#                            lagspa = log(chispa_df$lagspa))

chispa_df$lagspa <- chispa_df$lagspa * 1000
# WLSE
wlse_spa <- get_estimate_variospa(chispa_df, weights = "exp", summary = TRUE)

c1 <- wlse_spa[[1]]
beta1 <- wlse_spa[[2]]
alpha1 <- wlse_spa[[3]]


# Temporal chi WLSE ############################################################

tmax <- 10
q_no0_temp <- 0.95 # quantile for temporal chi

# Temporal chi with spatial lag fixed at 0
# compute chiplot of every site with itself but lagged in time
chimat_dtlag <- temporal_chi(comephore, quantile = q_no0_temp, tmax = tmax,
                             mean = FALSE, zeros = FALSE)

chi_df_dt <- data.frame(chimat_dtlag)
colnames(chi_df_dt) <- c(0:tmax) # temporal lags from 1 to tmax
rownames(chi_df_dt) <- c(1:nsites) # stations
# remove lag 0
chi_df_dt <- chi_df_dt[, -1] # remove lag 0

# Mean of chi
chimat_dt_mean <- temporal_chi(comephore, tmax, quantile = q_no0_temp,
                               mean = TRUE, zeros = FALSE)
df_chi <- data.frame(lag = c(0:tmax), chi = chimat_dt_mean)

# df_chi_not0 <- df_chi[df_chi$lag > 0, ]
wlse_temp <- get_estimate_variotemp(df_chi, weights = "exp", summary = TRUE)
c2 <- as.numeric(wlse_temp[[1]])
beta2 <- as.numeric(wlse_temp[[2]])
alpha2 <- as.numeric(wlse_temp[[3]])

# Result WLSE
df_result <- data.frame(beta1 =  beta1,
                        beta2 = beta2,
                        alpha1 = alpha1,
                        alpha2 = alpha2)

colnames(df_result) <- c("beta1", "beta2", "alpha1", "alpha2")

# CHOOSE EXTREME EPISODE FOR R-PARETO ##########################################

q <- 0.97 # quantile

# get central site from sites_coords
comephore_subset <- comephore[rownames(comephore) >= "2008-01-01", ]
set_st_excess <- get_spatiotemp_excess(comephore_subset, quantile = q,
                                      remove_zeros = TRUE)
# remove date column if exists
if ("date" %in% colnames(comephore_subset)) {
  comephore_subset <- comephore_subset[, -c(which(colnames(comephore_subset) == "date"))]
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
min_spatial_dist <- 5 # in km
delta <- 30 # in hours
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
list_results <- mclapply(1:length(s0_list), function(i) {
  s0 <- s0_list[i]
  s0_coords <- df_coords[s0, ]
  u <- u_list[i]
  episode <- list_episodes[[i]]
  ind_t0_ep <- 0 # index of t0 in the episode
  lags <- get_conditional_lag_vectors(df_coords, s0_coords, ind_t0_ep,
                                tau_vect, latlon = FALSE)
  excesses <- empirical_excesses_rpar(episode,
                                  threshold = u, # !!!!!!!
                                  df_lags = lags, t0 = ind_t0_ep)
  lags$tau <- lags$tau # convert tau from hours to seconds
  list(lags = lags, excesses = excesses)
}, mc.cores = detectCores() - 1)

list_lags <- lapply(list_results, `[[`, "lags")
list_excesses <- lapply(list_results, `[[`, "excesses")

# ADD ADVECTION ESTIMATES ######################################################
adv_filename <- paste0(data_folder, "comephore/adv_estim/oldies/advection_results_q",
                       q * 100, "_delta", delta, "_dmin", min_spatial_dist,
                       ".csv")
adv_df <- read.csv(adv_filename, sep = ",")

# get only matching episodes from selected_points
matching_indices <- match(selected_points$t0_date, adv_df$t0)
# remove NA indices
matching_indices <- matching_indices[!is.na(matching_indices)]
# get only matching rows
adv_df <- adv_df[matching_indices, ]
rownames(adv_df) <- NULL  # reset row names

selected_episodes <- selected_points
selected_episodes$adv_x <- rep(NA, nrow(selected_episodes))
selected_episodes$adv_y <- rep(NA, nrow(selected_episodes))

# get adv values for each episode according to the t0_date
for (i in 1:nrow(selected_episodes)) {
  t0_date <- selected_episodes$t0_date[i]
  adv_row <- adv_df[adv_df$t0 == t0_date, ]
  if (nrow(adv_row) == 0) {
    print(i)
    print(paste("No advection data found for t0_date =", t0_date))
  } else {
    # if there are multiple rows, take the first one
    adv_row <- adv_row[1, ]
    adv_x <- adv_row$mean_dx_mps
    adv_y <- adv_row$mean_dy_mps
    selected_episodes$adv_x[i] <- adv_x
    selected_episodes$adv_y[i] <- adv_y
  }
}

head(selected_episodes)
nrow(selected_episodes)
# Remove episodes with NA advection values
ind_NA_adv <- which(is.na(selected_episodes$adv_x) | is.na(selected_episodes$adv_y))
selected_episodes_nona <- selected_episodes[-ind_NA_adv, ]

wind_df <- selected_episodes[, c("adv_x", "adv_y")]
colnames(wind_df) <- c("vx", "vy")
length(wind_df$vx) # should be the same as number of episodes

# convert m/s to km/h
wind_df$vx <- selected_episodes$adv_x * 3.6  # convert to km/h
wind_df$vy <- selected_episodes$adv_y * 3.6  # convert to km/h
# OPTIMIZATION #################################################################

ind_NA <- which(is.na(wind_df$vx))
ind_NA <- ind_NA_adv
wind_opt <- wind_df
if (any(ind_NA > 0)) {
  # remove these episodes
  wind_opt <- wind_df[-ind_NA, ]
  episodes_opt <- list_episodes[-ind_NA]
  lags_opt <- list_lags[-ind_NA]
  excesses_opt <- list_excesses[-ind_NA]
} else {
  wind_opt <- wind_df
  episodes_opt <- list_episodes
  lags_opt <- list_lags
  excesses_opt <- list_excesses
}

# essayer de mettre en metres maybe
# changer init à 0.5

# fixer les etas à 1


init_param <- c(beta1, beta2, alpha1, alpha2, 1, 1)
result <- optim(par = init_param, fn = neg_ll_composite,
        list_lags = lags_opt, list_episodes = episodes_opt,
        list_excesses = excesses_opt, hmax = 7,
        wind_df = wind_opt,
        latlon = FALSE,
        distance = "euclidean",
        method = "L-BFGS-B",
        lower = c(1e-08, 1e-08, 1e-08, 1e-08, 1e-08, 1e-08),
        upper = c(10, 10, 1.999, 1.999, 2, 2),
        control = list(maxit = 10000,
                      trace = 1),
        hessian = FALSE)

# CIs using Jackknife blocks
n_total <- length(episodes_opt)
block_size <- 15
n_blocks <- floor(n_total / block_size)
blocks <- split(1:(block_size * n_blocks), rep(1:n_blocks, each = block_size))

n_cores <- detectCores() - 1

jack_estimates_list <- mclapply(1:n_blocks, function(i) {
  cat("Bloc", i, "over", n_blocks, "\n")
  
  exclude_idx <- blocks[[i]]
  
  jack_episodes <- episodes_opt[-exclude_idx]
  jack_lags <- lags_opt[-exclude_idx]
  jack_excesses <- excesses_opt[-exclude_idx]
  jack_wind <- wind_opt[-exclude_idx, , drop = FALSE]
  
  res <- tryCatch({
    optim(par = init_param, fn = neg_ll_composite,
          list_lags = jack_lags, list_episodes = jack_episodes,
          list_excesses = jack_excesses, hmax = 7,
          wind_df = jack_wind,
          latlon = FALSE,
          directional = TRUE,
          fixed_eta1 = FALSE,
          fixed_eta2 = FALSE,
          method = "L-BFGS-B",
          lower = c(1e-08, 1e-08, 1e-08, 1e-08, 1e-08, 1e-08),
          upper = c(10, 10, 1.999, 1.999, 10, 10),
          control = list(maxit = 10000),
          hessian = FALSE)
  }, error = function(e) NULL)
  
  if (!is.null(res)) {
    return(res$par)
  } else {
    return(rep(NA, length(init_param)))
  }
}, mc.cores = n_cores)

jack_estimates <- do.call(rbind, jack_estimates_list)

jack_estimates <- na.omit(jack_estimates)
n_eff <- nrow(jack_estimates)

filename <- paste0(data_folder, "comephore/optim_results_all_estimates.csv")
write.csv(jack_estimates, filename, row.names = FALSE)

jack_mean <- colMeans(jack_estimates)
jack_var <- (n_eff - 1) / n_eff * colSums((jack_estimates - matrix(rep(jack_mean, each = n_eff), ncol=6))^2)
jack_se <- sqrt(jack_var)

z <- qnorm(0.975)
lower_ci <- result$par - z * jack_se
upper_ci <- result$par + z * jack_se

jackknife_block_results <- data.frame(
  Parameter = c("beta1", "beta2", "alpha1", "alpha2", "eta1", "eta2"),
  Estimate = result$par,
  StdError = jack_se,
  CI_lower = lower_ci,
  CI_upper = upper_ci
)

filename <- paste0(data_folder, "comephore/optim_results_q", q * 100, "_delta", delta, "_dmin", min_spatial_dist, ".csv")
write.csv(jackknife_block_results, filename, row.names = FALSE)


