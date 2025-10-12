# remove objects from the environment
rm(list = ls())

# clear console
cat("\014")

# Load libraries and set theme
source("./script/load_libraries.R")

# Get all files in the folder "R"
functions_folder <- "./R"
files <- list.files(functions_folder, full.names = TRUE)
# # load all functions in files
invisible(lapply(files, function(f) source(f, echo = FALSE)))
library(latex2exp)


# get rain data from omsev
filename_omsev <- paste0(data_folder,
                         "omsev/omsev_5min/rain_mtp_5min_2019_2024.RData")

load(filename_omsev)

# get location of each rain gauge
filename_loc <- paste0(data_folder, "omsev/loc_rain_gauges.csv")
location_gauges <- read.csv(filename_loc)
location_gauges$Station <- c("iem", "mse", "poly", "um", "cefe", "cnrs",
                             "crbm", "archiw", "archie", "um35", "chu1",
                             "chu2", "chu3", "chu4", "chu5", "chu6", "chu7",
                             "cines", "brives", "hydro")

# save rain as csv
filename_rain <- paste0(data_folder, "omsev/omsev_5min/rain_mtp_5min_2019_2024_cleaned.csv")
# write.csv(rain, file = filename_rain, row.names = FALSE)

rain_omsev <- read.csv(filename_rain)
head(rain_omsev)
# begin in 2020-01-01
# rain_omsev <- rain_omsev[rain_omsev$dates >= "2020-01-01", ]

# put dates as rownames
rownames(rain_omsev) <- rain_omsev$dates
rain <- rain_omsev[-1] # remove dates column
# rain$mse
# Get distances matrix
dist_mat <- get_dist_mat(location_gauges)
df_dist <- reshape_distances(dist_mat)

# remove cines, hydro, brives
# colnames(rain)
# rain <- rain[, !(colnames(rain) %in% c("cines", "hydro", "brives"))]

# location_gauges <- location_gauges[location_gauges$Station != "cines" &
#                                    location_gauges$Station != "hydro" &
#                                    location_gauges$Station != "brives", ]
dist_mat <- get_dist_mat(location_gauges)
df_dist <- reshape_distances(dist_mat)
site_names <- location_gauges$Station
sites_names <- colnames(rain)

################################################################################
# WLSE results -----------------------------------------------------------------
################################################################################
foldername <- paste0(data_folder, "omsev/WLSE/")
df_spa_result <- read.csv(paste0(foldername, "wlse_spa_summary.csv"))
df_temp_result <- read.csv(paste0(foldername, "wlse_temp_summary.csv"))
df_result_all <- merge(df_spa_result, df_temp_result, by = NULL)

# select one row
df_result <- df_result_all[df_result_all$q_spa == q &
                             df_result_all$q_temp == q, ]

beta1 <- df_result$beta1
beta2 <- df_result$beta2
alpha1 <- df_result$alpha1
alpha2 <- df_result$alpha2

# estimates
wlse_omsev <- c(beta1, beta2, alpha1, alpha2) # m / 5 min

################################################################################
# VARIOGRAM --------------------------------------------------------------------
################################################################################

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

# Define parameter grids
q_values <- c(0.9, 0.92, 0.95, 0.97, 0.98)
delta_values <- c(7, 12, 15, 18, 24)
min_spatial_dist_values <- c(500, 750, 1000, 1200)

results_list <- list()
counter <- 1
tau_vect <- 0:10
for (q in q_values) {
    for (delta in delta_values) {
        for (min_spatial_dist in min_spatial_dist_values) {
            cat("Running for q =", q, "delta =", delta, "min_spatial_dist =", min_spatial_dist, "\n")
            
            set_st_excess <- get_spatiotemp_excess(rain, quantile = q, remove_zeros = TRUE)
            list_s <- set_st_excess$list_s
            list_t <- set_st_excess$list_t
            list_u <- set_st_excess$list_u
            
            s0t0_set <- get_s0t0_pairs(grid_coords_m, rain,
                                       min_spatial_dist = min_spatial_dist,
                                        episode_size = delta,
                                        set_st_excess = set_st_excess,
                                        n_max_episodes = 10000,
                                         latlon = FALSE)
            
            selected_points <- s0t0_set %>%
                mutate(t0_date = as.POSIXct(t0_date, format = "%Y-%m-%d %H:%M:%S", tz = "UTC"))
            
            n_episodes <- length(selected_points$s0)
            t0_list <- selected_points$t0
            s0_list <- selected_points$s0
            
            list_episodes_points <- get_extreme_episodes(selected_points, rain,
                                                         episode_size = episode_size, unif = FALSE,
                                                         beta = 0)
            
            list_episodes <- list_episodes_points$episodes
            selected_points$t0_date <- as.POSIXct(selected_points$t0_date, format = "%Y-%m-%d %H:%M:%S", tz = "UTC")
            selected_points$t0_date_rounded <- ceiling_date(selected_points$t0_date, "hour")
            
            datetimes <- unique(selected_points$t0_date)
            datetimes_hour <- unique(selected_points$t0_date_rounded)
            
            # Save datetime lists to csv
            datetime_filename <- paste(data_folder, "/omsev/t0_episodes_q", q * 100,
                                                                 "_delta", delta, "_dmin", min_spatial_dist,
                                                                 ".csv", sep = "")
            write.csv(data.frame(t0_date = datetimes_hour), datetime_filename, row.names = FALSE)
            datetime_filename <- paste(data_folder, "/omsev/t0_5min_episodes_q", q * 100,
                                                                 "_delta", delta, "_dmin", min_spatial_dist,
                                                                 ".csv", sep = "")
            write.csv(data.frame(t0_date = datetimes), datetime_filename, row.names = FALSE)
            
            selected_episodes <- selected_points

            V_episodes <- data.frame(
                v_x = rep(0, n_episodes),
                v_y = rep(0, n_episodes)
            )
            colnames(V_episodes) <- c("vx", "vy")
            
            # If selected_episodes is defined elsewhere, you may need to update this part
            s0_list <- selected_episodes$s0
            list_episodes_points <- get_extreme_episodes(selected_episodes, rain,
                                                episode_size = delta, unif = FALSE)
            list_episodes <- list_episodes_points$episodes
            u0_list <- selected_episodes$u_s0
            tmax <- max(tau_vect)
            df_coords <- as.data.frame(grid_coords_m)
            
            list_results <- mclapply(1:length(s0_list), function(i) {
                s0 <- s0_list[i]
                col_s0 <- which(colnames(rain) == s0)
                s0_coords <- df_coords[col_s0, ]
                episode <- list_episodes[[i]]
                ind_t0_ep <- 0
                lags <- get_conditional_lag_vectors(df_coords, s0_coords, ind_t0_ep,
                                                                                        tau_vect, latlon = FALSE)
                lags$hnorm <- lags$hnorm
                u <- u0_list[i]
                excesses <- empirical_excesses_rpar(episode, quantile = u, threshold = TRUE,
                                                                                        df_lags = lags, t0 = ind_t0_ep)
                lags$tau <- lags$tau
                list(lags = lags, excesses = excesses)
            }, mc.cores = detectCores() - 1)
            
            list_lags <- lapply(list_results, `[[`, "lags")
            list_excesses <- lapply(list_results, `[[`, "excesses")
            
            hmax <- max(dist_mat)
            init_params_omsev <- wlse_omsev
            
            result <- optim(par = init_params_omsev, fn = neg_ll_composite,
                                            list_lags = list_lags, list_episodes = list_episodes,
                                            list_excesses = list_excesses, hmax = hmax,
                                            wind_df = V_episodes,
                                            latlon = FALSE,
                                            distance = "lalpha",
                                            method = "L-BFGS-B",
                                            lower = c(1e-08, 1e-08, 1e-08, 1e-08),
                                            upper = c(10, 10, 1.999, 1.999),
                                            control = list(maxit = 10000, trace = 1),
                                            hessian = F)
            
            results_list[[counter]] <- list(
                q = q,
                delta = delta,
                min_spatial_dist = min_spatial_dist,
                optim_result = result
            )
            counter <- counter + 1
        }
    }
}
