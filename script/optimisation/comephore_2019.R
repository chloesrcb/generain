###############################################################################
# INITIALISATION
###############################################################################
rm(list = ls())
cat("\014")

# Load libraries
source("./script/load_libraries.R")
library(sf)
library(dplyr)
library(ggplot2)
library(gridExtra)
library(latex2exp)
library(lubridate)
library(knitr)
library(kableExtra)
library(parallel)
library(reshape2)

###############################################################################
# LOAD FUNCTIONS
###############################################################################
functions_folder <- "./R"
files <- list.files(functions_folder, full.names = TRUE)
invisible(lapply(files, function(f) source(f, echo = FALSE)))

###############################################################################
# LOAD DATA
###############################################################################
filename_com <- paste0(data_folder, "comephore/comephore_2008_2024_5km.csv")
comephore_raw <- read.csv(filename_com, sep = ",")

filename_loc <- paste0(data_folder, "comephore/loc_px_zoom_5km.csv")
loc_px <- read.csv(filename_loc, sep = ",")

# Filtrage des pixels présents dans les données
loc_px <- loc_px[loc_px$pixel_name %in% colnames(comephore_raw)[-1], ]
rownames(loc_px) <- NULL

# Conversion date
df_comephore <- as.data.frame(comephore_raw)
df_comephore$date <- as.POSIXct(df_comephore$date,
                                format = "%Y-%m-%d %H:%M:%S", tz = "UTC")

# Filtrage temporel
df_comephore <- df_comephore[
  df_comephore$date >= as.POSIXct("2008-01-01 00:00:00", tz = "UTC"),
]

# Index
rownames(df_comephore) <- format(df_comephore$date, "%Y-%m-%d %H:%M:%S")
comephore <- df_comephore[-1]

###############################################################################
# GESTION DES COORDONNÉES
###############################################################################
nsites <- nrow(loc_px)
sites_coords <- loc_px[, c("Longitude", "Latitude")]
rownames(sites_coords) <- loc_px$pixel_name
sites_coords_sf <- st_as_sf(sites_coords, coords = c("Longitude", "Latitude"), crs = 4326)
sites_coords_sf <- st_transform(sites_coords_sf, crs = 2154)

coords_m <- st_coordinates(sites_coords_sf)
grid_coords_km <- sites_coords
grid_coords_m  <- sites_coords
grid_coords_m$x_m  <- coords_m[, "X"] - min(coords_m[, "X"])
grid_coords_m$y_m  <- coords_m[, "Y"] - min(coords_m[, "Y"])
grid_coords_km$x_km <- grid_coords_m$x_m / 1000
grid_coords_km$y_km <- grid_coords_m$y_m / 1000

filename <- paste0(data_folder, "comephore/grid/grid_coords_5km.csv")
write.csv(grid_coords_km, file = filename, row.names = TRUE)

grid_coords_m <- grid_coords_m[, c("x_m", "y_m")]
colnames(grid_coords_m) <- c("Longitude", "Latitude")
dist_mat <- get_dist_mat(grid_coords_m, latlon = FALSE)

df_dist <- reshape_distances(dist_mat)
df_dist_km <- df_dist
df_dist_km$value <- round(df_dist$value / 1000, 1)

###############################################################################
# CHARGEMENT RÉSULTATS WLSE
###############################################################################
hmax <- 7
tmax <- 10
foldername <- paste0(data_folder, "/comephore/WLSE/")

df_spa_result <- read.csv(paste0(foldername, "wlse_spa_summary.csv"))
df_temp_result <- read.csv(paste0(foldername, "wlse_temp_summary.csv"))
df_result_all <- merge(df_spa_result, df_temp_result, by = NULL)
df_result_all <- df_result_all[, c("q_spa", "q_temp", "beta1", "alpha1", "beta2", "alpha2")]

df_result_all <- df_result_all %>%
  mutate(across(c(beta1, alpha1, beta2, alpha2), ~round(as.numeric(.), 4)))

df_result <- df_result_all[df_result_all$q_spa == 0.95 & df_result_all$q_temp == 0.95, ]

###############################################################################
# CONFIGURATION GLOBALE (années / suffixes fichiers)
###############################################################################
first_ts <- as.POSIXct(rownames(comephore)[1], tz = "UTC")
starting_year <- year(first_ts)
starting_year_suffix <- if (starting_year == 2008) "" else paste0("_from", starting_year)

###############################################################################
# SÉLECTION DES ÉPISODES
###############################################################################
configs <- expand.grid(
  q = c(0.95, 0.97),
  min_spatial_dist = c(5, 7),
  delta = c(30)
)

episode_counts <- data.frame()

for (i in seq_len(nrow(configs))) {
  q <- configs$q[i]
  min_spatial_dist <- configs$min_spatial_dist[i]
  delta <- configs$delta[i]

  comephore_subset <- comephore[rownames(comephore) >= "2008-01-01", ]

  set_st_excess <- get_spatiotemp_excess(comephore_subset, quantile = q, remove_zeros = TRUE)
  list_s <- set_st_excess$list_s
  list_t <- set_st_excess$list_t
  list_u <- set_st_excess$list_u

  s0t0_set <- get_s0t0_pairs(grid_coords_km, comephore_subset,
                             min_spatial_dist = min_spatial_dist,
                             episode_size = delta,
                             set_st_excess = set_st_excess,
                             n_max_episodes = 10000,
                             latlon = FALSE)

  selected_points <- s0t0_set

  filename <- paste0(data_folder, "comephore/episodes/s0t0_pairs/s0t0_pairs_q",
                     q * 100, "_delta", delta, "_dmin", min_spatial_dist, starting_year_suffix, ".csv")
  write.csv(selected_points, file = filename, row.names = FALSE)

  datetimes <- as.POSIXct(selected_points$t0_date, tz = "UTC")
  datetimes_df <- data.frame(datetime = unique(datetimes))
  filename_dt <- paste0(data_folder, "comephore/episodes/t0_episodes/t0_episodes",
                        "_q", q * 100, "_delta", delta, "_dmin", min_spatial_dist,
                        starting_year_suffix, ".csv")
  write.csv(datetimes_df, file = filename_dt, row.names = FALSE)

  u_min <- min(selected_points$u_s0)
  u_max <- max(selected_points$u_s0)

  episode_counts <- rbind(episode_counts, data.frame(
    q = q, min_spatial_dist = min_spatial_dist, delta = delta,
    n_episodes = nrow(selected_points), u_min = u_min, u_max = u_max
  ))
}

foldername <- paste0(data_folder, "comephore/episodes/")
write.csv(episode_counts,
          paste0(foldername, "episodes_selection_summary_zoom5km", starting_year_suffix, ".csv"),
          row.names = FALSE)

###############################################################################
# ANALYSE D’UN ÉPISODE (RECHARGEMENT)
###############################################################################
q <- 0.95
min_spatial_dist <- 5
delta <- 30

# Rechargement cohérent
filename <- paste0(data_folder, "comephore/episodes/t0_episodes/t0_episodes_q",
                   q * 100, "_delta", delta, "_dmin", min_spatial_dist, starting_year_suffix, ".csv")
datetimes_df <- read.csv(filename, sep = ",")
datetimes_df$datetime <- ifelse(
  nchar(datetimes_df$datetime) == 10,
  paste0(datetimes_df$datetime, " 00:00:00"),
  datetimes_df$datetime
)
datetimes_df$datetime <- as.POSIXct(datetimes_df$datetime, tz = "UTC")

# Recharger les paires
file_s0t0 <- paste0(data_folder, "comephore/episodes/s0t0_pairs/s0t0_pairs_q",
                    q * 100, "_delta", delta, "_dmin", min_spatial_dist, starting_year_suffix, ".csv")
selected_points <- read.csv(file_s0t0)
selected_points$t0_date <- as.POSIXct(selected_points$t0_date, tz = "UTC")

list_episodes_points <- get_extreme_episodes(selected_points, comephore,
                                     episode_size = delta, unif = FALSE,
                                     beta = 0)

list_episodes <- list_episodes_points$episodes
episode <- list_episodes[[1]]
s0_list <- selected_points$s0
u_list <- selected_points$u_s0
tau_vect <- 0:10
tmax <- max(tau_vect)
df_coords <- as.data.frame(grid_coords_km)
# Compute the lags and excesses for each conditional point
list_results <- mclapply(1:length(s0_list), function(i) {
  s0 <- s0_list[i]
  s0_coords <- df_coords[s0, ]
  # t0 <- t0_list[i]
  u <- u_list[i]
  episode <- list_episodes[[i]]
  ind_t0_ep <- 0 # index of t0 in the episode
  lags <- get_conditional_lag_vectors(df_coords, s0_coords, ind_t0_ep,
                                tau_vect, latlon = FALSE)
  # lags$hnorm <- lags$hnorm / 1000 # convert to km
  excesses <- empirical_excesses_rpar(episode,
                                  threshold = u, # !!!!!!!
                                  df_lags = lags, t0 = ind_t0_ep)
  lags$tau <- lags$tau # convert tau from hours to seconds
  list(lags = lags, excesses = excesses)
}, mc.cores = detectCores() - 1)

list_lags <- lapply(list_results, `[[`, "lags")
list_excesses <- lapply(list_results, `[[`, "excesses")


###############################################################################
# AJOUT DES ESTIMATIONS D’ADVECTION
###############################################################################
adv_filename <- paste0(data_folder, "comephore/adv_estim/advection_results_q",
                       q * 100, "_delta", delta, "_dmin", min_spatial_dist, starting_year_suffix, ".csv")
adv_df <- read.csv(adv_filename, sep = ",")
adv_df$t0 <- as.POSIXct(adv_df$t0, tz = "UTC")

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
    selected_episodes$adv_x[i] <- adv_row$mean_dx_mps[1]
    selected_episodes$adv_y[i] <- adv_row$mean_dy_mps[1]
  }
}

# Nettoyage des NA
ind_NA_adv <- which(is.na(selected_episodes$adv_x) | is.na(selected_episodes$adv_y))
selected_episodes_nona <- selected_episodes[-ind_NA_adv, ]

wind_df <- data.frame(
  vx = selected_episodes$adv_x * 3.6,  # m/s -> km/h
  vy = selected_episodes$adv_y * 3.6
)

###############################################################################
# OPTIMISATION
###############################################################################
ind_NA <- ind_NA_adv
wind_opt <- if (length(ind_NA) > 0) wind_df[-ind_NA, ] else wind_df

beta1 <- df_result$beta1
beta2 <- df_result$beta2
alpha1 <- df_result$alpha1
alpha2 <- df_result$alpha2
init_param <- c(beta1, beta2, alpha1, alpha2, 1, 1)

result <- optim(
  par = init_param,
  fn = neg_ll_composite,
  list_lags = list_lags,
  list_episodes = list_episodes,
  list_excesses = list_excesses,
  hmax = 7,
  wind_df = wind_opt,
  latlon = FALSE,
  distance = "lalpha",
  method = "L-BFGS-B",
  lower = c(1e-08, 1e-08, 1e-08, 1e-08, 1e-08, 1e-08),
  upper = c(10, 10, 1.999, 1.999, 10, 10),
  control = list(maxit = 10000, trace = 1),
  hessian = TRUE
)



result <- optim(
  par = init_param,
  fn = neg_ll_composite,
  list_lags = list_lags,
  list_episodes = list_episodes,
  list_excesses = list_excesses,
  hmax = 7,
  wind_df = wind_opt,
  latlon = FALSE,
  distance = "euclidean",
  method = "L-BFGS-B",
  lower = c(1e-08, 1e-08, 1e-08, 1e-08, 1e-08, 1e-08),
  upper = c(10, 10, 1.999, 1.999, 10, 10),
  control = list(maxit = 10000, trace = 1),
  hessian = TRUE
)
