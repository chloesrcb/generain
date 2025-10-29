###############################################################################
# INITIALISATION
###############################################################################
rm(list = ls())
cat("\014")

muse <- FALSE

if (muse) {
  folder_muse <- "/home/serrec/work_rainstsimu/rpareto/"
  setwd(folder_muse)
  source("load_libraries.R")
  im_folder <- "./images"
  source("config_com.R")
  data_folder <- "./data/"
  ncores <- 27
} else {
  source("./script/load_libraries.R")
  source("./script/optimisation/config_com.R")
  ncores <- parallel::detectCores() - 1
}

###############################################################################
# CHARGEMENT DES FONCTIONS
###############################################################################
functions_folder <- "./R"
files <- list.files(functions_folder, full.names = TRUE)
invisible(lapply(files, function(f) source(f, echo = FALSE)))

library(latex2exp)
library(sf)
library(dplyr)
library(ggplot2)
library(gridExtra)
library(lubridate)
library(parallel)
library(numDeriv)

###############################################################################
# PARAMÈTRES
###############################################################################
eta_type <- if (!is.na(fixed_eta1) && !is.na(fixed_eta2)) {
  "fixed_eta"
} else if (!is.na(fixed_eta1)) {
  "fixed_eta1"
} else if (!is.na(fixed_eta2)) {
  "fixed_eta2"
} else {
  "free_eta"
}

distance_type <- "lalpha"

# Création du dossier résultat
foldername_res <- file.path(
  paste0(data_folder, "comephore/optim_results"),
  distance_type,
  eta_type
)
if (!dir.exists(foldername_res)) {
  dir.create(foldername_res, recursive = TRUE)
  message("Created folder: ", foldername_res)
}

###############################################################################
# CHARGEMENT DES DONNÉES
###############################################################################
filename_com <- paste0(data_folder, "comephore/comephore_2008_2024_5km.csv")
comephore_raw <- read.csv(filename_com, sep = ",")

filename_loc <- paste0(data_folder, "comephore/loc_px_zoom_5km.csv")
loc_px <- read.csv(filename_loc, sep = ",")

# Filtrage des pixels cohérents
loc_px <- loc_px[loc_px$pixel_name %in% colnames(comephore_raw)[-1], ]
rownames(loc_px) <- NULL

df_comephore <- as.data.frame(comephore_raw)
colnames(df_comephore)[1] <- "date"
df_comephore$date <- as.POSIXct(df_comephore$date, format = "%Y-%m-%d %H:%M:%S", tz = "UTC")

df_comephore <- df_comephore[
  df_comephore$date >= as.POSIXct("2008-01-01 00:00:00", tz = "UTC"),
]

rownames(df_comephore) <- format(df_comephore$date, "%Y-%m-%d %H:%M:%S")
comephore <- df_comephore[-1]
head(comephore)
###############################################################################
# COORDONNÉES ET DISTANCES
###############################################################################
sites_coords <- loc_px[, c("Longitude", "Latitude")]
rownames(sites_coords) <- loc_px$pixel_name

sites_coords_sf <- st_as_sf(sites_coords, coords = c("Longitude", "Latitude"), crs = 4326)
sites_coords_sf <- st_transform(sites_coords_sf, crs = 2154)

coords_m <- st_coordinates(sites_coords_sf)
grid_coords_km <- sites_coords
grid_coords_m <- sites_coords
grid_coords_m$x_m <- coords_m[, "X"] - min(coords_m[, "X"])
grid_coords_m$y_m <- coords_m[, "Y"] - min(coords_m[, "Y"])
grid_coords_km$x_km <- grid_coords_m$x_m / 1000
grid_coords_km$y_km <- grid_coords_m$y_m / 1000
# remove original lon/lat columns
grid_coords_km <- grid_coords_km[, c("x_km", "y_km")]
colnames(grid_coords_km) <- c("Longitude", "Latitude")
rownames(grid_coords_km) <- rownames(sites_coords)
head(grid_coords_km)

###############################################################################
# WLSE SPATIAL & TEMPOREL
###############################################################################
foldername <- paste0(data_folder, "/comephore/WLSE/")
df_spa_result <- read.csv(paste0(foldername, "wlse_spa_summary.csv"))
df_temp_result <- read.csv(paste0(foldername, "wlse_temp_summary.csv"))
df_result_all <- merge(df_spa_result, df_temp_result, by = NULL)
df_result_all <- df_result_all[, c("q_spa", "q_temp", "beta1", "alpha1", "beta2", "alpha2")]
df_result_all <- df_result_all %>%
  mutate(across(c(beta1, alpha1, beta2, alpha2), ~ round(as.numeric(.), 4)))

df_result <- df_result_all[df_result_all$q_spa == q & df_result_all$q_temp == q, ]
beta1 <- df_result$beta1
beta2 <- df_result$beta2
alpha1 <- df_result$alpha1
alpha2 <- df_result$alpha2

###############################################################################
# DÉTERMINATION DE L’ANNÉE DE DÉBUT (pour cohérence des noms de fichiers)
###############################################################################
first_ts <- as.POSIXct(rownames(comephore)[1], tz = "UTC")
starting_year <- year(first_ts)
starting_year_suffix <- if (starting_year == 2008) "" else paste0("_from", starting_year)

###############################################################################
# EXTRACTION DES ÉPISODES EXTRÊMES
###############################################################################
comephore_subset <- comephore[rownames(comephore) >= as.POSIXct("2008-01-01", tz = "UTC"), ]
set_st_excess <- get_spatiotemp_excess(data = comephore_subset, quantile = q, remove_zeros = TRUE)

list_s <- set_st_excess$list_s
list_t <- set_st_excess$list_t
list_u <- set_st_excess$list_u

for (i in seq_along(list_s)) {
  s0 <- list_s[[i]]
  t0 <- list_t[[i]][1]
  u_s0 <- list_u[[i]][1]
  if (comephore_subset[t0, s0] <= u_s0) {
    stop(paste("Excess is not above threshold for s =", s0, "and t =", t0))
  }
}

s0t0_set <- get_s0t0_pairs(grid_coords_km, comephore_subset,
                            min_spatial_dist = min_spatial_dist,
                            episode_size = delta,
                            set_st_excess = set_st_excess,
                            n_max_episodes = 10000,
                            latlon = FALSE)
selected_points <- s0t0_set
head(selected_points)
# > head(selected_points)
#        s0    t0             t0_date  u_s0
#    <char> <int>              <char> <num>
# 1:    p79  3022 2008-05-05 22:00:00 13.57
# 2:   p147  3022 2008-05-05 22:00:00 15.00
# 3:   p181  3022 2008-05-05 22:00:00 14.00
# 4:   p251  3022 2008-05-05 22:00:00 14.00
# 5:    p19  3276 2008-05-16 12:00:00 15.00
# 6:    p56  3291 2008-05-17 03:00:00 15.00
###############################################################################
# ÉPISODES ET LAGS
###############################################################################
list_episodes_points <- get_extreme_episodes(selected_points, comephore_subset,
                                             episode_size = delta, unif = FALSE, beta = 0)
list_episodes <- list_episodes_points$episodes
s0_list <- selected_points$s0
t0_list <- selected_points$t0
u_list <- selected_points$u_s0

tau_vect <- 0:10
df_coords <- as.data.frame(grid_coords_km)

list_results <- parallel::mclapply(1:length(s0_list), function(i) {
  s0 <- s0_list[i]
  s0_coords <- df_coords[s0, ]
  u <- u_list[i]
  episode <- list_episodes[[i]]
  ind_t0_ep <- 0
  lags <- get_conditional_lag_vectors(df_coords, s0_coords, ind_t0_ep, tau_vect, latlon = FALSE)
  excesses <- empirical_excesses_rpar(episode, threshold = u, df_lags = lags, t0 = ind_t0_ep)
  list(lags = lags, excesses = excesses)
}, mc.cores = ncores)

list_lags <- lapply(list_results, `[[`, "lags")
list_excesses <- lapply(list_results, `[[`, "excesses")
# i <- 10
# lags <- list_lags[[i]]
# excesses <- list_excesses[[i]]
# episodes <- list_episodes[[i]]
# u <- u_list[[i]]
# sum(episodes[2,]>u)  # should be equal to sum(excesses$kij)
# head(lags)
# head(excesses)
# sum(excesses$kij[lags$tau == 1])  # should be equal to nrow(episodes) - 1
###############################################################################
# AJOUT DES ADVECTIONS
###############################################################################
adv_filename <- paste0(data_folder, "comephore/adv_estim/advection_results_q",
                       q * 100, "_delta", delta, "_dmin", min_spatial_dist,
                       starting_year_suffix, ".csv")
adv_df <- read.csv(adv_filename, sep = ",")
adv_df$t0 <- as.POSIXct(adv_df$t0, tz = "UTC")
selected_points$t0_date <- as.POSIXct(selected_points$t0_date, tz = "UTC")

adv_merged <- merge(selected_points[, c("t0_date")],
                    adv_df, by.x = "t0_date", by.y = "t0",
                    all.x = TRUE, sort = FALSE)
selected_episodes <- selected_points
selected_episodes$adv_x <- adv_merged$mean_dx_kmh
selected_episodes$adv_y <- adv_merged$mean_dy_kmh

wind_df <- data.frame(vx = selected_episodes$adv_x,
                      vy = selected_episodes$adv_y)

stopifnot(length(list_episodes) == nrow(wind_df))
head(wind_df)
# > head(wind_df)
#            vx         vy
# 1  0.84658645  0.7716187
# 2  0.84658645  0.7716187
# 3  0.84658645  0.7716187
# 4  0.84658645  0.7716187
# 5 -0.03622333  0.3916635
# 6  1.74948968 -1.0378787
# > tail(wind_df)
#               vx          vy
# 1058  0.73110030 -0.39961678
# 1059  0.73110030 -0.39961678
# 1060 -0.98791199 -0.03985671
# 1061 -0.05205483  0.15432273
# 1062 -0.05205483  0.15432273
# 1063 -0.05205483  0.15432273
###############################################################################
# OPTIMISATION
###############################################################################
message("Starting optimization...")

if (eta_type == "free_eta") {
  init_param <- c(beta1, beta2, alpha1, alpha2, 1, 1)
  result <- optim(
    par = init_param, fn = neg_ll_composite,
    list_lags = list_lags, list_episodes = list_episodes,
    list_excesses = list_excesses, hmax = 7,
    wind_df = wind_df, latlon = FALSE,
    distance = distance_type, method = "L-BFGS-B",
    lower = rep(1e-08, 6),
    upper = c(10, 10, 1.999, 1.999, 10, 10),
    control = list(maxit = 10000, trace = 1),
    hessian = FALSE
  )
} else if (eta_type == "fixed_eta") {
  init_param <- c(beta1, beta2, alpha1, alpha2)
  result <- optim(
    par = init_param, fn = neg_ll_composite_fixed_eta,
    list_lags = list_lags, list_episodes = list_episodes,
    list_excesses = list_excesses, hmax = 7,
    wind_df = wind_df, latlon = FALSE,
    distance = distance_type,
    fixed_eta1 = fixed_eta1, fixed_eta2 = fixed_eta2,
    method = "L-BFGS-B",
    lower = rep(1e-08, 4),
    upper = c(10, 10, 1.999, 1.999),
    control = list(maxit = 10000, trace = 1),
    hessian = FALSE
  )
}
# 0.3 0.7, 0.4, 0.7, 1.4, 5.4

if (result$convergence != 0) stop("Optimization did not converge")
message("Optimization converged")

###############################################################################
# GRADIENT & SAUVEGARDE DES RÉSULTATS
###############################################################################
objfun <- function(p) {
  if (eta_type == "fixed_eta") {
    neg_ll_composite_fixed_eta(p, list_episodes, list_excesses, list_lags,
                               wind_df = wind_df, latlon = FALSE,
                               distance = distance_type,
                               fixed_eta1 = fixed_eta1, fixed_eta2 = fixed_eta2,
                               hmax = 7)
  } else {
    neg_ll_composite(p, list_episodes, list_excesses, list_lags,
                     wind_df = wind_df, latlon = FALSE,
                     distance = distance_type, hmax = 7)
  }
}

grad_vec <- grad(objfun, result$par)
print(grad_vec)

result_df <- data.frame(
  beta1 = result$par[1],
  beta2 = result$par[2],
  alpha1 = result$par[3],
  alpha2 = result$par[4],
  eta1 = ifelse(length(result$par) >= 5, result$par[5], fixed_eta1),
  eta2 = ifelse(length(result$par) >= 6, result$par[6], fixed_eta2),
  nll = result$value,
  grad_beta1 = grad_vec[1],
  grad_beta2 = grad_vec[2],
  grad_alpha1 = grad_vec[3],
  grad_alpha2 = grad_vec[4],
  grad_eta1 = ifelse(length(grad_vec) >= 5, grad_vec[5], NA),
  grad_eta2 = ifelse(length(grad_vec) >= 6, grad_vec[6], NA)
)

filename <- paste0(foldername_res, "/results_q",
                   q * 100, "_delta", delta, "_dmin", min_spatial_dist,
                   starting_year_suffix, ".csv")
write.csv(result_df, filename, row.names = FALSE)
cat("Results saved to", filename, "\n")
