# remove objects from the environment
rm(list = ls())

# clear console
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
library(parallel)

# LOAD DATA ####################################################################
filename_com <- paste0(data_folder, "comephore/rebuild_clean/comephore_2008_2025_within5km.csv")
comephore_raw <- read.csv(filename_com, sep = ",")
filename_loc <- paste0(data_folder, "comephore/rebuild_clean/coords_pixels_within5km.csv")
loc_px <- read.csv(filename_loc, sep = ",")


loc_px <- loc_px[loc_px$pixel_name %in% colnames(comephore_raw)[-1], ]
rownames(loc_px) <- NULL

df_comephore <- as.data.frame(comephore_raw)
df_comephore$date <- as.POSIXct(df_comephore$date,
                                format = "%Y-%m-%d %H:%M:%S", tz = "UTC")
rownames(df_comephore) <- format(df_comephore$date, "%Y-%m-%d %H:%M:%S")
comephore <- df_comephore[-1]

# DISTANCE AND COORDS ##########################################################
nsites <- nrow(loc_px)
sites_coords <- loc_px[, c("Longitude", "Latitude")]
rownames(sites_coords) <- loc_px$pixel_name
sites_coords_sf <- st_as_sf(sites_coords, coords = c("Longitude", "Latitude"),
                            crs = 2154)
# sites_coords_sf <- st_transform(sites_coords_sf, crs = 2154)
coords_m <- st_coordinates(sites_coords_sf)
grid_coords_km <- as.data.frame(coords_m / 1000)
colnames(grid_coords_km) <- c("Longitude", "Latitude")
rownames(grid_coords_km) <- rownames(sites_coords)

head(comephore)
# remove column "dates"
comephore <- comephore[, which(colnames(comephore) != "date")]
# CHOOSE EXTREME EPISODE FOR R-PARETO ##########################################
set_st_excess <- get_spatiotemp_excess(data = comephore, quantile = q,
                                      remove_zeros = TRUE)

list_s <- set_st_excess$list_s
list_t <- set_st_excess$list_t
list_u <- set_st_excess$list_u

for (i in seq_along(list_s)) {
  s0 <- list_s[[i]]
  t0 <- list_t[[i]][1]
  u_s0 <- list_u[[i]][1]
  if (comephore[t0, s0] <= u_s0) {
    stop("Excess is not above threshold for s =", s0, " and t =", t0)
  }
}




library(data.table)
library(dplyr)
library(ggplot2)

dmins <- c(3, 4, 5, 6, 7, 8, 9, 10)  # in km

episode_pair_stats <- function(sel, sites_coords, dmin_km,
                 delta_steps,
                 latlon = FALSE, beta = 0) {

  if (nrow(sel) < 2) {
  return(list(
    n_episodes = nrow(sel),
    n_pairs = 0,
    n_pairs_dtlt_delta = 0,
    p_close_space_50km = NA_real_,
    p_close_space_120km = NA_real_,
    p_close_space_50km_given_dtlt_delta = NA_real_,
    p_conflict_dmin_delta = NA_real_,
    spatial_d10 = NA_real_,
    temporal_d10 = NA_real_
  ))
  }

  dist_matrix <- get_dist_mat(sites_coords, latlon = latlon)

  s <- sel$s0
  t <- sel$t0

  idx <- combn(seq_along(s), 2)
  s1 <- s[idx[1, ]]; s2 <- s[idx[2, ]]
  t1 <- t[idx[1, ]]; t2 <- t[idx[2, ]]

  spatial_d <- dist_matrix[cbind(s1, s2)]

  temporal_steps <- abs(t1 - t2)
  temporal_hours <- temporal_steps

  delta_eff_steps <- delta_steps + 2 * beta
  is_dt_close <- temporal_steps < delta_eff_steps

  p_close_50_all <- mean(spatial_d < 50, na.rm = TRUE)
  p_close_120_all <- mean(spatial_d < 120, na.rm = TRUE)

  n_pairs_dt <- sum(is_dt_close, na.rm = TRUE)
  p_close_50_given_dt <- if (n_pairs_dt == 0) NA_real_ else
  mean((spatial_d < 50)[is_dt_close], na.rm = TRUE)

  p_conflict <- mean((spatial_d < dmin_km) & is_dt_close, na.rm = TRUE)

  list(
  n_episodes = nrow(sel),
  n_pairs = length(spatial_d),
  n_pairs_dtlt_delta = n_pairs_dt,
  p_close_space_50km = p_close_50_all,
  p_close_space_120km = p_close_120_all,
  p_close_space_50km_given_dtlt_delta = p_close_50_given_dt,
  p_conflict_dmin_delta = p_conflict,
  spatial_d10 = as.numeric(quantile(spatial_d, 0.10, na.rm = TRUE)),
  temporal_d10 = as.numeric(quantile(temporal_hours, 0.10, na.rm = TRUE))
  )
}

episode_size <- 24
set_st_excess <- get_spatiotemp_excess(comephore, quantile = q, 
                     remove_zeros = TRUE)
# res <- lapply(dmins, function(dm) {
#   sel <- get_s0t0_pairs(
#   sites_coords = grid_coords_km,
#   data = comephore,
#   min_spatial_dist = dm,
#   episode_size = episode_size,
#   set_st_excess = set_st_excess,
#   n_max_episodes = 10000,
#   latlon = FALSE,
#   beta = 0
#   )

#   st <- episode_pair_stats(
#   sel = sel,
#   sites_coords = grid_coords_km,
#   dmin_km = dm,
#   delta_steps = episode_size,
#   latlon = FALSE,
#   beta = 0
#   )

#   data.frame(
#   dmin_km = dm,
#   n_episodes = st$n_episodes,
#   n_pairs_dtlt_delta = st$n_pairs_dtlt_delta,
#   p_close_space_50km_all = st$p_close_space_50km,
#   p_close_space_50km_given_dtlt_delta = st$p_close_space_50km_given_dtlt_delta,
#   p_conflict_dmin_delta = st$p_conflict_dmin_delta
#   )
# })


res <- lapply(dmins, function(dm) {
  sel <- get_s0t0_pairs(
    sites_coords = grid_coords_km,
    data = comephore,
    min_spatial_dist = dm,
    episode_size = episode_size,
    set_st_excess = set_st_excess,
    n_max_episodes = 10000,
    latlon = FALSE,
    beta = 0
  )

  data.frame(
    dmin_km = dm,
    n_episodes = nrow(sel)
  )
})

df_tradeoff <- bind_rows(res)

pA <- ggplot(df_tradeoff, aes(x = dmin_km, y = n_episodes)) +
  geom_line(size = 1.1, color = btfgreen) +
  geom_point(size = 2, color = btfgreen) +
  geom_vline(xintercept = 5, linetype = "dashed", color = "red") +
  theme_minimal() +
  labs(
  x = expression(d[min]~"(km)"),
  y = "Number of selected episodes"
  ) + 
  btf_theme

foldername <- paste0(im_folder, "/optim/comephore/choice_config/")
if (!dir.exists(foldername)) {
  dir.create(foldername, recursive = TRUE)
}
filename <- paste0(foldername, "tradeoff_dmin_episodes.png")
ggsave(filename, plot = pA, width = 7, height = 5, units = "in", dpi = 300)

dmin_fixed <- 5
delta_grid <- c(12, 15, 20, 24, 30, 36, 38, 40, 45, 48)  # in hours

set_st_excess <- get_spatiotemp_excess(comephore, quantile = q,
                     remove_zeros = TRUE)

res <- lapply(delta_grid, function(delta_steps) {
  sel <- get_s0t0_pairs(
    sites_coords = grid_coords_km,
    data = comephore,
    min_spatial_dist = dmin_fixed,
    episode_size = delta_steps,
    set_st_excess = set_st_excess,
    n_max_episodes = 10000,
    latlon = FALSE,
    beta = 0
  )

  data.frame(
    delta_steps = delta_steps,
    n_episodes = nrow(sel)
  )
})


df_tradeoff_delta <- bind_rows(res)


p_delta <- ggplot(df_tradeoff_delta, aes(x = delta_steps, y = n_episodes)) +
  geom_line(size = 1.1, color = btfgreen) +
  geom_point(size = 2, color = btfgreen) +
  geom_vline(xintercept = 24, linetype = "dashed", color = "red") +
  theme_minimal() +
  labs(
  x = expression(delta~" (hours)"),
  y = "Number of selected episodes"
  ) +
  btf_theme

print(p_delta)

ggsave(paste0(foldername, "tradeoff_delta_episodes_dmin5.png"),
     plot = p_delta, width = 7, height = 5, units = "in", dpi = 300)




