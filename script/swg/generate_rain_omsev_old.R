
# remove objects from the environment
rm(list = ls())

# clear console
cat("\014")

# Load libraries and set theme
source("./script/load_libraries.R")

# library(generain)
library(ggplot2)
library(reshape2)
library(animation)
library(RandomFields)
library(RandomFieldsUtils)
library(dplyr)
library(tidyr)
library(lubridate)
library(purrr)


# Get all files in the folder "R"
functions_folder <- "./R"
files <- list.files(functions_folder, full.names = TRUE)
# load all functions in files
invisible(lapply(files, function(f) source(f, echo = FALSE)))
library(latex2exp)

# get rain data
filename_rain <- paste0(data_folder, "omsev/omsev_5min/rain_mtp_5min_2019_2024_cleaned.csv")
rain_omsev <- read.csv(filename_rain)
head(rain_omsev)

# put dates as rownames
rownames(rain_omsev) <- rain_omsev$dates
rain <- rain_omsev[-1] # remove dates column

# get location of each rain gauge
filename_loc <- paste0(data_folder, "omsev/loc_rain_gauges.csv")
location_gauges <- read.csv(filename_loc)
location_gauges$Station <- c("iem", "mse", "poly", "um", "cefe", "cnrs",
                             "crbm", "archiw", "archie", "um35", "chu1",
                             "chu2", "chu3", "chu4", "chu5", "chu6", "chu7",
                             "cines", "brives", "hydro")


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
grid_omsev <- as.data.frame(coords_m / 1000)
colnames(grid_omsev) <- c("Longitude", "Latitude")
rownames(grid_omsev) <- rownames(sites_coords)
# params
q <- 0.95
min_spatial_dist <- 1200 # in meters
delta <- 12  # in time steps (5 min each)

# in rain remove when all data are NA
set_st_excess <- get_spatiotemp_excess(rain, quantile = q, remove_zeros = TRUE)

# Spatio-temporal neighborhood parameters
s0t0_set <- get_s0t0_pairs(grid_omsev, rain,
                            min_spatial_dist = min_spatial_dist,
                            episode_size = delta,
                            set_st_excess = set_st_excess,
                            n_max_episodes = 10000,
                            latlon = FALSE)

selected_points <- s0t0_set
selected_points <- selected_points %>%
  mutate(t0_date = as.POSIXct(t0_date, format = "%Y-%m-%d %H:%M:%S", tz = "UTC"))

n_episodes <- length(selected_points$s0)
t0_list <- selected_points$t0
s0_list <- selected_points$s0

list_episodes_points <- get_extreme_episodes(selected_points, rain,
                                     episode_size = delta, unif = FALSE)
list_episodes <- list_episodes_points$episodes

# Round to the next hour
selected_points$t0_date_rounded <- ceiling_date(selected_points$t0_date, "hour")

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

params_vario <- list(
  beta1 = 0.1,
  beta2 = 0.9,
  alpha1 = 0.5,
  alpha2 = 0.6
)


# compute probability margins parameters
p0_list <- apply(rain, 2, function(x) {
  mean(x == 0, na.rm = TRUE)
})
params_margins <- list(
  xi    = rep(0.25, nrow(grid_omsev)),
  sigma = rep(0.6, nrow(grid_omsev)),
  kappa = rep(0.3, nrow(grid_omsev)),
  p0   = p0_list
)
episodes <- vector("list", n_episodes)
u_s <- apply(rain, 2, function(x) {
  quantile(x[x > 0], probs = q, na.rm = TRUE)
})

# simulate all episodes
times <- 0:(delta - 1)
set.seed(1423555)
for (i in seq_len(n_episodes)) {
  cat("Simulation épisode", i)
  s0 <- s0_list[[i]]
  coord_s0 <- grid_omsev[rownames(grid_omsev) == s0, ]
  coord_vec <- as.numeric(coord_s0)
  adv_vect <- c(V_episodes$vx[i], V_episodes$vy[i])
  episodes[[i]] <- sim_episode_coords(
    params_vario = params_vario,
    params_margins = params_margins,
    coords = grid_omsev,
    times = times,
    adv = adv_vect,
    t0 = 0,
    s0 = coord_vec,
    u_s = u_s
  )
  
}

episodes_df <- vector("list", n_episodes)

for (i in seq_len(n_episodes)) {
  ep <- episodes[[i]]           # ta matrice simulée [site × time]
  ep_df <- as.data.frame(t(ep)) # transpose pour avoir sites en colonnes
  colnames(ep_df) <- rownames(grid_omsev)
  s0 <- s0_list[[i]]
  # Récupère la date de départ de l’épisode réel
  t0_date <- selected_points$t0_date[i]
  
  # Crée la séquence temporelle (une date toutes les 5 minutes ici ? à ajuster selon ton pas)
  time_step <- difftime(list_episodes[[1]] %>% rownames() %>% .[2],
                        list_episodes[[1]] %>% rownames() %>% .[1],
                        units = "mins") %>% as.numeric()
  
  if (is.na(time_step)) time_step <- 5  # fallback par défaut si non défini
  
  times_seq <- t0_date + seq(0, delta - 1) * time_step * 60
  rownames(ep_df) <- format(times_seq, "%Y-%m-%d %H:%M:%S")
  
  episodes_df[[i]] <- ep_df
}

all_episodes_df <- do.call(rbind, episodes_df)
head(all_episodes_df)
starting_date <- min(rain_omsev$dates)
ending_date <- max(rain_omsev$dates)

full_dates <- seq(
  from = as.POSIXct(starting_date, tz = "UTC"),
  to   = as.POSIXct(ending_date, tz = "UTC"),
  by   = paste(time_step, "mins")
)


start_all <- as.POSIXct(rownames(rain)[1], tz = "UTC")
end_all   <- as.POSIXct("2025-01-31 23:59:00", tz = "UTC")
step_min  <- 5

full_dates <- seq(from = start_all, to = end_all, by = paste(step_min, "mins"))

episodes_df_dated <- map2(
  episodes_df, seq_along(episodes_df),
  ~{
    ep   <- .x
    i    <- .y
    t0   <- selected_points$t0_date[i]
    delta <- nrow(ep)
    times_seq <- seq(from = t0, by = paste(step_min,"mins"), length.out = delta)
    ep$time <- times_seq
    ep
  }
)

comb <- bind_rows(episodes_df_dated)
comb$time <- as.POSIXct(comb$time, tz = "UTC")
site_cols <- colnames(rain)

missing_cols <- setdiff(site_cols, names(comb))
if (length(missing_cols) > 0) {
  comb[missing_cols] <- 0
}

comb <- comb %>% select(time, all_of(site_cols))

comb_agg <- comb %>%
  group_by(time) %>%
  summarise(across(everything(), ~ sum(.x, na.rm = TRUE)), .groups = "drop")

final_df <- tibble(time = full_dates) %>%
  left_join(comb_agg, by = "time") %>%
  mutate(across(all_of(site_cols), ~ replace_na(.x, 0)))

stopifnot(identical(colnames(final_df)[-1], site_cols))
print(head(final_df, 10))


x11()

# plot real vs simulated for a site
site <- "poly"
real_rain_site <- rain[, site]
simu_rain_site <- final_df[[site]]
real_rain_df <- data.frame(
  time = as.POSIXct(rownames(rain), tz = "UTC"),
  rainfall = real_rain_site
)
simu_rain_df <- data.frame(
  time = final_df$time,
  rainfall = simu_rain_site
)


first_date <- min(real_rain_df$time, na.rm = TRUE)
first_date_sim <- min(final_df$time, na.rm = TRUE)
# remove 0 to plot 
real_rain_df <- real_rain_df %>%
  filter(rainfall > 0 & time >= first_date_sim & time <= max(final_df$time, na.rm = TRUE))
simu_rain_df <- simu_rain_df %>%
  filter(rainfall > 0 & time >= first_date & time <= max(real_rain_df$time, na.rm = TRUE))

ggplot() +
  geom_point(data = real_rain_df, aes(x = time, y = rainfall), color = "red", alpha = 0.5) +
  geom_point(data = simu_rain_df, aes(x = time, y = rainfall), color = "blue", alpha = 0.5) +
  labs(title = paste("Real vs Simulated Rainfall at Site", site),
       x = "Time",
       y = "Rainfall Intensity") +
  theme_minimal() +
  scale_color_manual(name = "Data Type", 
             values = c("Real" = "red", "Simulated" = "blue"),
             labels = c("Real", "Simulated")) +
    guides(color = guide_legend(override.aes = list(alpha = 1)))

episode_all <- final_df[, -1]

i <- 100
s0 <- s0_list[[i]]
coords_s0 <- location_gauges[location_gauges$Station == s0, ]
coords_s0 <- coords_s0[, c("Longitude", "Latitude")]
colnames(coords_s0) <- c("x", "y")
create_generator_gif_points(
  simu_df = episodes_df[[i]],
  coords = location_gauges,
  outfile = paste0("simulation_omsev_ep", i, ".gif"),
  interval = 0.5,
  forcedtemp = 30,
  s0 = coords_s0
)


s0 <- s0_list[[i]]
coords_s0 <- location_gauges[location_gauges$Station == s0, ]
coords_s0 <- coords_s0[, c("Longitude", "Latitude")]
colnames(coords_s0) <- c("x", "y")
create_generator_gif_points(
  simu_df = list_episodes[[i]],
  coords = location_gauges,
  outfile = paste0("real_omsev_ep", i, ".gif"),
  interval = 0.5,
  forcedtemp = 30,
  s0 = coords_s0
)
