
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

library(sf)
library(dplyr)
library(ggplot2)
library(gridExtra)

# load("workspace.RData")
# library(generain)

# LOAD DATA ####################################################################
filename_com <- paste0(data_folder, "comephore/zoom_5km.csv")
comephore_raw <- read.csv(filename_com, sep = ",")
filename_loc <- paste0(data_folder, "comephore/loc_px_zoom_5km.csv")
loc_px <- read.csv(filename_loc, sep = ",")

df_comephore <- comephore_raw
head(df_comephore)

# Take only data after 2007
colnames(df_comephore)[1] <- "date"
df_comephore <- df_comephore[df_comephore$date >= "2008-01-01", ]

# put date in index
rownames(df_comephore) <- df_comephore$date
comephore <- df_comephore[-1] # remove dates column


# # get wind data
filename_wind <- paste0(data_folder, "wind/data_gouv/wind_mtp.csv")
wind_mtp <- read.csv(filename_wind)

# Convert datetime to POSIXct
wind_mtp$datetime <- as.POSIXct(wind_mtp$datetime,
                                format = "%Y-%m-%d %H:%M:%S", tz = "UTC")

# Apply function to the DD column
wind_mtp$cardDir <- sapply(wind_mtp$DD, convert_to_cardinal)
wind_mtp$cardDir <- as.character(wind_mtp$cardDir)  # Ensure it's character
wind_mtp$cardDir[is.na(wind_mtp$DD)] <- NA
# Check if NA values are properly handled
summary(wind_mtp)

head(wind_mtp$cardDir)

# DISTANCE AND COORDS ##########################################################

# Get distances matrix
# dist_mat <- get_dist_mat(loc_px)

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

df_dist <- reshape_distances(dist_mat)
df_dist_km <- df_dist
df_dist_km$value <- ceiling(round(df_dist$value, 1)) / 1000


# CHOOSE EXTREME EPISODE FOR R-PARETO ##########################################

q <- 0.998 # quantile
# get central site from sites_coords

set_st_excess <- get_spatiotemp_excess(comephore, q)

list_s <- set_st_excess$list_s
unique(unlist(list_s))
list_t <- set_st_excess$list_t

comephore[list_t[[1]], list_s[[1]]]

list_u <- set_st_excess$list_u
list_u[[1]]

# Spatio-temporal neighborhood parameters
min_spatial_dist <- 5 # in km
delta <- 30 # in hours
episode_size <- delta # size of the episode
s0t0_set <- get_s0t0_pairs(grid_coords_km, comephore,
                            min_spatial_dist = min_spatial_dist,
                            episode_size = episode_size,
                            set_st_excess = set_st_excess,
                            n_max_episodes = 10000,
                            latlon = FALSE)

selected_points <- s0t0_set

# Threshold histogram
df_threshold <- data.frame(u_s0 = selected_points$u_s0)
breaks <- seq(floor(min(df_threshold$u_s0)), ceiling(max(df_threshold$u_s0)), by = 0.1)

ggplot(df_threshold, aes(x = u_s0)) +
  geom_histogram(breaks = breaks, fill = btfgreen, alpha = 0.5, color = "#5f5d5d") +
  btf_theme +
  xlab(TeX(paste0("Threshold for quantile $q = ", q, "$"))) +
  ylab("Count")
filename <- paste(im_folder, "optim/comephore/3km_threshold_histogram_q",
                  q * 1000, "_delta", delta, "_dmin", min_spatial_dist,
                  ".png", sep = "")
ggsave(filename, width = 20, height = 15, units = "cm")



n_episodes <- length(selected_points$s0)
print(n_episodes)
length(unique(selected_points$s0)) # can be same s0
length(unique(selected_points$t0)) # never same t0?
print(min(selected_points$u_s0)) # min threshold
t0_list <- selected_points$t0
s0_list <- selected_points$s0
list_episodes_points <- get_extreme_episodes(selected_points, comephore,
                                     episode_size = episode_size, unif = FALSE,
                                     beta = 0)

list_episodes <- list_episodes_points$episodes
episode <- list_episodes[[1]]
u <- u_list[[1]] 
# check the episode
head(episode)
library(ggplot2)
library(reshape2)  # for melting wide data to long format

# Convert matrix to data frame
index <- 10
sort(t0_list)
which(t0_list == t0_list[index])
episode_test <- list_episodes[[index]]
df_episode <- as.data.frame(episode_test)
df_episode$Time <- 0:(nrow(df_episode) - 1)  # Add a time column
s0_list[index]
u_episode <- selected_points$u_s0[index]
# t0_episode <- t0_list[index]
# Convert from wide to long format
df_long <- melt(df_episode, id.vars = "Time")
head(df_long)
colnames(df_long) <- c("Time", "Pixel", "Value")
ggplot(df_long, aes(x = Time, y = Value, group = Pixel)) +
  geom_line(color = btfgreen) +
  geom_hline(yintercept = u_episode, color = "red", linetype = "dashed") +
  theme_minimal()

# filename <- paste(im_folder, "optim/comephore/extreme_episode", index, "_min", min_spatial_dist,
#                   "km_max", tmax, "h_delta_", delta, ".png", sep = "")
# # filename <- "test.png"
# ggsave(filename, width = 20, height = 15, units = "cm")


list_episodes_unif_points <- get_extreme_episodes(selected_points, comephore,
                                      episode_size = episode_size, unif = TRUE)

list_episodes_unif <- list_episodes_unif_points$episodes

s0_list <- selected_points$s0
t0_list <- selected_points$t0
u_list <- selected_points$u_s0

library(parallel)

tau_vect <- 0:10
tmax <- max(tau_vect)
df_coords <- as.data.frame(grid_coords_m)

# ADD WIND DATA ################################################################

compute_wind_episode_multi <- function(episode, wind_df, s0_name, u) {
  # Récupérer timestamps
  timestamps <- as.POSIXct(rownames(episode), tz = "UTC")
  episode$timestamp <- timestamps
  
  # Filtrer pour dépassements
  excess_rows <- episode[episode[[s0_name]] > u, ]
  
  if (nrow(excess_rows) == 0) {
    return(data.frame(
      n_exces = 0,
      FF_mean = NA, DD_mean = NA, 
      FXY_max = NA, DXY_max = NA,
      cardDir_max_gust = NA
    ))
  }
  
  # Extraire vent à ces timestamps
  wind_subset <- wind_df %>% filter(datetime %in% excess_rows$timestamp)
  
  if (nrow(wind_subset) == 0) {
    return(data.frame(
      n_exces = nrow(excess_rows),
      FF_mean = NA, DD_mean = NA, 
      FXY_max = NA, DXY_max = NA,
      cardDir_max_gust = NA
    ))
  }
  
  # Moyenne direction corrigée (circulaire)
  theta <- wind_subset$DD * pi / 180
  mean_theta <- atan2(mean(sin(theta), na.rm = TRUE), mean(cos(theta), na.rm = TRUE))
  mean_theta_deg <- (mean_theta * 180 / pi) %% 360
  
  max_idx <- which.max(wind_subset$FXY)
  
  return(data.frame(
    n_exces = nrow(excess_rows),
    FF_mean = mean(wind_subset$FF, na.rm = TRUE),
    DD_mean = mean_theta_deg,
    FXY_max = max(wind_subset$FXY, na.rm = TRUE),
    DXY_max = wind_subset$DXY[max_idx],
    cardDir_max_gust = wind_subset$cardDir[max_idx]
  ))
}

wind_per_episode <- Map(compute_wind_episode_multi, list_episodes, s0_list, u_list,
                        MoreArgs = list(wind_df = wind_mtp))

wind_ep_1 <- compute_wind_episode_multi(list_episodes[[3]], wind_mtp, s0_list[3], u_list[3])

wind_per_episode <- Map(compute_wind_episode, list_episodes, s0_list, u_list,
               MoreArgs = list(wind_df = wind_mtp, speed_time = 0))

wind_ep_df <- do.call(rbind, wind_per_episode)
head(wind_ep_df)

wind_ep_df[wind_ep_df$n_exces > 1,]


compute_wind_episode_all <- function(episode, episode_id, wind_df) {
  # Vérifie que les rownames sont bien des timestamps
  timestamps <- as.POSIXct(rownames(episode),
                           format = "%Y-%m-%d %H:%M:%S", tz = "UTC")

  # Subset le vent sur ces timestamps
  wind_subset <- wind_df %>%
    dplyr::filter(datetime %in% timestamps) %>%
    dplyr::mutate(episode_id = episode_id,
                  time_rel = seq_along(datetime))  # Index relatif dans l'épisode

  return(wind_subset)
}


all_wind_distributions <- lapply(seq_along(list_episodes_unif), function(i) {
  compute_wind_episode_all(list_episodes_unif[[i]], episode_id = i, wind_mtp)
})

df_wind_episodes <- do.call(rbind, all_wind_distributions)


ep_id <- 30
df_plot <- df_wind_episodes %>%
  filter(episode_id == ep_id)

ggplot(df_plot, aes(x = factor(cardDir))) +
  geom_bar(fill = btfgreen) +
  coord_polar(start = -pi/8) +
  theme_minimal() +
  labs(title = paste("Episode", ep_id),
       x = "Direction", y = "Frequence")


wind_stats <- df_wind_episodes %>%
  group_by(episode_id) %>%
  summarise(mean_FF = mean(FF, na.rm = TRUE),
            dominant_cardDir = names(sort(table(cardDir), decreasing = TRUE))[1],
            .groups = "drop")


df_plot <- df_wind_episodes %>%
    filter(episode_id == ep_id) %>%
  mutate(
    u = -FF * sin(DD * pi / 180),
    v = -FF * cos(DD * pi / 180),
    t = seq_along(datetime)  # index temporel
  )

ggplot(df_plot, aes(x = t)) +
  geom_segment(aes(xend = t + u, y = 0, yend = v),
               arrow = arrow(length = unit(0.2, "cm")),
               color = "steelblue") +
  coord_fixed() +
  theme_minimal() +
  labs(
    title = paste("Episode", ep_id),
    x = "Time (index)",
    y = "Component N-S"
  )

filename <- paste(im_folder, "wind/episode/", 
                  "wind_vector_episode_", ep_id, ".png", sep = "")
ggsave(filename, width = 20, height = 15, units = "cm")

ep_id <- 500
episode <- list_episodes[[ep_id]]
s0 <- s0_list[ep_id]  # Choose the first site for the example
u <- u_list[ep_id]  # Threshold for the episode

df_rain <- data.frame(
  datetime = as.POSIXct(rownames(episode), format = "%Y-%m-%d %H:%M:%S", tz = "UTC"),
  rain = episode[[s0]]
)

df_wind_episode <- df_wind_episodes %>%
  filter(episode_id == ep_id)

df_plot <- df_rain %>%
  left_join(df_wind_episode, by = "datetime")

# DD in radians
df_plot <- df_plot %>%
    mutate(DD = DD * pi / 180)

df_long <- df_plot %>%
  pivot_longer(cols = c("rain", "FF", "DD"),
               names_to = "variable", values_to = "value") %>%
  mutate(variable = recode(variable,
                           "rain" = "Rain in s0 (mm/h)",
                           "FF"   = "Wind speed (m/s)",
                           "DD"   = "Wind direction (radians)"))

# Plot
ggplot(df_long, aes(x = datetime, y = value, color = variable)) +
  geom_line(size = 1) +
  geom_abline(intercept = u, slope = 0, color = "red", linetype = "dashed") +
  scale_color_manual(values = c(
    "Rain in s0 (mm/h)" = "dodgerblue",
    "Wind speed (m/s)" = "darkgreen",
    "Wind direction (radians)" = "#9c5447"
  )) +
  theme_minimal() +
  labs(
    title = paste0("Episode ", ep_id),
    x = "Time",
    y = "Value",
    color = "Variable"
  ) 

# Save plot
filename <- paste(im_folder, "wind/episode/", 
                  "wind_episode_", ep_id, ".png", sep = "")
ggsave(filename, width = 20, height = 15, units = "cm")






compute_wind_episode_vector_mean <- function(episode, wind_df) {
  # Récupérer les timestamps de l'épisode
  timestamps <- as.POSIXct(rownames(episode), tz = "UTC")
  episode$timestamp <- timestamps
  
  # Extraire les données de vent correspondant aux timestamps
  wind_subset <- wind_df %>% filter(datetime %in% episode$timestamp)
  
  if (nrow(wind_subset) == 0) {
    return(data.frame(
      n_obs = 0,
      FF_vector_mean = NA,
      DD_vector_mean = NA,
      FXY_max = NA,
      DXY_max = NA,
      cardDir_max_gust = NA
    ))
  }

  # Convertir direction en radians
  theta_rad <- wind_subset$DD * pi / 180
  
  # Calcul des composantes u (Est) et v (Nord)
  u <- wind_subset$FF * cos(theta_rad)
  v <- wind_subset$FF * sin(theta_rad)
  
  # Moyenne des composantes
  u_mean <- mean(u, na.rm = TRUE)
  v_mean <- mean(v, na.rm = TRUE)
  
  # Force et direction du vent moyen résultant
  FF_vector_mean <- sqrt(u_mean^2 + v_mean^2)
  DD_vector_mean <- (atan2(v_mean, u_mean) * 180 / pi) %% 360

  # Rafale maximale
  max_idx <- which.max(wind_subset$FXY)
  
  return(data.frame(
    n_obs = nrow(wind_subset),
    FF_vector_mean = FF_vector_mean,
    DD_vector_mean = DD_vector_mean,
    FXY_max = max(wind_subset$FXY, na.rm = TRUE),
    DXY_max = wind_subset$DXY[max_idx],
    cardDir_max_gust = wind_subset$cardDir[max_idx]
  ))
}

wind_per_episode <- Map(compute_wind_episode_vector_mean, list_episodes,
                        MoreArgs = list(wind_df = wind_mtp))

wind_ep_df <- do.call(rbind, wind_per_episode)
head(wind_ep_df)

wind_ep_df[wind_ep_df$n_exces > 1,]

