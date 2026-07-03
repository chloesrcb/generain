rm(list = ls())

# clear console
cat("\014")

# Load libraries and set theme
source("./script/load_libraries.R")

# library(generain)
library(animation)
library(tidyr)
library(lubridate)
library(purrr)
library(viridis)
library(magick)
library(sf)
library(units)
library(ggspatial)
library(dplyr)
library(tidyr)
library(ggplot2)
library(forcats)

# load the functions file
source("./R/spg.R", echo = FALSE)
source("./R/distances.R", echo = FALSE)
library(latex2exp)

################################################################################
# DATA AND COORDS
################################################################################

# get location of each rain gauge
filename_loc <- paste0(data_folder, "omsev/loc_rain_gauges.csv")
location_gauges <- read.csv(filename_loc)

dist_mat <- get_dist_mat(location_gauges)
df_dist <- reshape_distances(dist_mat)


sites_names <- colnames(rain)
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

rownames(grid_coords_m) <- rownames(sites_coords)


# Grid inside Verdanson basin
file_bv <- file.path(data_folder, "geometry/verdanson_basin.geojson")
bassin_geom <- st_read(file_bv, quiet = TRUE) %>%
  st_transform(2154) %>%
  st_make_valid()

bassin_union <- st_union(bassin_geom)

# Rain gauges
sites_sf <- location_gauges %>%
  st_as_sf(coords = c("Longitude", "Latitude"), crs = 4326)

sites_l93 <- st_transform(sites_sf, 2154) # convert to Lambert-93
cell_m <- 100 # pixel size in meters

# Regular square grid covering the basin bbox
grid_raw <- st_make_grid(
  bassin_union,
  cellsize = cell_m,
  square = TRUE
) %>%
  st_as_sf()

grid_raw$pixel_id <- paste0("pixel_", seq_len(nrow(grid_raw)))

# Keep pixels whose centroid is inside the basin
grid_cent_raw <- st_centroid(grid_raw)
inside_basin <- lengths(st_within(grid_cent_raw, bassin_union)) > 0
grid_square_l93 <- grid_raw[inside_basin, ]
grid_cent_l93   <- grid_cent_raw[inside_basin, ]
grid_poly_l93 <- st_intersection(grid_square_l93, bassin_union)

grid_poly_l93$pixel_id <- grid_square_l93$pixel_id
grid_cent_l93$pixel_id <- grid_square_l93$pixel_id

grid_xy <- st_coordinates(grid_cent_l93)

grid_df_l93 <- data.frame(
  x = grid_xy[, 1],
  y = grid_xy[, 2]
)
rownames(grid_df_l93) <- grid_cent_l93$pixel_id

# Conditioning pixel near the centre of the basin
centre_l93 <- st_centroid(bassin_union)

s0_pixel_id <- "pixel_2100" # referenced s0 pixel id

s0_poly_l93 <- grid_poly_l93 %>%
  filter(pixel_id == s0_pixel_id)

# Check grid
p_grid <- ggplot() +
  geom_sf(data = grid_poly_l93, fill = "grey95", color = "grey75", linewidth = 0.25) +
  geom_sf(data = bassin_geom, fill = NA, color = "black", linewidth = 0.8) +
  geom_sf(data = sites_l93, color = "red", size = 2) +
  geom_sf(data = s0_poly_l93, fill = NA, color = "blue", linewidth = 1.2) +
  coord_sf(expand = FALSE) +
  theme_void() +
  theme(
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA)
  )

print(p_grid)

ggsave(
  filename = paste0(im_folder, "swg/omsev/grid_pixels_inside_basin_l93.png"),
  plot = p_grid,
  width = 8,
  height = 6,
  dpi = 300,
  bg = "white"
)

#############################################################################################
# SIMULATION
#############################################################################################

# Parameters for simulation
u_emp <- 1 # threshold for simulation
nT <- 12 # time steps
steps <- 0:(nT - 1)

# marginal parameters for simulation
params_margins_common <- list(
  p0 = 0.989,
  xi = 0.244,
  sigma = 0.536,
  kappa = 0.308
)

# variogram parameters from optimization results
params_m5min <- list(
  beta1 = 0.48,
  beta2 = 0.77,
  alpha1 = 0.12,
  alpha2 = 0.65
)

# random advection for one episode
adv_5m <- c(100, -200)

# Simulate random episode
sim_episode <- sim_episode_grid_m5(
  params_vario = params_m5min,
  params_margins_common = params_margins_common,
  coords = grid_df_l93,
  steps = steps,
  adv_m5 = adv_5m,
  t0 = 0,
  s0_pixel_id = s0_pixel_id,
  u_emp = u_emp,
  seed = 2026,
  cell_m = 100
)

head(sim_episode)


###############################################################################################
# Ploting the simulated episode
###############################################################################################
threshold <- u_emp
positive_vals <- sim_episode[is.finite(sim_episode) & sim_episode > 0]

fill_max <- if (length(positive_vals) > 0) {
  quantile(positive_vals, 0.995, na.rm = TRUE)
} else {
  threshold
}

fill_max <- max(fill_max, 2 * threshold)
fill_limits <- c(0, fill_max)

episode_idx <- "new"
dir_frames <- file.path(im_folder, paste0("swg/omsev/frames_episode_", episode_idx))
dir.create(dir_frames, recursive = TRUE, showWarnings = FALSE)

cent_df <- data.frame(
  pixel_id = grid_poly_l93$pixel_id,
  x = grid_xy[, 1],
  y = grid_xy[, 2]
)

# adv_5m <- round(adv_5m, 3)
dx_m <- adv_5m[1]
dy_m <- adv_5m[2]

for (tt in seq_len(nT)) {

  df_t <- data.frame(
    pixel_id = rownames(sim_episode),
    rain = sim_episode[, tt]
  )

  map_t <- grid_poly_l93 %>%
    left_join(df_t, by = "pixel_id") %>%
    mutate(rain_plot = rain)

  tmp <- df_t %>%
    left_join(cent_df, by = "pixel_id") %>%
    mutate(w = pmax(rain, 0))

  has_rain <- any(tmp$w > threshold / 4, na.rm = TRUE)

  p <- ggplot() +
    geom_sf(data = bassin_geom, fill = NA, color = "black", linewidth = 0.8) +

    geom_sf(data = map_t, aes(fill = rain_plot), color = NA) +
    geom_sf(data = s0_poly_l93, fill = NA, color = "#ee8686", linewidth = 1.2) +
    geom_sf(data = sites_l93, color = "white", size = 3) +
    geom_sf(data = sites_l93, color = "grey40", size = 2) +
    scale_fill_gradientn(
      colours = c("white", "#dbe9f6", "#9ecae1", "#4a90d9", "#08519c"),
      values = scales::rescale(c(0, 0.2, threshold, 2 * threshold, fill_limits[2])),
      limits = fill_limits,
      oob = scales::squish,
      na.value = "transparent",
      name = paste0("Rainfall (mm/5min)\nThreshold u = ", round(threshold, 2))
    ) +
    labs(
      title = paste0(
        "t = ", tt
      )
    ) +
    theme_void() +
    theme(
  plot.title = element_text(size = 18, hjust = 0.5),
  plot.background = element_rect(fill = "transparent", color = NA),
  panel.background = element_rect(fill = "transparent", color = NA),
  legend.background = element_rect(fill = "transparent", color = NA),
  legend.key = element_rect(fill = "transparent", color = NA),
  legend.box.background = element_rect(fill = "transparent", color = NA),
  legend.text = element_text(size = 14)
)
    coord_sf(expand = FALSE)

  if (has_rain) {
    # put arrow on the maximal value of rain
    # mx <- tmp$x[which.max(tmp$w)]
    # my <- tmp$y[which.max(tmp$w)]
    # mx_pt_l93 <- st_sfc(st_point(c(mx, my)), crs = 2154) %>% st_sf()
    
    # put arrow on the barycenter of rain
    bx <- sum(tmp$x * tmp$w, na.rm = TRUE) / sum(tmp$w, na.rm = TRUE)
    by <- sum(tmp$y * tmp$w, na.rm = TRUE) / sum(tmp$w, na.rm = TRUE)

    bary_pt_l93 <- st_sfc(st_point(c(bx, by)), crs = 2154) %>% st_sf()
    
    arrow_line_l93 <- st_sfc(
      st_linestring(rbind(
        c(bx, by),
        c(bx + dx_m, by + dy_m)
      )),
      crs = 2154
    ) %>%
      st_sf()

    p <- p +
      geom_sf(data = bary_pt_l93, color = "white", size = 2) +
      geom_sf(data = bary_pt_l93, color = "red", size = 1) +
      geom_sf(
        data = arrow_line_l93,
        arrow = arrow(length = unit(0.35, "cm")),
        color = "white",
        linewidth = 2
      ) +
      geom_sf(
        data = arrow_line_l93,
        arrow = arrow(length = unit(0.3, "cm")),
        color = "red",
        linewidth = 1
      )
  }

  out_png <- file.path(dir_frames, sprintf("frame_%03d.png", tt))

  ggsave(
    out_png,
    plot = p,
    width = 7,
    height = 6,
    dpi = 150,
    bg = "transparent"
  )
}

