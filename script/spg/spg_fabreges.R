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
library(sf)
library(dplyr)
library(lwgeom)

file_b <- file.path(data_folder, "geometry/Mesh_Cournonsec+Fabregues.geojson")

mesh <- st_read(file_b, quiet = TRUE) %>%
  st_transform(2154) %>%
  st_make_valid()
mesh_union <- st_union(mesh) %>%
  st_make_valid()

# Envelope of the bassin
bassin_outline <- st_concave_hull(
  mesh_union,
  ratio = 0.15,
  allow_holes = FALSE
) %>%
  st_make_valid()

cell_m <- 100 # resolution of the grid in meters

margin <- 5 * cell_m # add a margin around the bassin

# Create a grid of square pixels over the bassin
bbox_b <- st_bbox(bassin_outline)
bbox_b_expanded <- bbox_b
bbox_b_expanded[c("xmin", "ymin")] <- bbox_b[c("xmin", "ymin")] - margin
bbox_b_expanded[c("xmax", "ymax")] <- bbox_b[c("xmax", "ymax")] + margin
bbox_sfc <- st_as_sfc(bbox_b_expanded)
# grid
grid_raw <- st_make_grid(
  bbox_sfc,
  cellsize = cell_m,
  square = TRUE
) %>%
  st_as_sf() 

# Select only the pixels that intersect with the bassin
grid_square_l93 <- grid_raw[
  lengths(st_intersects(grid_raw, bassin_outline)) > 0,
] %>%
  mutate(pixel_id = paste0("pixel_", row_number()))

# select the conditioning pixel
s0_pixel_id <- grid_square_l93$pixel_id[1400]
s0_geom <- grid_square_l93[grid_square_l93$pixel_id == s0_pixel_id, ]

ggplot() +
  geom_sf(data = grid_square_l93,
          fill = "grey95", color = "grey75", linewidth = 0.15) +
  geom_sf(data = s0_geom,
          fill = NA, color = "red", linewidth = 1.5) +
  geom_sf(data = st_as_sf(bassin_outline),
          fill = NA, color = "black", linewidth = 1) +
  coord_sf(expand = FALSE) +
  theme_void()

# ggsave(
#   filename = paste0(im_folder, "swg/fabregues_grid_pixels_inside_basin_l93.png"),
#   plot = p_grid,
#   width = 8,
#   height = 6,
#   dpi = 300,
#   bg = "white"
# )

# get the centroids of the pixels
grid_cent_l93 <- st_centroid(grid_square_l93)
grid_cent_l93$pixel_id <- grid_square_l93$pixel_id
grid_xy <- st_coordinates(grid_cent_l93)
# coordinates of the pixels in a data frame
grid_df_l93 <- data.frame(
  x = grid_xy[, 1],
  y = grid_xy[, 2]
)
rownames(grid_df_l93) <- grid_cent_l93$pixel_id
grid_df_l93$pixel_id <- grid_cent_l93$pixel_id

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
adv_5m <- c(100, 200)

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
  seed = NULL,
  cell_m = 100
)



# To SW2D format --------------------------------------------------------------
dt_sec <- 300  # 5 minutes

# mm/5min -> m/s
sim_episode_mps <- sim_episode
sim_episode_mps <- sim_episode_mps * 1e-3 / dt_sec

mat <- sim_episode_mps
mat <- sim_episode_mps
stations <- rownames(mat)
times <- seq(0, by = 300, length.out = ncol(mat))  # 5 min = 300s

df_sw2d <- data.frame(
  t = times,
  t(mat) 
)

file_out <- "precipitation_time_series.txt"

N <- nrow(mat)

folder_out <- paste0(data_folder, "spg/fabregues/")
if (!dir.exists(folder_out)) {
  dir.create(folder_out, recursive = TRUE)
}

con <- file(paste0(folder_out, file_out), open = "w")

writeLines("#Precipitation time series", con)
writeLines("#Precipitation stations", con)
writeLines(as.character(N), con)
writeLines("#", con)
writeLines("#t + N times (<Tab> P)", con)

apply(df_sw2d, 1, function(row) {
  writeLines(paste(row, collapse = "\t"), con)
})

close(con)

# get the coordinates of the stations and save them to a CSV file
write.csv(grid_df_l93,
          paste0(folder_out, "precipitation_station_coordinates.csv"),
          row.names = FALSE)

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

str_adv <- paste0("adv_", adv_5m[1], "_", adv_5m[2])
episode_idx <- str_adv
dir_frames <- file.path(im_folder, paste0("swg/fabregues/frames_episode_", episode_idx))
dir.create(dir_frames, recursive = TRUE, showWarnings = FALSE)

cent_df <- data.frame(
  pixel_id = grid_square_l93$pixel_id,
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

  map_t <- grid_square_l93 %>%
    left_join(df_t, by = "pixel_id") %>%
    mutate(rain_plot = rain)

  tmp <- df_t %>%
    left_join(cent_df, by = "pixel_id") %>%
    mutate(w = pmax(rain, 0))

  has_rain <- any(tmp$w > threshold / 4, na.rm = TRUE)

  p <- ggplot() +
    geom_sf(data = bassin_geom, fill = NA, color = "black", linewidth = 0.8) +

    geom_sf(data = map_t, aes(fill = rain_plot), color = NA) +
    geom_sf(data = s0_geom, fill = NA, color = "#ee8686", linewidth = 1.2) +
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





