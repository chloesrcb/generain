rm(list = ls())

# clear console
cat("\014")

library(data.table)
library(latex2exp)
library(lubridate)
library(fuzzyjoin)
library(grid)

muse <- FALSE

if (muse) {
  # Get the muse folder
  folder_muse <- "/home/serrec/work_rainstsimu/rpareto/"
  setwd(folder_muse)
  # Load libraries and set theme
  source("load_libraries.R")
  im_folder <- "./images"
  source("config_omsev.R")
  data_folder <- "./data/"
  ncores <- 27
} else {
  # Load libraries and set theme
  source("./script/load_libraries.R")
  source("./script/optimisation/config_omsev.R")
  ncores <- detectCores() - 1
}

# Get all files in the folder "R"
functions_folder <- "./R"
files <- list.files(functions_folder, full.names = TRUE)
# # load all functions in files
invisible(lapply(files, function(f) source(f, echo = FALSE)))

filename_loc <- paste0(data_folder,
                           "omsev/loc_rain_gauges.csv")
# get location of each rain gauge
location_gauges <- read.csv(filename_loc)

# get rain data
filename_rain <- paste0(data_folder, "omsev/omsev_5min/rain_mtp_5min_2019_2025.csv")
rain_omsev <- read.csv(filename_rain)
head(rain_omsev)
tail(rain_omsev)
# put dates as rownames
rownames(rain_omsev) <- rain_omsev$dates
# remove dates after january 2025
# rain_omsev <- rain_omsev[rownames(rain_omsev) < "2025-01-31 00:00:00", ]
rain <- rain_omsev[-1] # remove dates column

# remove cines, hydro, brives
# rain <- rain[, !(colnames(rain) %in% c("cines", "hydro", "brives"))]
# location_gauges <- location_gauges[location_gauges$Station != "cines" &
#                                    location_gauges$Station != "hydro" &
#                                    location_gauges$Station != "brives", ]
dist_mat <- get_dist_mat(location_gauges)
df_dist <- reshape_distances(dist_mat)

sites_names <- colnames(rain)

sites_coords <- location_gauges[, c("Longitude", "Latitude")]

rownames(sites_coords) <- location_gauges$Station
sites_coords_sf <- st_as_sf(sites_coords, coords = c("Longitude", "Latitude"),
                            crs = 4326)
sites_coords_sf <- st_transform(sites_coords_sf, crs = 2154)
coords_m <- st_coordinates(sites_coords_sf)
grid_coords_m <- as.data.frame(coords_m)
colnames(grid_coords_m) <- c("Longitude", "Latitude")
rownames(grid_coords_m) <- rownames(sites_coords)
grid_coords_km <- as.data.frame(coords_m / 1000)
colnames(grid_coords_km) <- c("Longitude", "Latitude")
rownames(grid_coords_km) <- rownames(sites_coords)


