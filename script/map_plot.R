# remove objects from the environment
rm(list = ls())

# clear console
cat("\014")

# Load libraries and set theme
source("./script/load_libraries.R")
library(ggplot2)
library(ggrepel)
library(dplyr)
library(sf)

# LOAD DATA ####################################################################
# Meteo France COMEPHORE data
filepath <- paste0(data_folder, "/comephore/")
loc_px <- read.csv(paste0(filepath, "coords_pixels.csv"), sep=",")
colnames(loc_px) <- c("pixel_name", "Longitude", "Latitude")

# Get distances matrix
dist_mat_com <- get_dist_mat(loc_px)
df_dist_com <- reshape_distances(dist_mat_com)

# OMSEV HSM rainfall data
location_gauges <- read.csv(paste0(data_folder, "/omsev/loc_rain_gauges.csv"), sep=";")
# location_gauges <- read.csv(paste0(data_folder, "/omsev/pluvio_mtp_loc_till_2022.csv"))
location_gauges$codestation <- c("gimel", "archie", "archiw", "cefe", "chu1",
                                "chu2", "chu3", "chu4", "chu5", "chu6", "chu7",
                                "cines", "cnrs",  "crbm", "hydro", "iem",
                                "mse", "poly", "brives", "um25", "um35")
# location_gauges$y <- -location_gauges$y
# Convert distances lambert93 to WGS84
location_gauges <- st_as_sf(location_gauges, coords = c("y", "x"), crs = 27572)
location_gauges <- st_transform(location_gauges, crs = 4326)

location_gauges$y <- -location_gauges$y

# Créer un objet spatial en Lambert II étendu
location_sf <- st_as_sf(location_gauges, coords = c("x", "y"), crs = 2154)

# Transformer en WGS84
location_wgs84 <- st_transform(location_sf, crs = 4326)
# Get lat lon coordinates
loc_gauges <- data.frame(
  codestation = location_gauges$codestation,
  Longitude = st_coordinates(location_gauges)[, 1],
  Latitude = st_coordinates(location_gauges)[, 2]
)


library(generain)
# Get distances matrix
dist_mat_hsm <- get_dist_mat(location_gauges)
df_dist_hsm <- reshape_distances(dist_mat_hsm)

# PLOT PIXELS AND SITES ########################################################

pixel_size <- 1000  # 1 km x 1 km

# Convert CRS (WGS84 -> Lambert-93)
loc_sf <- st_as_sf(loc_px, coords = c("Longitude", "Latitude"), crs = 2154) 

# Function to create a square polygon around a point (Lambert-93)
create_square <- function(point, size) {
  x <- st_coordinates(point)[1]
  y <- st_coordinates(point)[2]
  half_size <- size / 2
  coords <- matrix(c(
    x - half_size, y - half_size,  # Lower-left
    x + half_size, y - half_size,  # Lower-right
    x + half_size, y + half_size,  # Upper-right
    x - half_size, y + half_size,  # Upper-left
    x - half_size, y - half_size   # Close the square
  ), ncol = 2, byrow = TRUE)

  st_polygon(list(coords))
}

# Create square polygons around each point
pixel_polygons <- loc_sf %>%
  rowwise() %>%
  mutate(geometry = list(create_square(geometry, pixel_size))) %>%
  ungroup() %>%
  st_as_sf()

# Transform to WGS84 (EPSG:4326)
st_crs(pixel_polygons) <- 2154
pixel_polygons <- st_transform(pixel_polygons, crs = 4326)


# Get montpellier map for the background
mtp_geo <- st_read("./data/geometry/VilleMTP_MTP_SousQuartiers.geojson")
mtp_geo <- st_transform(mtp_geo, crs = 4326)

# Plot
ggplot() +
  geom_sf(data = mtp_geo, fill = "lightgrey", alpha = 0.5) +
  geom_sf(data = pixel_polygons, fill = NA, color = btfgreen, size = 3) +
  geom_point(data = location_gauges, aes(x = Longitude, y = Latitude),
             color = "#d13737", size = 2, shape = 16) +
  geom_text_repel(data = location_gauges, aes(x = Longitude, y = Latitude, label = codestation), 
                  color = "#d13737", size = 2) +
  labs(title = "Pixels of 1km x 1km",
       x = "Longitude", y = "Latitude") +
  theme_minimal()
# save plot 
filename <- paste0(im_folder, "/map/pixelsCOM_sitesOMSEV.png")
ggsave(filename, width = 10, height = 10, dpi = 300)

ggplot() +
  geom_sf(data = mtp_geo, fill = "lightgrey", alpha = 0.5) +
  geom_sf(data = pixel_polygons, fill = NA, color = btfgreen, size = 3) +
  geom_point(data = location_gauges, aes(x = Longitude, y = Latitude),
             color = "#d13737", size = 2.5, shape = 19) +
  geom_text_repel(data = location_gauges, aes(x = Longitude, y = Latitude, label = codestation), 
                  color = "#d13737", size = 4) +
  labs(title = "Pixels of 1km x 1km",
       x = "Longitude", y = "Latitude") +
  theme_minimal() +
  coord_sf(xlim = c(3.82, 3.9), ylim = c(43.6, 43.67))
# save plot 
filename <- paste0(im_folder, "/map/pixelsCOM_sitesOMSEV_zoom.png")
ggsave(filename, width = 10, height = 10, dpi = 300)
