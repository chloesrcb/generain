# remove objects from the environment
rm(list = ls())

# clear console
cat("\014")

# Load libraries and set theme
source("./script/load_libraries.R")
source("./script/pinnEV.R")

# LOAD DATA ####################################################################
# Meteo France COMEPHORE data
comephore_raw <- read.csv("./data/comephore/inside_mtp.csv", sep = ",")
loc_px <- read.csv("./data/comephore/loc_pixels_mtp.csv", sep = ",")

# Get the date after 2007
comephore_raw$date <- as.Date(comephore_raw$date, format = "%Y-%m-%d %H:%M:%S")
df_comephore <- comephore_raw[comephore_raw$date >= "2008-01-01", ]
rain_com <- df_comephore
colnames(rain_com)[1] <- "date"

# Get distances matrix
dist_mat_com <- get_dist_mat(loc_px)
df_dist_com <- reshape_distances(dist_mat_com)

# OMSEV HSM rainfall data
load("./data/PluvioMontpellier_1min/rain_mtp_5min_2019_2022.RData")
rain_hsm <- rain.all5[c(1, 6:ncol(rain.all5))]
location_gauges <- read.csv("./data/PluvioMontpellier_1min/pluvio_mtp_loc.csv")
location_gauges$codestation <- c("iem", "mse", "poly", "um", "cefe", "cnrs",
                                 "crbm", "archiw", "archie", "um35", "chu1",
                                 "chu2", "chu3", "chu4", "chu5", "chu6", "chu7")

# Get distances matrix
dist_mat_hsm <- get_dist_mat(location_gauges)
df_dist_hsm <- reshape_distances(dist_mat_hsm)

# GET CORRESPONDING PIXELS ####################################################
# Compute distance between each rain gauge and each pixel
distances <- as.data.frame(matrix(NA, nrow = nrow(location_gauges),
                            ncol = nrow(loc_px)))

for (i in 1:nrow(location_gauges)) {
  for (j in 1:nrow(loc_px)) {
    distances[i, j] <- distHaversine(
      c(location_gauges$Longitude[i], location_gauges$Latitude[i]),
      c(loc_px$Longitude[j], loc_px$Latitude[j])
    )
  }
}

# Assign column and row names for easier interpretation
colnames(distances) <- loc_px$pixel_name
rownames(distances) <- location_gauges$codestation

# Find the closest pixel for each rain gauge
closest_pixels <- apply(distances, 1,
                        function(x) colnames(distances)[which.min(x)])

# Add the closest pixel to location_gauges
location_gauges$closest_pixel <- closest_pixels
location_gauges <- location_gauges[-1]
# Save the result or inspect
print(location_gauges)

location_gauges$coord_x_px <- loc_px$Longitude[match(closest_pixels,
                                                     loc_px$pixel_name)]
location_gauges$coord_y_px <- loc_px$Latitude[match(closest_pixels,
                                                    loc_px$pixel_name)]

# GET DATES INTERVAL ###########################################################
min_date <- min(rain_hsm$dates)
max_date <- max(rain_hsm$dates)

extended_min_date <- min_date - 3600 * 2  # Subtract 3600 * 2 seconds (2 hours)
extended_max_date <- max_date + 3600 * 2 # Add 3600 * 2 seconds (2 hours)

filtered_comephore <- rain_com %>%
    dplyr::filter(date >= extended_min_date & date <= extended_max_date)
