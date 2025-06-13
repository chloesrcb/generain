# remove objects from the environment
rm(list = ls())

# clear console
cat("\014")

# Load libraries and set theme
source("./script/load_libraries.R")
source("./script/downscaling/pinnEV.R")
library(sf)
library(sp)
library(geosphere)

# LOAD DATA ####################################################################
# Meteo France COMEPHORE data
filename_com <- paste0(data_folder, "comephore/inside_mtp.csv")
filename_loc <- paste0(data_folder, "comephore/loc_pixels_mtp.csv")
comephore_raw <- read.csv(filename_com, sep = ",")
loc_px <- read.csv(filename_loc, sep = ",")
colnames(loc_px) <- c("pixel_name", "Longitude", "Latitude")
# ncol(comephore_raw)
# Convert loc_px to an sf object with Lambert-93 CRS
# loc_sf <- st_as_sf(loc_px, coords = c("Longitude", "Latitude"), crs = 4326)

# # Transform to WGS84 (EPSG:4326)
# loc_wgs84 <- st_transform(loc_sf, 4326)

# # Extract transformed coordinates and add them back to the dataframe
# loc_px$Longitude <- st_coordinates(loc_wgs84)[, 1]  # X = Longitude
# loc_px$Latitude <- st_coordinates(loc_wgs84)[, 2]   # Y = Latitude

# View updated dataframe
head(loc_px)

# Get the date after 2007
comephore_raw$date <- as.POSIXct(comephore_raw$date, format = "%Y-%m-%d %H:%M:%S", tz = "GMT")
df_comephore <- comephore_raw[comephore_raw$date >= "2008-01-01", ]
rain_com <- df_comephore
colnames(rain_com)[1] <- "date"

# Get distances matrix
dist_mat_com <- get_dist_mat(loc_px)
df_dist_com <- reshape_distances(dist_mat_com)

# OMSEV HSM rainfall data
filename_rain <- paste0(data_folder,
                        "omsev/omsev_5min/rain_mtp_5min_2019_2022.RData")
load(filename_rain)
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

print(location_gauges)

# GET DATES INTERVAL ###########################################################
min_date <- min(rain_hsm$dates)
max_date <- max(rain_hsm$dates)

# Extend the dates by 2 hours
extended_min_date <- min_date - 3600 * 2  # Subtract 3600 * 2 seconds (2 hours)
extended_max_date <- max_date + 3600 * 2 # Add 3600 * 2 seconds (2 hours)

# Ensure the extended dates are in Date format (if rain_com$dates is Date)
extended_min_date <- as.POSIXct(extended_min_date, tz = "GMT")
extended_max_date <- as.POSIXct(extended_max_date, tz = "GMT")

# Filter the COMEPHORE data to less dates
filtered_comephore <- rain_com %>%
  dplyr::filter(date >= extended_min_date & date <= extended_max_date)
# filtered_comephore$date[1]
# GET TABLES OF RESPONSE and PREDICTORS ########################################

# Initialize the sparse grid
grid_df <- data.frame(
    time = character(),
    lon_Y = numeric(),
    lat_Y = numeric(),
    lon_X = numeric(),
    lat_X = numeric(),
    Y_obs = numeric()
)

# Add the X1, X2, ..., X27 columns
for (i in 1:27) {
    grid_df[[paste0("X", i)]] <- numeric()
}

colnames(grid_df)
output_file <- "./data/downscaling_table.csv"
# Create a function to write chunks to a file
write_to_file <- function(data, file_path, append = FALSE) {
    write.table(
        data, file_path, row.names = FALSE, 
        col.names = !append, sep = ";", append = append
    )
}

write_to_file(grid_df, output_file, append = FALSE)

# Remove rows with all NA values
rain_hsm_nona <- rain_hsm[!apply(rain_hsm[, -which(names(rain_hsm) == "dates")],
                                          1, function(x) all(is.na(x))), ]


# Create the CSV if it doesn't exist
if (!file.exists(output_file)) {
    write.csv(grid_df, output_file, row.names = FALSE, sep = ";")
}

# Pre-calculate the filtered comephore dates once outside the loops
# filtered_comephore$date <- as.POSIXct(filtered_comephore$date, 
#                                     format = "%Y-%m-%d %H:%M:%S", tz = "GMT")

filtered_comephore$date[1]
# Initialize a list to hold rows temporarily
chunk_list <- list()
chunk_size <- 1000  # Number of rows to accumulate before writing

# Loop over the rows in rain_hsm (time-based)
for (t in 1:nrow(rain_hsm_nona)) {
  t_row <- rain_hsm_nona[t, ]
  t_date <- as.POSIXct(t_row$dates, format = "%Y-%m-%d %H:%M:%S %Z",
                      tz = "GMT")

  # Loop over the rows in location_gauges (site-based)
  for (i in 1:nrow(location_gauges)) {
    site <- location_gauges[i, ]
    Y_obs <- t_row[[site$codestation]]

    if (is.na(Y_obs)) next  # Skip if Y_obs is NA

    # Get closest pixel coordinates
    coord_x_Xcenter <- site$coord_x_px
    coord_y_Xcenter <- site$coord_y_px

    # Distances between X1 and each Comephore pixel (Using distHaversine)
    distances <- distHaversine(cbind(coord_x_Xcenter, coord_y_Xcenter),
                              cbind(loc_px$Longitude, loc_px$Latitude))

    # 3 km around X1 i.e., 9 pixels
    cube_bounds <- loc_px[distances <= 1500, ]
    # nrow(cube_bounds)
    if (nrow(cube_bounds) == 0) next  # Skip if no relevant grid points

    # Temporal bounds (±1.5 hours) i.e., 3 observations by pixel
    t_start <- t_date - 3600 * 1.5
    t_end <- t_date + 3600 * 1.5

    # Find the rows in filtered_comephore that match the temporal bounds
    date_indices <- which(filtered_comephore$date >= t_start &
                          filtered_comephore$date <= t_end)
    if (length(date_indices) == 0) next  # Skip if no relevant observations

    # Filter out the date-specific rows
    cube_obs <- filtered_comephore[date_indices, ]

    # Now filter the relevant pixel names
    pixel_names <- cube_bounds$pixel_name
    relevant_cols <- match(pixel_names, colnames(cube_obs))
    relevant_cols <- relevant_cols[!is.na(relevant_cols)]  # Remove NAs
    if (length(relevant_cols) == 0) next  # Skip if no matching columns

    # Extract the observations for these relevant pixels
    X_obs <- as.vector(as.matrix(cube_obs[, relevant_cols]))

    # Ensure we have exactly 27 values for X_obs
    if (length(X_obs) < 27) {
        X_obs <- c(X_obs, rep(NA, 27 - length(X_obs)))
    }

    # Convert X_obs into named columns X_1 to X_27
    X_obs_named <- setNames(as.list(X_obs), paste0("X", 1:27))

    # Add row to chunk list
    chunk_list[[length(chunk_list) + 1]] <- data.frame(
        time = t_date,
        lon_Y = site$Longitude,
        lat_Y = site$Latitude,
        lon_X = coord_x_Xcenter,
        lat_X = coord_y_Xcenter,
        Y_obs = Y_obs,
        X_obs_named  # Add X_1 to X_27 as separate columns
    )

    # Write chunk to file if chunk size is reached
    if (length(chunk_list) >= chunk_size) {
        write_to_file(do.call(rbind, chunk_list), output_file,
                              append = TRUE)
        chunk_list <- list()  # Clear the chunk
    }
  }
}

# Write remaining rows
if (length(chunk_list) > 0) {
  write_to_file(do.call(rbind, chunk_list), output_file, append = TRUE)
}


# First rows
# print(head(output_df))

# Write the final output_df to a new CSV file
# write.csv(output_df, "output_dataset.csv", row.names = FALSE)
