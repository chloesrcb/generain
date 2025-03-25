# remove objects from the environment
rm(list = ls())

# clear console
cat("\014")

# Load libraries and set theme
source("./script/load_libraries.R")
source("./script/pinnEV.R")
library(sf)
library(sp)
library(geosphere)

# LOAD DATA ####################################################################
# Meteo France COMEPHORE data
comephore_raw <- read.csv("./data/comephore/zoom_3km.csv", sep = ",")
loc_px <- read.csv("./data/comephore/loc_px_zoom_3km.csv", sep = ",")
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



# GET DATA ####################################################################

output_file <- "./data/downscaling_table.csv"

df_Y_X_raw <- read.csv(output_file, sep = ";", header = TRUE)
head(df_Y_X_raw)
# colnames(df_Y_X_raw) <- c("time", "lon_Y", "lat_Y", "lon_X", "lat_X", "Y_obs",
#                           paste0("X", 1:27))

# Rename X_obs to X
# colnames(df_Y_X_raw)[which(colnames(df_Y_X_raw) == "X_obs")] <- "X"
nrow(df_Y_X_raw)
# head(df_Y_X_raw)

# first line
df_Y_X_raw[1,]
df_Y_X <- df_Y_X_raw

# Keep time as a POSIXct object
df_Y_X$time <- as.POSIXct(df_Y_X$time,
                         format = "%Y-%m-%d %H:%M:%S", tz = "GMT")

library(lubridate)
df_Y_X$hour <- hour(df_Y_X$time)
df_Y_X$minute <- minute(df_Y_X$time)
df_Y_X$day <- day(df_Y_X$time)
df_Y_X$month <- month(df_Y_X$time)
df_Y_X$year <- year(df_Y_X$time)
head(df_Y_X)


library(keras)
library(tensorflow)
library(reticulate)
py_version <- "3.9.18"
# path_to_python <- reticulate::install_python(version=py_version)
# path_to_python <- "/home/cserreco/.pyenv/versions/3.9.18/bin/python3.9"
# reticulate::virtualenv_remove("pinnEV_env")
# Create a virtual envionment 'pinnEV_env' with Python 3.9.18. Install tensorflow
# within this environment.
# reticulate::virtualenv_create(envname = 'pinnEV_env',
#                               python=path_to_python,
#                               version=py_version)

path<- paste0(reticulate::virtualenv_root(),"/pinnEV_env/bin/python")
Sys.setenv(RETICULATE_PYTHON = "/home/cserreco/.virtualenvs/pinnEV_env/bin/python")
# Sys.setenv(RETICULATE_PYTHON = path) # Set Python interpreter
Sys.setenv(RETICULATE_LOG_LEVEL = "DEBUG")
tf_version="2.13.1"
reticulate::use_virtualenv("pinnEV_env", required = T)
# tensorflow::install_tensorflow(method="virtualenv", envname="pinnEV_env",
#                                version=tf_version) #Install version of tensorflow in virtual environment
# keras::install_keras(method = c("virtualenv"), envname = "pinnEV_env",version=tf_version)

# keras::is_keras_available() #Check if keras is available

# tf$constant("Hello TensorFlow!")
# py_run_string("import tensorflow as tf; print(tf.__version__)")
#Install spektral 1.3.0 - this is for the graph convolutional neural networks. Remove all code hereafter if not necessary.
# reticulate::virtualenv_install("pinnEV_env",
#                                packages = "spektral", version="1.3.0")



df_downscaling_all <- df_Y_X
# Remove rows with Y_obs = 0 for model EGPD of positive rainfall
df_downscaling_nn <- df_downscaling_all[df_downscaling_all$Y_obs > 0, ]
nrow(df_downscaling_nn)

n_observations <- nrow(df_downscaling_nn)
n_sites <- 17
n_predictors <- ncol(df_downscaling_nn) - 2

# remove the first column
df_downscaling_nn <- df_downscaling_nn[, -1]

# Standardize data
df_standardized <- df_downscaling_nn %>%
  mutate(across(where(is.numeric), scale))
head(df_standardized)

df_standardized <- df_downscaling_nn %>%
  # Convert time variables to cyclical features
  mutate(
    hour_sin = sin(2 * pi * hour / 24),
    hour_cos = cos(2 * pi * hour / 24),
    minute_sin = sin(2 * pi * minute / 60),
    minute_cos = cos(2 * pi * minute / 60),
    day_sin = sin(2 * pi * day / 31),  # Max 31 days
    day_cos = cos(2 * pi * day / 31),
    month_sin = sin(2 * pi * (month-1) / 12),
    month_cos = cos(2 * pi * (month-1) / 12)
  ) %>%
  select(-hour, -minute, -day, -month, -year) %>%  # Remove original columns
  # Standardize all numerical predictors (except cyclical time features)
  mutate(across(c(starts_with("X"), "Y_obs", "lon_X",
                        "lon_Y", "lat_X", "lat_Y"), scale))
head(df_standardized)
head(df_downscaling_nn)


X_all <- df_standardized[, -c(which(names(df_standardized) == "Y_obs"),
                           which(names(df_standardized) == "time"))]
# Convert into a 3D array
X_all <- array(unlist(X_all), dim = c(n_observations, n_sites, n_predictors))
Y_obs <- df_standardized$Y_obs

# Check for missing values
sum(is.na(X_all))
# sum(is.na(Y_obs))

# Identify rows with missing values in X_all
rows_with_na <- apply(X_all, 1, function(row) any(is.na(row)))

# Remove corresponding rows from X_all and Y_obs
X_all_clean <- X_all[!rows_with_na, , , drop = FALSE]  # Keep the same dimensions as the original array
Y_obs_clean <- Y_obs[!rows_with_na]
n_obs <- length(Y_obs_clean)
dim(X_all_clean)

# Flatten X_all_clean into a 2D array for glmnet
X_all_clean_flat <- array(X_all_clean, dim = c(n_obs, n_sites * n_predictors))
dim(X_all_clean_flat)

# Fit a LASSO model to select the most important predictors
library(glmnet)
fit_lasso <- cv.glmnet(as.matrix(X_all_clean_flat), Y_obs_clean, alpha = 1)


# Coefficients pour lambda.min (le meilleur lambda)
coef_lambda_min <- coef(fit_lasso, s = "lambda.min")
print(coef_lambda_min)
# Coefficients pour lambda.1se (le lambda à 1 écart-type)
coef(fit_lasso, s = "lambda.1se")

selected_predictors <- which(coef(fit_lasso, s = "lambda.min") != 0)
# Exclude intercept (index 1)
selected_predictors <- selected_predictors[selected_predictors > 1] - 1
X_all_clean_selected_flat <- X_all_clean_flat[, selected_predictors]


# Retransform the selected predictors into a 3D array
X_all_clean_selected <- array(X_all_clean_selected_flat, dim = c(n_obs, n_sites, length(selected_predictors)))
head(X_all_clean_selected)


# # Hyperparameter tuning for the neural network
# library(caret)

# grid <- expand.grid(
#   batch.size = c(10, 20),
#   widths = list(c(4, 2), c(8, 4, 2)),
#   n.ep = c(100, 200)
# )

# results <- apply(grid, 1, function(params) {
#   fit <- eGPD.NN.train(
#     Y.train, Y.valid, X.s, X.k, type = "MLP",
#     n.ep = params$n.ep, batch.size = params$batch.size,
#     widths = params$widths, seed = 1
#   )
#   return(fit$val_loss)  # Retourne la perte de validation
# })

# best_params <- grid[which.min(results), ]



X.nn.k <- X_all_clean # Nn predictors for kappa
X.nn.s <- X_all_clean # Nn predictors for sigma
length(Y_obs_clean)
sum(is.na(Y_obs_clean))

set.seed(123)
n_obs <- length(Y_obs_clean)
train_idx <- sample(seq_len(n_obs), size = 0.8 * n_obs)
Y.train <- Y_obs_clean
Y.valid <- Y_obs_clean
Y.train[-train_idx] <- -1e10
Y.valid[train_idx] <- -1e10
Y.train <- array(Y.train, dim = c(length(Y.train), 1))
Y.valid <- array(Y.valid, dim = c(length(Y.valid), 1))

# # NEW
# train_idx <- sample(seq_len(n_obs), size = 0.8 * n_obs)
# X.train <- X_all_clean_selected[train_idx, , ]
# X.valid <- X_all_clean_selected[-train_idx, , ]
# Y.train <- Y_obs_clean[train_idx]
# Y.valid <- Y_obs_clean[-train_idx]
# Y.train <- array(Y.train, dim = c(length(Y.train), 1))
# Y.valid <- array(Y.valid, dim = c(length(Y.valid), 1))

# X.nn.k <- X_all_clean_selected  # Nn predictors for kappa
# X.nn.s <- X_all_clean_selected # Nn predictors for sigma
# List of predictors
X.s <- list("X.nn.s" = X.nn.s)
X.k <- list("X.nn.k" = X.nn.k)


length(Y.train)

# initialize the neural network model EGPD fitting
# get censoring
censores <- seq(0.25, 0.35, 0.001)
rain_new <- rain_hsm[-1] # remove the first column date
# rain_cp <- rain_new[, c(5,6)]
df_score_ohsm <- choose_censore(rain_new, censores, n_samples = 100)

params_subrain <- get_egpd_estimates(rain_new,
                      left_censoring = 0.3)


mean_kappa <- mean(params_subrain$kappa)
mean_sigma <- mean(params_subrain$sigma)
mean_xi <- mean(params_subrain$xi)

# Train the neural network model
NN.fit <- eGPD.NN.train(Y.train, Y.valid, X.s, X.k, type = "MLP",
      n.ep = 200, batch.size = 10, init.scale = mean_sigma,
      init.kappa = mean_kappa, init.xi = mean_xi, widths = c(4, 2), seed = 1)
reticulate::py_last_error()
NN.fit
# we can change widths to reduce complexity (before c(6, 3))

# Predict with the trained model
out <- eGPD.NN.predict(X.s = X.s, X.k = X.k, NN.fit$model)
# out$pred.xi
# out$pred.kappa
# Results
print("sigma linear coefficients: "); print(round(out$pred.sigma, 5))
print("kappa linear coefficients: "); print(round(out$pred.kappa,2))

head(out$pred.sigma)

sigma_site1 <- out$pred.sigma[, 1]
kappa_site1 <- out$pred.kappa[, 1]

year_idx <- which(df_downscaling_nn$year == 2019)
sigma_site1_2019 <- sigma_site1[year_idx]
kappa_site1_2019 <- kappa_site1[year_idx]
print("Sigma for site 1 in 2019: ")
print(sigma_site1_2019)
