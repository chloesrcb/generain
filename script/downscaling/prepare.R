# Clear workspace and console
rm(list = ls())
cat("\014")

muse <- TRUE

if (muse) {
  folder_muse <- "/home/serrec/work_rainstsimu/downscaling"
  setwd(folder_muse)
  source("utils.R")
  source("config.R")
  path_to_python <- "/home/serrec/.pyenv/versions/3.9.18/bin/python3.9"
} else {
  source("./script/load_libraries.R")
  source("./script/downscaling/pinnEV.R")
  # If you have config.R, source it; else define data_folder manually
  if (file.exists("./script/config.R")) source("./script/config.R")
  path_to_python <- "/home/cserreco/.pyenv/versions/3.9.18/bin/python3.9"
}


get_dist_mat <- function(df_lonlat) {
  # df_lonlat must have cols Longitude, Latitude
  n <- nrow(df_lonlat)
  out <- matrix(0, n, n)
  for (i in 1:n) {
    out[i, ] <- geosphere::distHaversine(
      cbind(df_lonlat$Longitude[i], df_lonlat$Latitude[i]),
      cbind(df_lonlat$Longitude, df_lonlat$Latitude)
    )
  }
  out
}

reshape_distances <- function(mat) {
  # Optional helper, not strictly used below
  as.data.frame(as.table(mat))
}

add_time_features <- function(df, time_col = "time") {
  tt <- df[[time_col]]
  df %>%
    mutate(
      hour = hour(tt),
      minute = minute(tt),
      day = day(tt),
      month = month(tt),
      year = year(tt),
      hour_sin   = sin(2*pi*hour/24),
      hour_cos   = cos(2*pi*hour/24),
      minute_sin = sin(2*pi*minute/60),
      minute_cos = cos(2*pi*minute/60),
      day_sin    = sin(2*pi*day/31),
      day_cos    = cos(2*pi*day/31),
      month_sin  = sin(2*pi*(month-1)/12),
      month_cos  = cos(2*pi*(month-1)/12)
    ) %>%
    select(-hour, -minute, -day, -month, -year)
}

standardize_3d_by_train <- function(X, train_t_idx) {
  # X: array (T, S, d)
  Xs <- X
  d <- dim(X)[3]
  mus <- numeric(d)
  sds <- numeric(d)
  for (j in 1:d) {
    mus[j] <- mean(X[train_t_idx,,j], na.rm = TRUE)
    sds[j] <- sd(X[train_t_idx,,j], na.rm = TRUE)
    Xs[,,j] <- (X[,,j] - mus[j]) / (sds[j] + 1e-8)
  }
  list(X_scaled = Xs, mu = mus, sd = sds)
}

# BUILD DOWNSCALING TABLE (CSV)

filename_com <- paste0(data_folder, "comephore/comephore_2008_2024_5km.csv")
filename_loc_px <- paste0(data_folder, "comephore/loc_px_zoom_5km.csv")
filename_rain   <- paste0(data_folder, "omsev/omsev_5min/rain_mtp_5min_2019_2022.RData")
filename_loc <- paste0(data_folder, "omsev/loc_rain_gauges.csv")
output_file <- paste0(data_folder, "downscaling/downscaling_table.csv")


dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)

# Load COMEPHORE
comephore_raw <- read.csv(filename_com, sep = ",")
loc_px <- read.csv(filename_loc_px, sep = ",")
colnames(loc_px) <- c("pixel_name", "Longitude", "Latitude")

comephore_raw$date <- as.POSIXct(comephore_raw$date, format = "%Y-%m-%d %H:%M:%S", tz = "GMT")
df_comephore <- comephore_raw %>% filter(date >= as.POSIXct("2008-01-01", tz="GMT"))
rain_com <- df_comephore
colnames(rain_com)[1] <- "date"

# Load OMSEV
load(filename_rain) # loads rain.all5
rain_hsm <- rain.all5[c(1, 6:ncol(rain.all5))]
colnames(rain_hsm)[1] <- "dates"
rain_hsm$dates <- as.POSIXct(rain_hsm$dates, tz="GMT")

# Load gauge locations
# get location of each rain gauge
location_gauges <- read.csv(filename_loc)
location_gauges$Station <- c("iem", "mse", "poly", "um", "cefe", "cnrs",
                             "crbm", "archiw", "archie", "um35", "chu1",
                             "chu2", "chu3", "chu4", "chu5", "chu6", "chu7",
                             "cines", "brives", "hydro")

# Closest pixel per gauge
distances_gp <- matrix(NA_real_, nrow = nrow(location_gauges), ncol = nrow(loc_px))
for (i in 1:nrow(location_gauges)) {
  distances_gp[i,] <- geosphere::distHaversine(
    cbind(location_gauges$Longitude[i], location_gauges$Latitude[i]),
    cbind(loc_px$Longitude, loc_px$Latitude)
  )
}
colnames(distances_gp) <- loc_px$pixel_name
rownames(distances_gp) <- location_gauges$codestation
closest_pixels <- apply(distances_gp, 1, function(x) colnames(distances_gp)[which.min(x)])

location_gauges$closest_pixel <- closest_pixels
location_gauges$coord_x_px <- loc_px$Longitude[match(closest_pixels, loc_px$pixel_name)]
location_gauges$coord_y_px <- loc_px$Latitude[match(closest_pixels, loc_px$pixel_name)]

# Restrict COMEPHORE to OMSEV date range ±2h
min_date <- min(rain_hsm$dates, na.rm = TRUE)
max_date <- max(rain_hsm$dates, na.rm = TRUE)
extended_min_date <- min_date - 3600*2
extended_max_date <- max_date + 3600*2

filtered_comephore <- rain_com %>%
  filter(date >= extended_min_date & date <= extended_max_date)

# Create empty CSV with correct header
grid_df <- data.frame(
  time   = character(),
  station= character(),
  lon_Y  = numeric(),
  lat_Y  = numeric(),
  lon_X  = numeric(),
  lat_X  = numeric(),
  Y_obs  = numeric(),
  stringsAsFactors = FALSE
)
for (i in 1:27) grid_df[[paste0("X", i)]] <- numeric()

write.table(grid_df, output_file,
            row.names = FALSE, col.names = TRUE, sep = ";", append = FALSE)

# Remove rows where ALL gauges are NA
rain_hsm_nona <- rain_hsm[!apply(rain_hsm[, -which(names(rain_hsm) == "dates")],
                                 1, function(x) all(is.na(x))), ]

# Chunk writer
write_to_file <- function(data, file_path, append = FALSE) {
  write.table(
    data, file_path,
    row.names = FALSE, col.names = !append,
    sep = ";", append = append
  )
}
# =========================
# PARALLEL BUILD OF DOWNSCALING TABLE
# =========================


library(future)

cpus <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "1"))
if (is.na(cpus) || cpus < 1) cpus <- 1

# cpus <- max(1, cpus - 1)
plan(multisession, workers = cpus)

options(
  parallelly.maxWorkers.localhost = cpus,
  future.fork.enable = FALSE
)

# Create a stable station id column
if (!("station" %in% names(location_gauges))) {
  if ("codestation" %in% names(location_gauges)) {
    location_gauges$station <- location_gauges$codestation
  } else if ("Station" %in% names(location_gauges)) {
    location_gauges$station <- location_gauges$Station
  } else {
    stop("No station id column found in location_gauges.")
  }
}

# Precompute spatial neighborhood pixels for each station (once!)
radius_m <- 1500

pixel_in_radius_by_station <- lapply(seq_len(nrow(location_gauges)), function(i) {
  cx <- location_gauges$coord_x_px[i]
  cy <- location_gauges$coord_y_px[i]

  dpx <- geosphere::distHaversine(
    cbind(cx, cy),
    cbind(loc_px$Longitude, loc_px$Latitude)
  )
  pix <- loc_px$pixel_name[dpx <= radius_m]
  # keep exactly 9 pixels if you want a fixed 3x3 neighborhood:
  # pix <- pix[order(dpx[dpx <= radius_m])][1:min(9, sum(dpx <= radius_m))]
  pix
})
names(pixel_in_radius_by_station) <- location_gauges$station

# Precompute temporal indices in COMEPHORE for each OMSEV timestamp (once!)
# Sort COMEPHORE times
filtered_comephore <- filtered_comephore %>% arrange(date)
com_times <- filtered_comephore$date

# Colonnes stations dispo dans rain_hsm_nona (tout sauf dates)
rain_stations <- setdiff(names(rain_hsm_nona), "dates")

# On garde uniquement les stations présentes dans rain
location_gauges <- location_gauges %>%
  dplyr::filter(station %in% rain_stations)

stopifnot(nrow(location_gauges) > 0)

# OMSEV times
hsm_times <- rain_hsm_nona$dates
hsm_times <- as.POSIXct(hsm_times, tz="GMT")

# Function: indices for ±1.5h using binary search (fast)
get_time_window_idx <- function(t_date, com_times, half_window_secs = 5400) {
  t_start <- t_date - half_window_secs
  t_end   <- t_date + half_window_secs
  # findInterval gives position of last <= value
  left  <- findInterval(t_start, com_times) + 1
  right <- findInterval(t_end, com_times)
  if (left > right) integer(0) else seq.int(left, right)
}

time_indices_list <- lapply(hsm_times, get_time_window_idx, com_times = com_times)

# Parallel over timestamps t
process_one_time <- function(t_idx) {
  print(paste0("Processing time index ", t_idx, " / ", nrow(rain_hsm_nona)))
  t_row  <- rain_hsm_nona[t_idx, ]
  t_date <- as.POSIXct(t_row$dates, tz="GMT")

  idx <- time_indices_list[[t_idx]]
  if (length(idx) == 0) return(NULL)

  cube_obs <- filtered_comephore[idx, , drop = FALSE]
  cube_mat <- as.matrix(cube_obs[, -1, drop = FALSE])
  colnames(cube_mat) <- colnames(cube_obs)[-1]

  out_list <- vector("list", length = nrow(location_gauges))
  k <- 0

  for (i in seq_len(nrow(location_gauges))) {
    st <- location_gauges$station[i]
    y  <- t_row[[st]]

    if (is.na(y)) next

    pix <- pixel_in_radius_by_station[[st]]
    if (length(pix) == 0) next

    pix_cols <- match(pix, colnames(cube_mat))
    pix_cols <- pix_cols[!is.na(pix_cols)]

    if (length(pix_cols) == 0) next

    X_obs <- as.vector(cube_mat[, pix_cols, drop=FALSE])

    # enforce exactly 27 values
    if (length(X_obs) < 27) X_obs <- c(X_obs, rep(NA_real_, 27 - length(X_obs)))
    if (length(X_obs) > 27) X_obs <- X_obs[1:27]

    k <- k + 1
    out_list[[k]] <- data.frame(
      time    = t_date,
      station = st,
      lon_Y   = location_gauges$Longitude[i],
      lat_Y   = location_gauges$Latitude[i],
      lon_X   = location_gauges$coord_x_px[i],
      lat_X   = location_gauges$coord_y_px[i],
      Y_obs   = as.numeric(y),
      setNames(as.list(X_obs), paste0("X", 1:27)),
      stringsAsFactors = FALSE
    )
  }

  out_list <- Filter(Negate(is.null), out_list)
  if (length(out_list) == 0) return(NULL)
  dplyr::bind_rows(out_list)
}

message("SLURM_CPUS_PER_TASK = ", Sys.getenv("SLURM_CPUS_PER_TASK"))
print(future::nbrOfWorkers())

# Run in parallel (returns list of data.frames)
res_list <- future.apply::future_lapply(seq_len(nrow(rain_hsm_nona)), process_one_time)

# Bind + write once (much faster than writing chunks from each worker)
res_df <- dplyr::bind_rows(res_list)

# Write CSV
dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)
write.table(res_df, output_file, row.names = FALSE, col.names = TRUE, sep=";")
message("Parallel downscaling_table.csv written to: ", output_file)

# Optional: free workers
future::plan(future::sequential)
