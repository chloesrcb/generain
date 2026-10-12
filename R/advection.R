# =============================================================================
# FAST PARALLEL ADVECTION FUNCTIONS
# =============================================================================
#
# This version avoids lubridate objects inside parallel workers.
# Time offsets are handled as numeric seconds, which is more robust with future.
# =============================================================================


# =============================================================================
# HELPER: CONVERT q TO THE VALUE USED IN FILENAMES
# =============================================================================

format_q_for_filename <- function(q) {

  if (q <= 1) {
    return(round(q * 100))
  }

  return(round(q))
}


# =============================================================================
# HELPER: ROUND POSIXct TO THE NEXT HOUR
# =============================================================================

ceiling_hour_base <- function(x) {

  x_num <- as.numeric(as.POSIXct(x, tz = "UTC"))

  rounded_num <- ceiling(x_num / 3600) * 3600

  as.POSIXct(
    rounded_num,
    origin = "1970-01-01",
    tz = "UTC"
  )
}


# =============================================================================
# COMPUTE RAINFALL BARYCENTERS FOR ALL TIME STEPS
# =============================================================================

compute_all_barycenters <- function(
    data,
    rain_threshold = 0
) {

  required_columns <- c(
    "dates",
    "rain",
    "Longitude",
    "Latitude"
  )

  missing_columns <- setdiff(
    required_columns,
    names(data)
  )

  if (length(missing_columns) > 0) {
    stop(
      "Missing columns in rainfall data: ",
      paste(missing_columns, collapse = ", ")
    )
  }


  centers <- data |>
    dplyr::group_by(dates) |>
    dplyr::summarise(

      max_rain = if (
        all(is.na(rain))
      ) {
        NA_real_
      } else {
        max(rain, na.rm = TRUE)
      },

      total_rain = sum(
        rain,
        na.rm = TRUE
      ),

      x_num = sum(
        rain * Longitude,
        na.rm = TRUE
      ),

      y_num = sum(
        rain * Latitude,
        na.rm = TRUE
      ),

      .groups = "drop"
    ) |>
    dplyr::filter(
      !is.na(max_rain),
      max_rain >= rain_threshold,
      total_rain > 0
    ) |>
    dplyr::mutate(
      x_c = x_num / total_rain,
      y_c = y_num / total_rain
    ) |>
    dplyr::filter(
      is.finite(x_c),
      is.finite(y_c)
    ) |>
    dplyr::select(
      dates,
      x_c,
      y_c
    ) |>
    dplyr::arrange(dates)


  # Force dates to POSIXct explicitly
  centers$dates <- as.POSIXct(
    centers$dates,
    origin = "1970-01-01",
    tz = "UTC"
  )

  return(centers)
}


# =============================================================================
# COMPUTE ADVECTION FOR ONE EPISODE
# =============================================================================

compute_episode_advection <- function(
    t0,
    pre_seconds,
    post_seconds,
    centers,
    meters_per_degree,
    cos_lat
) {

  # Force t0 to a scalar POSIXct object
  t0 <- as.POSIXct(
    as.numeric(t0),
    origin = "1970-01-01",
    tz = "UTC"
  )


  # ---------------------------------------------------------------------------
  # Define the episode time window using numeric seconds
  # ---------------------------------------------------------------------------

  t_start <- t0 - as.numeric(pre_seconds)
  t_end <- t0 + as.numeric(post_seconds)


  # ---------------------------------------------------------------------------
  # Keep only barycenters inside the episode window
  # ---------------------------------------------------------------------------

  centers_ep <- centers[
    centers$dates >= t_start &
      centers$dates <= t_end,
    ,
    drop = FALSE
  ]


  # At least two valid time steps are required
  if (nrow(centers_ep) < 2) {
    return(
      c(
        mean_dx_kmh = 0,
        mean_dy_kmh = 0
      )
    )
  }


  # ---------------------------------------------------------------------------
  # Compute episode duration directly from numeric timestamps
  # ---------------------------------------------------------------------------

  duration_hours <- (
    as.numeric(centers_ep$dates[nrow(centers_ep)]) -
      as.numeric(centers_ep$dates[1])
  ) / 3600


  if (
    !is.finite(duration_hours) ||
      duration_hours <= 0
  ) {
    return(
      c(
        mean_dx_kmh = 0,
        mean_dy_kmh = 0
      )
    )
  }


  # ---------------------------------------------------------------------------
  # Compute total barycenter displacement
  # ---------------------------------------------------------------------------

  dx_deg <- (
    centers_ep$x_c[nrow(centers_ep)] -
      centers_ep$x_c[1]
  )

  dy_deg <- (
    centers_ep$y_c[nrow(centers_ep)] -
      centers_ep$y_c[1]
  )


  # ---------------------------------------------------------------------------
  # Convert displacement to km/h
  # ---------------------------------------------------------------------------

  dx_kmh <- (
    dx_deg *
      meters_per_degree *
      cos_lat
  ) / 1000 / duration_hours


  dy_kmh <- (
    dy_deg *
      meters_per_degree
  ) / 1000 / duration_hours


  return(
    c(
      mean_dx_kmh = dx_kmh,
      mean_dy_kmh = dy_kmh
    )
  )
}


# =============================================================================
# PREPARE OMSEV AND COMEPHORE DATA
# =============================================================================

prepare_advection_data <- function(
    data_folder = "./data/"
) {

  message("Loading OMSEV and COMEPHORE data...")


  required_packages <- c(
    "dplyr",
    "tidyr",
    "readr"
  )

  missing_packages <- required_packages[
    !vapply(
      required_packages,
      requireNamespace,
      FUN.VALUE = logical(1),
      quietly = TRUE
    )
  ]

  if (length(missing_packages) > 0) {
    stop(
      "Missing R packages: ",
      paste(missing_packages, collapse = ", ")
    )
  }


  # ---------------------------------------------------------------------------
  # OMSEV rainfall
  # ---------------------------------------------------------------------------

  filename_omsev <- file.path(
    data_folder,
    "omsev",
    "omsev_5min",
    "rain_mtp_5min_2019_2025.csv"
  )

  if (!file.exists(filename_omsev)) {
    stop(
      "OMSEV rainfall file not found: ",
      filename_omsev
    )
  }

  rain <- utils::read.csv(
    filename_omsev,
    check.names = FALSE
  )

  if (!"dates" %in% names(rain)) {
    stop(
      "Column 'dates' is missing from OMSEV rainfall file."
    )
  }

  rain$dates <- as.POSIXct(
    rain$dates,
    tz = "UTC"
  )


  # ---------------------------------------------------------------------------
  # OMSEV gauge coordinates
  # ---------------------------------------------------------------------------

  filename_loc_gauges <- file.path(
    data_folder,
    "omsev",
    "loc_rain_gauges.csv"
  )

  if (!file.exists(filename_loc_gauges)) {
    stop(
      "OMSEV gauge location file not found: ",
      filename_loc_gauges
    )
  }

  location_gauges <- utils::read.csv(
    filename_loc_gauges,
    check.names = FALSE
  )


  required_gauge_columns <- c(
    "Station",
    "Longitude",
    "Latitude"
  )

  missing_gauge_columns <- setdiff(
    required_gauge_columns,
    names(location_gauges)
  )

  if (length(missing_gauge_columns) > 0) {
    stop(
      "Missing columns in OMSEV gauge location file: ",
      paste(missing_gauge_columns, collapse = ", ")
    )
  }


  # ---------------------------------------------------------------------------
  # Convert OMSEV data to long format
  # ---------------------------------------------------------------------------

  rain_long <- rain |>
    tidyr::pivot_longer(
      cols = -dates,
      names_to = "Station",
      values_to = "rain"
    ) |>
    dplyr::left_join(
      location_gauges,
      by = "Station"
    )


  positive_rain <- rain_long$rain[
    is.finite(rain_long$rain) &
      rain_long$rain > 0
  ]

  if (length(positive_rain) == 0) {
    stop(
      "No positive rainfall value found in OMSEV data."
    )
  }

  rain_threshold <- min(
    positive_rain,
    na.rm = TRUE
  )


  message("Precomputing OMSEV barycenters...")

  centers_omsev <- compute_all_barycenters(
    data = rain_long,
    rain_threshold = rain_threshold
  )


  rm(
    rain,
    rain_long
  )

  gc(verbose = FALSE)


  # ---------------------------------------------------------------------------
  # COMEPHORE rainfall
  # ---------------------------------------------------------------------------

  filename_com <- file.path(
    data_folder,
    "comephore",
    "rebuild_clean",
    "comephore_2008_2025_within5km.csv"
  )

  if (!file.exists(filename_com)) {
    stop(
      "COMEPHORE rainfall file not found: ",
      filename_com
    )
  }

  comephore_raw <- utils::read.csv(
    filename_com,
    check.names = FALSE
  )


  # Detect date column
  if ("datetime" %in% names(comephore_raw)) {

    comephore_raw$datetime <- as.POSIXct(
      comephore_raw$datetime,
      tz = "UTC"
    )

    comephore_raw <- comephore_raw |>
      dplyr::filter(
        datetime >= as.POSIXct(
          "2019-01-01",
          tz = "UTC"
        )
      ) |>
      dplyr::rename(
        date = datetime
      )

  } else if ("date" %in% names(comephore_raw)) {

    comephore_raw$date <- as.POSIXct(
      comephore_raw$date,
      tz = "UTC"
    )

    comephore_raw <- comephore_raw |>
      dplyr::filter(
        date >= as.POSIXct(
          "2019-01-01",
          tz = "UTC"
        )
      )

  } else {

    stop(
      "No 'datetime' or 'date' column found in COMEPHORE file."
    )
  }


  # ---------------------------------------------------------------------------
  # COMEPHORE pixel coordinates
  # ---------------------------------------------------------------------------

  filename_loc_pixels <- file.path(
    data_folder,
    "comephore",
    "rebuild_clean",
    "coords_pixels_within5km.csv"
  )

  if (!file.exists(filename_loc_pixels)) {
    stop(
      "COMEPHORE pixel location file not found: ",
      filename_loc_pixels
    )
  }

  loc_px <- utils::read.csv(
    filename_loc_pixels,
    check.names = FALSE
  )


  required_pixel_columns <- c(
    "pixel_name",
    "Longitude",
    "Latitude"
  )

  missing_pixel_columns <- setdiff(
    required_pixel_columns,
    names(loc_px)
  )

  if (length(missing_pixel_columns) > 0) {
    stop(
      "Missing columns in COMEPHORE pixel location file: ",
      paste(missing_pixel_columns, collapse = ", ")
    )
  }


  loc_px <- loc_px |>
    dplyr::select(
      pixel_name,
      Longitude,
      Latitude
    )


  # ---------------------------------------------------------------------------
  # Convert COMEPHORE data to long format
  # ---------------------------------------------------------------------------

  com_long <- comephore_raw |>
    tidyr::pivot_longer(
      cols = -date,
      names_to = "Station",
      values_to = "rain"
    ) |>
    dplyr::rename(
      dates = date
    ) |>
    dplyr::left_join(
      loc_px,
      by = c(
        "Station" = "pixel_name"
      )
    )


  message("Precomputing COMEPHORE barycenters...")

  centers_com <- compute_all_barycenters(
    data = com_long,
    rain_threshold = 0
  )


  rm(
    comephore_raw,
    com_long
  )

  gc(verbose = FALSE)


  # ---------------------------------------------------------------------------
  # Geographic conversion constants
  # ---------------------------------------------------------------------------

  meters_per_degree <- 111320


  cos_lat_omsev <- cos(
    mean(
      location_gauges$Latitude,
      na.rm = TRUE
    ) *
      pi / 180
  )


  cos_lat_com <- cos(
    mean(
      loc_px$Latitude,
      na.rm = TRUE
    ) *
      pi / 180
  )


  message(
    "Advection data prepared: ",
    nrow(centers_omsev),
    " OMSEV time steps and ",
    nrow(centers_com),
    " COMEPHORE time steps."
  )


  return(
    list(
      centers_omsev = centers_omsev,
      centers_com = centers_com,
      meters_per_degree = meters_per_degree,
      cos_lat_omsev = cos_lat_omsev,
      cos_lat_com = cos_lat_com
    )
  )
}


# =============================================================================
# COMPUTE ADVECTION FOR ONE CONFIGURATION
# =============================================================================

compute_advection <- function(
    q,
    delta,
    dmin,
    adv_data,
    data_folder = "./data/",
    beta_lag = 0,
    save_output = TRUE
) {

  q_file <- format_q_for_filename(q)


  message(
    "Computing advection for q = ",
    q_file,
    ", delta = ",
    delta,
    ", dmin = ",
    dmin
  )


  # ---------------------------------------------------------------------------
  # Optional beta suffix
  # ---------------------------------------------------------------------------

  if (beta_lag == 0) {
    end_file <- ""
  } else {
    end_file <- paste0(
      "_beta",
      beta_lag
    )
  }


  # ---------------------------------------------------------------------------
  # Load episode start times
  # ---------------------------------------------------------------------------

  filename_t0 <- sprintf(
    "t0_5min_episodes_q%d_delta%d_dmin%d%s.csv",
    q_file,
    delta,
    dmin,
    end_file
  )


  filepath_t0 <- file.path(
    data_folder,
    "omsev",
    "random",
    filename_t0
  )


  if (!file.exists(filepath_t0)) {
    stop(
      "Episode file does not exist: ",
      filepath_t0
    )
  }


  t0_df <- readr::read_csv(
    filepath_t0,
    show_col_types = FALSE
  )


  if (ncol(t0_df) == 0) {
    stop(
      "Episode file is empty: ",
      filepath_t0
    )
  }


  names(t0_df)[1] <- "t0"


  # ---------------------------------------------------------------------------
  # Convert t0 to numeric UNIX timestamps
  # This avoids POSIXlt/POSIXct issues inside parallel workers.
  # ---------------------------------------------------------------------------

  t0_character <- as.character(
    t0_df$t0
  )

  date_only <- grepl(
    "^\\d{4}-\\d{2}-\\d{2}$",
    t0_character
  )

  t0_character[date_only] <- paste0(
    t0_character[date_only],
    " 00:00:00"
  )

  t0_posix <- as.POSIXct(
    t0_character,
    tz = "UTC"
  )


  if (anyNA(t0_posix)) {
    bad_rows <- which(
      is.na(t0_posix)
    )

    stop(
      "Could not parse t0 dates at row(s): ",
      paste(head(bad_rows, 20), collapse = ", ")
    )
  }


  # Numeric timestamps are safest for future workers
  t0_num <- as.numeric(
    t0_posix
  )


  # ---------------------------------------------------------------------------
  # Time windows expressed in seconds
  # ---------------------------------------------------------------------------

  pre_omsev_seconds <- 30 * 60
  post_omsev_seconds <- 30 * 60

  pre_com_seconds <- 60 * 60
  post_com_seconds <- 60 * 60


  # ---------------------------------------------------------------------------
  # OMSEV advection
  # ---------------------------------------------------------------------------

  message(
    "Computing OMSEV advection for ",
    length(t0_num),
    " episodes..."
  )


  omsev_list <- future.apply::future_lapply(
    X = seq_along(t0_num),
    FUN = function(i) {

      t0_i <- as.POSIXct(
        t0_num[i],
        origin = "1970-01-01",
        tz = "UTC"
      )


      result <- compute_episode_advection(
        t0 = t0_i,
        pre_seconds = pre_omsev_seconds,
        post_seconds = post_omsev_seconds,
        centers = adv_data$centers_omsev,
        meters_per_degree = adv_data$meters_per_degree,
        cos_lat = adv_data$cos_lat_omsev
      )


      data.frame(
        episode = i,
        t0_omsev = t0_i,
        mean_dx_kmh_omsev = unname(
          result["mean_dx_kmh"]
        ),
        mean_dy_kmh_omsev = unname(
          result["mean_dy_kmh"]
        )
      )
    },
    future.seed = TRUE,
    future.scheduling = 2
  )


  episode_advect_omsev <- dplyr::bind_rows(
    omsev_list
  )


  # ---------------------------------------------------------------------------
  # COMEPHORE advection
  # ---------------------------------------------------------------------------

  message(
    "Computing COMEPHORE advection for ",
    length(t0_num),
    " episodes..."
  )


  com_list <- future.apply::future_lapply(
    X = seq_along(t0_num),
    FUN = function(i) {

      # Round to the next hour using base R arithmetic
      t0_i <- as.POSIXct(
        ceiling(
          t0_num[i] / 3600
        ) * 3600,
        origin = "1970-01-01",
        tz = "UTC"
      )


      result <- compute_episode_advection(
        t0 = t0_i,
        pre_seconds = pre_com_seconds,
        post_seconds = post_com_seconds,
        centers = adv_data$centers_com,
        meters_per_degree = adv_data$meters_per_degree,
        cos_lat = adv_data$cos_lat_com
      )


      data.frame(
        episode = i,
        t0_comephore = t0_i,
        mean_dx_kmh_comephore = unname(
          result["mean_dx_kmh"]
        ),
        mean_dy_kmh_comephore = unname(
          result["mean_dy_kmh"]
        )
      )
    },
    future.seed = TRUE,
    future.scheduling = 2
  )


  episode_advect_com <- dplyr::bind_rows(
    com_list
  )


  # ---------------------------------------------------------------------------
  # Combine OMSEV and COMEPHORE
  # ---------------------------------------------------------------------------

  adv_df <- dplyr::left_join(
    episode_advect_omsev,
    episode_advect_com,
    by = "episode"
  )


  # Use OMSEV only when COMEPHORE returns a zero vector
  use_omsev <- (
    adv_df$mean_dx_kmh_comephore == 0 &
      adv_df$mean_dy_kmh_comephore == 0
  )


  adv_df$vx_final <- ifelse(
    use_omsev,
    adv_df$mean_dx_kmh_omsev,
    adv_df$mean_dx_kmh_comephore
  )


  adv_df$vy_final <- ifelse(
    use_omsev,
    adv_df$mean_dy_kmh_omsev,
    adv_df$mean_dy_kmh_comephore
  )


  adv_df$speed_kmh <- sqrt(
    adv_df$vx_final^2 +
      adv_df$vy_final^2
  )


  # ---------------------------------------------------------------------------
  # Save output
  # ---------------------------------------------------------------------------

  if (save_output) {

    out_dir <- file.path(
      data_folder,
      "omsev",
      "adv_estim",
      "2025_all",
      "combined_comephore_omsev",
      "random"
    )


    if (!dir.exists(out_dir)) {
      dir.create(
        out_dir,
        recursive = TRUE
      )
    }


    filename_out <- sprintf(
      "episode_advection_q%d_delta%d_dmin%d%s.csv",
      q_file,
      delta,
      dmin,
      end_file
    )


    output_path <- file.path(
      out_dir,
      filename_out
    )


    utils::write.csv(
      adv_df,
      output_path,
      row.names = FALSE
    )


    message(
      "Advection saved to: ",
      output_path
    )
  }


  return(adv_df)
}


# =============================================================================
# LOAD EXISTING ADVECTION OR COMPUTE IT
# =============================================================================

get_advection <- function(
    q,
    delta,
    dmin,
    data_folder = "./data/",
    adv_data = NULL,
    beta_lag = 0,
    force = FALSE
) {

  q_file <- format_q_for_filename(q)


  if (beta_lag == 0) {
    end_file <- ""
  } else {
    end_file <- paste0(
      "_beta",
      beta_lag
    )
  }


  adv_dir <- file.path(
    data_folder,
    "omsev",
    "adv_estim",
    "2025_all",
    "combined_comephore_omsev",
    "random"
  )


  adv_filename <- sprintf(
    "episode_advection_q%d_delta%d_dmin%d%s.csv",
    q_file,
    delta,
    dmin,
    end_file
  )


  adv_path <- file.path(
    adv_dir,
    adv_filename
  )


  # Load existing file unless recomputation is requested
  if (
    file.exists(adv_path) &&
      !force
  ) {

    message(
      "Advection already computed. Loading: ",
      adv_path
    )


    return(
      utils::read.csv(
        adv_path,
        sep = ","
      )
    )
  }


  message(
    "Advection must be computed."
  )


  # Prepare data only if needed
  if (is.null(adv_data)) {
    adv_data <- prepare_advection_data(
      data_folder = data_folder
    )
  }


  return(
    compute_advection(
      q = q,
      delta = delta,
      dmin = dmin,
      adv_data = adv_data,
      data_folder = data_folder,
      beta_lag = beta_lag,
      save_output = TRUE
    )
  )
}