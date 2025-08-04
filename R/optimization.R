# EXTREME EPISODES -------------------------------------------------------------

#' get_spatiotemp_excess function
#' 
#' This function computes the spatio-temporal excesses based on a quantile or
#' threshold. It returns a list of sites, times, and thresholds where the
#' excesses occur.
#' 
#' @param data The rainfall data as a matrix.
#' @param quantile The quantile value to compute the threshold.
#' @param threshold A fixed threshold value or vector.
#' 
#' @return A list containing three vectors: list_s (sites), list_t (times),
#' 
#' @export 
get_spatiotemp_excess <- function(data, quantile = NULL,
                                  threshold = NULL, remove_zeros = FALSE) {
  data <- as.matrix(data)
  site_names <- colnames(data)

  if (!is.null(quantile) && !is.null(threshold)) stop("Specify either 'quantile' or 'threshold', not both.")
  if (is.null(quantile) && is.null(threshold)) stop("Specify 'quantile' or 'threshold'.")

  if (!is.null(quantile)) {

    if (remove_zeros) {
        thresholds_by_site <- apply(data, 2, function(col) {
        col <- col[!is.na(col) & col > 0]  # Retire NA et les 0
        if (length(col) < 50) {
          return(NA)  # Trop peu de données non-nulles pour estimer un seuil fiable
        }
        quantile(col, probs = q, na.rm = TRUE)
      }) 
    }
    else {
      thresholds_by_site <- apply(data, 2, function(col) quantile(col, probs = quantile, na.rm = TRUE))
    }
  } else {
    if (length(threshold) == 1) {
      thresholds_by_site <- rep(threshold, ncol(data))
      names(thresholds_by_site) <- site_names
    } else if (length(threshold) == ncol(data)) {
      thresholds_by_site <- threshold
      names(thresholds_by_site) <- site_names
    } else {
      stop("Invalid threshold input.")
    }
  }

  list_s <- list()
  list_t <- list()
  list_u <- list()

  for (i in seq_len(ncol(data))) {
    s0 <- site_names[i]
    u <- thresholds_by_site[s0]
    t0 <- which(data[, i] > thresholds_by_site[s0])
    if (length(t0) > 0) {
      list_s <- c(list_s, rep(s0, length(t0)))
      list_t <- c(list_t, t0)
      list_u <- c(list_u, rep(u, length(t0)))
    }
  }
  return(list(list_s = list_s, list_t = list_t, list_u = list_u))
}


#' get_s0t0_pairs function
#' 
#' This function retrieves pairs of sites and times (s0, t0) based on spatio-temporal
#' excesses. It ensures that selected pairs are sufficiently distant in both space and time.
#' 
#' @param sites_coords The coordinates of the sites.
#' @param data The rainfall data as a matrix.
#' @param min_spatial_dist The minimum spatial distance between two episodes.
#' @param episode_size The temporal window size.
#' @param set_st_excess A list containing the spatio-temporal excesses with
#'                      list_s, list_t, and list_u.
#' @param n_max_episodes The maximum number of episodes to select. Default is 10000.
#' @param latlon A boolean value indicating if the coordinates are in latitude/longitude.
#' @return A data.table with selected pairs of sites and times (s0, t0) and their
#'         corresponding excess values u_s0.
#' 
#' @import data.table
#' @import geosphere
#' @export
get_s0t0_pairs <- function(sites_coords, data, min_spatial_dist,
                           episode_size, set_st_excess,
                           n_max_episodes = 10000, latlon = FALSE){

  data <- as.matrix(data)
  site_names <- colnames(data)

  list_s <- set_st_excess$list_s
  list_t <- set_st_excess$list_t
  list_u <- set_st_excess$list_u

  # Compute distance matrix
  if (latlon) {
    # dist_matrix <- as.matrix(distm(sites_coords[site_names, ], fun = distHaversine)) / 1000
    dist_matrix <- get_dist_mat(sites_coords, latlon = TRUE)
  } else {
    dist_matrix <- get_dist_mat(sites_coords, latlon = FALSE)
    dist_matrix <- ceiling(round(dist_matrix, 1) * 1000) / 1000
  }
  rownames(dist_matrix) <- site_names
  colnames(dist_matrix) <- site_names

  # Initialization
  selected_points <- data.table(s0 = character(), t0 = integer(), 
                                t0_date = character(), u_s0 = numeric())
  nb_episode <- 0

  # Sort excesses by time for consistent processing
  excess_dt <- data.table(
    s = unlist(list_s),
    t = unlist(list_t),
    t_dates = names(list_t),
    u = unlist(list_u)
  )
  setorder(excess_dt, t)

  for (i in seq_len(nrow(excess_dt))) {
    s0 <- excess_dt$s[i]
    t0 <- excess_dt$t[i]
    t0_date <- excess_dt$t_dates[i]
    u0 <- excess_dt$u[i]

    # Check if new point is valid in space-time relative to all previously selected points
    is_valid <- TRUE
    for (j in seq_len(nrow(selected_points))) {
      s_prev <- selected_points$s0[j]
      t_prev <- selected_points$t0[j]

      spatial_dist <- dist_matrix[s0, s_prev]
      temporal_dist <- abs(t0 - t_prev)

      if (spatial_dist < min_spatial_dist && temporal_dist < episode_size) {
        is_valid <- FALSE
        break
      }
    }

    if (is_valid) {
      selected_points <- rbind(selected_points, data.table(s0 = s0, t0 = t0, 
                                                t0_date = t0_date, u_s0 = u0))
      nb_episode <- nb_episode + 1
      if (nb_episode >= n_max_episodes) break
    }
  }

  return(selected_points)
}


#' select_extreme_episodes function (Optimized)
#'
#' This function selects extreme episodes based on a quantile or threshold,
#' with spatial and temporal exclusion. Optimized for performance.
#'
#' @param sites_coords The coordinates of the sites.
#' @param data The rainfall dataframe.
#' @param quantile The quantile value.
#' @param min_spatial_dist The minimum spatial distance between two episodes.
#' @param episode_size The temporal window size.
#' @param n_max_episodes Max number of episodes to select. Default 10000.
#' @param time_ext Temporal window extension. Default 0.
#' @param latlon TRUE if coordinates are lat/lon. Default TRUE.
#' @param spatial_window_radius Optional radius (km) for candidate sites.
#' @param central_site Optional name to center spatial window.
#' @param threshold Optional fixed threshold(s).
#' @return data.table with selected episodes: s0, t0, u_s0.
#'
#' @import data.table
#' @import geosphere
#' @export
select_extreme_episodes <- function(sites_coords, data, min_spatial_dist,
                                    episode_size, n_max_episodes = 10000,
                                    time_ext = 0, latlon = TRUE,
                                    quantile = NULL, threshold = NULL,
                                    spatial_window_radius = NULL,
                                    central_site = NULL) {
  data <- as.matrix(data)
  site_names <- colnames(data)

  # --- Threshold computation ---
  if (!is.null(quantile) && !is.null(threshold)) stop("Specify either 'quantile' or 'threshold', not both.")
  if (is.null(quantile) && is.null(threshold)) stop("Specify 'quantile' or 'threshold'.")

  if (!is.null(quantile)) {
    thresholds_by_site <- apply(data, 2, function(col) quantile(col, probs = quantile, na.rm = TRUE))
  } else {
    if (length(threshold) == 1) {
      thresholds_by_site <- rep(threshold, ncol(data))
      names(thresholds_by_site) <- site_names
    } else if (length(threshold) == ncol(data)) {
      thresholds_by_site <- threshold
      names(thresholds_by_site) <- site_names
    } else {
      stop("Invalid threshold input.")
    }
  }

  # --- Distance matrix ---
  if (latlon) {
    dist_matrix <- as.matrix(distm(sites_coords[site_names, ], fun = distHaversine)) / 1000
  } else {
    dist_matrix <- get_dist_mat(sites_coords, latlon = FALSE)
  }
  rownames(dist_matrix) <- site_names
  colnames(dist_matrix) <- site_names

  # --- Spatial window site filtering ---
  valid_site_mask <- rep(TRUE, length(site_names))
  if (!is.null(spatial_window_radius)) {
    if (!is.null(central_site)) {
      valid_sites <- names(which(dist_matrix[, central_site] <= spatial_window_radius))
    } else {
      valid_sites <- names(which(apply(dist_matrix, 2, function(col) any(col <= spatial_window_radius))))
    }
    valid_site_mask <- site_names %in% valid_sites
  }

  # --- Initialization ---
  selected_points <- data.table(s0 = character(), t0 = integer(), u_s0 = numeric())
  invalid_time_mask <- matrix(FALSE, nrow = nrow(data), ncol = ncol(data))
  nb_episode <- 0

  while (nb_episode < n_max_episodes) {
    exceed_mask <- sweep(data, 2, thresholds_by_site, FUN = ">")
    exceed_mask[invalid_time_mask] <- FALSE
    exceed_indices <- which(exceed_mask, arr.ind = TRUE)

    # Apply spatial mask
    exceed_indices <- exceed_indices[exceed_indices[, 2] %in% which(valid_site_mask), , drop = FALSE]
    if (nrow(exceed_indices) == 0) break

    # Sort by exceedance value (descending), keep top 200
    exceed_values <- data[exceed_indices]
    ord <- order(exceed_values, decreasing = TRUE)
    exceed_indices <- exceed_indices[ord[1:min(200, length(ord))], , drop = FALSE]

    selected <- FALSE

    for (i in seq_len(nrow(exceed_indices))) {
      t0 <- exceed_indices[i, 1]
      s0_idx <- exceed_indices[i, 2]
      s0 <- site_names[s0_idx]
      u_s0 <- thresholds_by_site[s0]

      if (nrow(selected_points) == 0) {
        selected_points <- rbind(selected_points, data.table(s0 = s0, t0 = t0, u_s0 = u_s0))
        nb_episode <- nb_episode + 1
        selected <- TRUE
        break
      }

      selected_sites <- selected_points$s0
      time_diffs <- abs(selected_points$t0 - t0)
      spatial_dists <- dist_matrix[selected_sites, s0]

      if (all(spatial_dists >= min_spatial_dist | time_diffs >= episode_size)) {
        selected_points <- rbind(selected_points, data.table(s0 = s0, t0 = t0, u_s0 = u_s0))
        nb_episode <- nb_episode + 1
        selected <- TRUE
        break
      } else {
        invalid_time_mask[t0, s0_idx] <- TRUE
      }
    }

    if (selected) {
      t_inf <- max(1, t0 - time_ext)
      t_sup <- min(nrow(data), t0 + episode_size + time_ext - 1)
      invalid_time_mask[t_inf:t_sup, s0_idx] <- TRUE
    } else {
      break
    }
  }

  return(selected_points)
}


#' select_max_extreme_episodes function
#'
#' This function selects extreme episodes based on the given data and a 
#' quantile. A spatio-temporal specific neighborhood is considered to avoid
#' selecting nearby episodes. Each episode is selected according to a 
#' r-Pareto process with r(X) = max(X(s, t)). Each episode will have a size of
#' 2 * delta.
#'
#' @param sites_coords The coordinates of the sites.
#' @param data The rainfall dataframe.
#' @param quantile The quantile value.
#' @param min_spatial_dist The minimum spatial distance between two episodes.
#' @param delta The temporal window size.
#' @param n_max_episodes The maximum number of episodes to select. Default is
#'                     10000.
#' @param time_ext The temporal window extension. Default is 0.
#' @param latlon A boolean value to indicate if the locations are in latitude
#'            and longitude. Default is TRUE.
#' @return A list of selected points (s0, t0) and the corresponding threshold 
#'         value u_s0 for each episode.
#'
#' @import data.table
#' @import geosphere
#'
#' @export
select_max_extreme_episodes <- function(sites_coords, data, quantile,
                                    min_spatial_dist, delta,
                                    n_max_episodes = 10000, time_ext = 0,
                                    latlon = TRUE) {
  # Convert data to a matrix for fast access
  data <- as.matrix(data)
  # n_sites <- ncol(data)

  # data_unif <- apply(data, 2, function(col) rank(col, na.last = "keep") /
  #                                               (length(col) + 1))
  threshold <- quantile(data, probs = quantile, na.rm = TRUE)
  # threshold <- 1
  # Extract site names
  site_names <- colnames(data)

  # Compute distance matrix (ensure named rows/columns)
  if (latlon){
    dist_matrix <- as.matrix(distm(sites_coords[site_names, ],
                                        fun = distHaversine)) / 1000
  } else {
    dist_matrix <- as.matrix(dist(sites_coords[site_names, ])) / 1000
  }
  rownames(dist_matrix) <- site_names
  colnames(dist_matrix) <- site_names

  # Get initial max values
  max_value <- max(na.omit(data))
  max_indices <- which(data == max_value, arr.ind = TRUE)

  # Store selected points and corresponding threshold
  selected_points <- data.table(s0 = character(), t0 = integer(),
                                u_s0 = numeric())

  # Logical mask for invalid times
  invalid_time_mask <- matrix(FALSE, nrow(data), ncol(data))

  nb_episode <- 0

  while (nb_episode < n_max_episodes && max_value > threshold) {
    best_candidate <- NULL

    for (idx in seq_len(nrow(max_indices))) {
      s0_candidate <- site_names[max_indices[idx, 2]]
      t0_candidate <- max_indices[idx, 1]
      u_s0 <- quantile(data[, s0_candidate], probs = quantile,
                              na.rm = TRUE)
      if (nrow(selected_points) == 0) { # First selection
        best_candidate <- data.table(s0 = s0_candidate, t0 = t0_candidate,
                                    u_s0 = u_s0)
        nb_episode <- nb_episode + 1
        break
      } else {
        # Convert site names to indices for distance lookup
        selected_sites <- selected_points$s0
        valid_indices <- which(site_names %in% selected_sites)

        # Compute distances only if there are previous selections
        if (length(valid_indices) > 0) {
          distances <- as.vector(dist_matrix[selected_sites, s0_candidate,
                                  drop = FALSE])
          time_differences <- abs(selected_points$t0 - t0_candidate)
          valid_pairs <- (distances >= min_spatial_dist) |
                                  (time_differences > 2 * delta + 1)
          # Validate the candidate
          if (all(valid_pairs)) {
            best_candidate <- data.table(s0 = s0_candidate, t0 = t0_candidate,
                                        u_s0 = u_s0)
            nb_episode <- nb_episode + 1
            # print(nb_episode)
            break
          } else {
            # Mark the candidate as invalid
            invalid_time_mask[t0_candidate, max_indices[idx, 2]] <- TRUE
          }
        }
      }
    }

    if (!is.null(best_candidate)) {
      # Store selected episode
      selected_points <- rbindlist(list(selected_points, best_candidate))

      # Remove nearby temporal data
      t_inf <- max(1, best_candidate$t0 - (delta) - time_ext)
      t_sup <- min(nrow(data), best_candidate$t0 + delta + time_ext)
      invalid_time_mask[t_inf:t_sup,
                        which(site_names == best_candidate$s0)] <- TRUE
    }

    # Update max selection by ignoring invalid positions
    masked_data <- data
    masked_data[invalid_time_mask] <- -Inf
    max_value <- max(masked_data)
    max_indices <- which(masked_data == max_value, arr.ind = TRUE)
  }

  return(selected_points)
}


#' check_intervals_overlap function
#'
#' This function checks if the time intervals overlap for all extreme episodes
#' associated with a given s0 value.
#'
#' @param s0_value The s0 value.
#' @param selected_points The selected points (s0, t0, u_s0).
#' @param delta The temporal window size.
#' @param beta The temporal window extension. Default is 0.
#'
#' @return TRUE if the intervals overlap, FALSE otherwise.
#'
#' @import data.table
#'
#' @export
check_intervals_overlap <- function(s0_value, selected_points, delta,
                                                                beta = 0) {
  # Filter for the given s0
  filtered_points <- selected_points[selected_points$s0 == s0_value, ]

  if (nrow(filtered_points) < 2) {
    return(FALSE)  # No overlap if only one interval exists
  }

  # Compute time intervals
  intervals <- data.table(
    t_inf = filtered_points$t0 - delta - beta,
    t_sup = filtered_points$t0 + delta + beta
  )

  setorder(intervals) # Ensure the data is sorted by time

  # Check for overlap
  for (i in 2:nrow(intervals)) {
    if (intervals$t_inf[i] <= intervals$t_sup[i - 1]) {
      return(TRUE)  # Overlap found
    }
  }

  return(FALSE)  # No overlap found
}


#' get_extreme_episodes function
#'
#' This function extracts the extreme episodes based on the selected points
#' and the data according to r-Pareto process definition with r(X) = X(s0, t0).
#'
#' @param selected_points The selected points (s0, t0).
#' @param data The rainfall data.
#' @param episode_size The temporal window size.
#' @param beta The temporal window extension. Default is 0.
#' @param unif A boolean value to indicate if we want to get the uniformized
#'            data episodes. Default is FALSE.
#'
#' @return The list of extreme episodes of size 2 * delta each.
#'
#' @export
get_extreme_episodes <- function(selected_points, data, episode_size, beta = 0,
                                unif = FALSE) {
  if (unif) { # Uniformize the data
    data <- apply(data, 2, function(col) rank(col, na.last = "keep") /
                                              (length(col) + 1))
  }

  # Store valid episodes and valid indices
  episodes <- list()
  valid_indices <- c()

  for (i in 1:nrow(selected_points)) {
    t0 <- selected_points$t0[i]
    # t_inf <- t0 - (delta) - beta
    # t_sup <- t0 + delta + 2 * beta
    t_inf <- max(1, t0 - beta)
    t_sup <- min(nrow(data), t0 + episode_size - 1 + beta)
    # print(i)
    # Check that the episode is the correct size (episode_size)
    episode_size_test <- t_sup - t_inf + 1
    if (episode_size_test == episode_size + 2 * beta) {
      episode <- data[t_inf:t_sup, , drop = FALSE] # Get the episode
      episodes <- append(episodes, list(episode))
      valid_indices <- c(valid_indices, i) # Mark this index as valid
    }
  }

  # Subset selected_points to keep only valid points
  new_selected_points <- selected_points[valid_indices, , drop = FALSE]

  # Print message if selected_points has changed
  if (nrow(new_selected_points) < nrow(selected_points)) {
    removed_count <- nrow(selected_points) - nrow(new_selected_points)
    print(paste("Warning: Removed", removed_count, "invalid selected points."))
  }

  # Return both valid episodes and updated selected points
  return(list(episodes = episodes, selected_points = new_selected_points))
}



# EXCESSES ---------------------------------------------------------------------

#' get_marginal_excess function
#'
#' This function calculates the number of marginal excesses above a given
#' quantile threshold.
#'
#' @param data_rain The rainfall data.
#' @param quantile The quantile value.
#' @param threshold A boolean value to indicate if the quantile variable is a
#'                threshold value and not a uniform quantile. Default is FALSE.
#' @param index The index of the site. Default is 1.
#' @param time The starting time index. Default is 1.
#'
#' @return The number of marginal excesses.
#'
#' @import stats
#'
#' @export
get_marginal_excess <- function(data_rain, quantile, threshold = FALSE,
                                index = 1, time = 1) {
  # shifted data
  Tmax <- nrow(data_rain)
  data <- data_rain[time:Tmax, ] # get the data with time lag

  Tobs <- nrow(data) # number of total observations
  if(!threshold) {
    rain_unif <- rank(data[, index]) / (Tobs + 1) # uniformize the data
  } else {
    rain_unif <- data[, index] # threshold data
  }
  marginal_excesses <- sum(rain_unif > quantile) # number of marginal excesses
  return(marginal_excesses)
}


#' empirical_excesses_rpar function (clean version)
#'
#' Computes empirical excesses for r-Pareto process based on thresholds or quantiles.
#'
#' @param data_rain Matrix or data.frame of rainfall data (rows: time, columns: sites).
#' @param df_lags Dataframe with columns: s1, s2 (site indices), tau (temporal lags).
#' @param thresholds NULL (uniform quantile, default), or a single threshold (numeric), or a numeric vector of thresholds per site.
#' @param quantile Quantile level to use if thresholds is NULL. Default = 0.9.
#' @param t0 Starting time index (integer). Default = 0.
#'
#' @return A data.table with lag values, number of excesses (kij), and number of observed times (Tobs).
#' @export
empirical_excesses_rpar <- function(data_rain, df_lags,
                                    thresholds = NULL, quantile = 0.9,
                                    t0 = 0) {
  # Init output
  excesses <- as.data.table(df_lags[, c("s1", "s2", "tau")])
  excesses[, `:=`(Tobs = 1L, kij = 0L)]

  unique_tau <- unique(df_lags$tau)
  ind_s1 <- df_lags$s1[1]  # site of interest

  for (tau in unique_tau) {
    df_tau <- df_lags[df_lags$tau == tau, ]

    for (i in seq_len(nrow(df_tau))) {
      ind_s2 <- as.numeric(df_tau$s2[i])
      X_s2 <- data_rain[, ind_s2, drop = FALSE]

      if (all(is.na(X_s2))) next
      X_s2 <- as.vector(na.omit(X_s2[, 1]))

      t_index <- t0 + 1 + tau
      if (t_index > length(X_s2)) next

      X_s_t <- X_s2[t_index]
      if (is.list(X_s_t)) X_s_t <- X_s_t[[1]]

      # Determine the threshold to use
      if (!is.null(thresholds)) {
        if (length(thresholds) == 1) {
          threshold_val <- thresholds
        } else {
          threshold_val <- thresholds[ind_s2]
          if (is.na(threshold_val)) next
        }
        exceed <- X_s_t > threshold_val
      } else {
        # Use uniform quantile
        rank_val <- rank(X_s2) / (length(X_s2) + 1)
        exceed <- rank_val[t_index] > quantile
      }
      excesses[s1 == ind_s1 & s2 == ind_s2 & tau == tau,
              kij := as.integer(exceed)]
    }
  }

  return(excesses)
}

# #' empirical_excesses_rpar function
# #'
# #' This function calculates the empirical excesses based on indicators above a
# #' quantile threshold for the r-Pareto process.
# #'
# #' @param data_rain The rainfall data.
# #' @param quantile The quantile value.
# #' @param df_lags The dataframe with spatial and temporal lag values.
# #' @param threshold A boolean value to indicate if the quantile variable is a
# #'                threshold value and not a uniform quantile. Default is FALSE.
# #' @param t0 The starting time. Default is 1.
# #'
# #' @return The empirical excesses dataframe with the number of excesses kij and
# #' the number of possible excesses Tobs, with the lag values.
# #'
# #' @export
# empirical_excesses_rpar <- function(data_rain, quantile, df_lags,
#                                     threshold = FALSE, t0 = 0) {
#   excesses <- as.data.table(df_lags[, c("s1", "s2", "tau")]) # copy the dataframe
#   excesses[, Tobs := 0L]
#   excesses$Tobs <- 1
#   excesses[, kij := 0L]
#   # excesses$kij <- NA
#   unique_tau <- unique(df_lags$tau) # unique temporal lags
#   ind_s1 <- df_lags$s1[1] # s0

#   for (tau in unique_tau) { # loop over temporal lags
#     df_h_t <- df_lags[df_lags$tau == tau, ] # get the dataframe for each tau lag

#     for (i in seq_len(nrow(df_h_t))) { # loop over each pair of sites
#       # get the indices of the sites
#       ind_s2 <- as.numeric(as.character(df_h_t$s2[i]))
#       # get the data for the second site
#       X_s2 <- data_rain[, c(ind_s2), drop = FALSE]
#       if (all(is.na(X_s2))) {
#         next # skip if all values are NA
#       }

#       X_s2 <- as.vector(na.omit(X_s2[, 1]))

#       # if (!threshold) {
#       #   X_s2 <- rank(X_s2) / (length(X_s2) + 1) # uniformize the data
#       # }
#       # shifted data
#       X_s_t <- X_s2[t0 + 1 + tau] # X_{s,t0 + tau}
#       if (is.list(X_s_t)) {
#         X_s_t <- X_s_t[[1]] # extract the value if it's a list
#       }
#       nmargin <- sum(X_s_t > quantile) # 0 or 1
#       excesses$kij[excesses$s1 == ind_s1
#                     & excesses$s2 == ind_s2
#                     & excesses$tau == tau] <- nmargin
#     }
#   }
#   return(excesses)
# }

empirical_excesses_rpar <- function(data_rain, thresholds = NULL, df_lags, t0 = 0) {
  excesses <- as.data.table(df_lags[, c("s1", "s2", "tau")])
  excesses[, Tobs := 1L]
  excesses[, kij := 0L]
  unique_tau <- unique(df_lags$tau)
  ind_s1 <- df_lags$s1[1]

  for (tau in unique_tau) {
    df_h_t <- df_lags[df_lags$tau == tau, ]

    for (i in seq_len(nrow(df_h_t))) {
      ind_s2 <- as.numeric(as.character(df_h_t$s2[i]))
      X_s2 <- data_rain[, c(ind_s2), drop = FALSE]
      if (all(is.na(X_s2))) next

      X_s2 <- as.vector(na.omit(X_s2[, 1]))
      if ((t0 + 1 + tau) > length(X_s2)) next
      X_s_t <- X_s2[t0 + 1 + tau]
      if (is.list(X_s_t)) X_s_t <- X_s_t[[1]]

      if (!is.null(thresholds)) {
        q_s2 <- thresholds[ind_s2]
        if (is.na(q_s2)) next
        nmargin <- as.integer(X_s_t > q_s2)
      } else {
        # fallback: uniformization or binary indicator
        nmargin <- as.integer(X_s_t > 0)  # ou autre logique si tu veux
      }

      excesses$kij[excesses$s1 == ind_s1 &
                   excesses$s2 == ind_s2 &
                   excesses$tau == tau] <- nmargin
    }
  }
  return(excesses)
}



#' empirical_excesses function
#'
#' This function calculates the empirical excesses based on indicators above a
#' quantile threshold.
#'
#' @param data_rain The rainfall data.
#' @param quantile The quantile value.
#' @param df_lags The dataframe with spatial and temporal lag values.
#' @param threshold A boolean value to indicate if the quantile variable is a
#'                 threshold value and not a uniform quantile. Default is FALSE.
#' @param type The type of the process, "rpareto" or "brownresnick".
#'             Default is "rpareto".
#' @param t0 The conditioning time, for r-pareto process. Default is 1.
#'
#' @return The empirical excesses dataframe with the number of excesses kij and
#' the number of possible excesses Tobs, with the lag values.
#'
#' @import tidyr
#'
#' @export
empirical_excesses <- function(data_rain, quantile, df_lags, threshold = FALSE,
                               type = "rpareto", t0 = 0) {
  if (type == "rpareto") {
    excesses <- empirical_excesses_rpar(data_rain, quantile, df_lags, threshold,
                t0)
  } else if (type == "brownresnick")  {
    excesses <- df_lags # copy the dataframe
    unique_tau <- unique(df_lags$tau) # unique temporal lags

    for (t in unique_tau) { # loop over temporal lags
      df_h_t <- df_lags[df_lags$tau == t, ] # get the dataframe for each tau lag

      for (i in seq_len(nrow(df_h_t))) { # loop over each pair of sites
        # get the indices of the sites
        ind_s2 <- as.numeric(as.character(df_h_t$s2[i]))
        ind_s1 <- df_h_t$s1[i]

        # get the data for the pair of sites
        rain_cp <- data_rain[, c(ind_s1, ind_s2), drop = FALSE]
        rain_cp <- as.data.frame(na.omit(rain_cp))
        colnames(rain_cp) <- c("s1", "s2")

        Tmax <- nrow(rain_cp) # number of total observations
        rain_nolag <- rain_cp$s1[1:(Tmax - abs(t))] # get the data without lag
        rain_lag <- rain_cp$s2[(1 + abs(t)):Tmax] # get the data with lag
        Tobs <- length(rain_nolag) # number of observations for the lagged pair
                                   # i.e. T - tau

        # transform the data in uniform data
        if (!threshold) {
          rain_unif <- cbind(rank(rain_nolag) / (Tobs + 1),
                            rank(rain_lag) / (Tobs + 1))
        } else {
          rain_unif <- cbind(rain_nolag, rain_lag)
        }

        # number of joint excesses
        joint_excesses <- sum(rain_unif[, 2] > quantile &
                              rain_unif[, 1] > quantile)

        # store the number of excesses and T - tau
        excesses$Tobs[excesses$s1 == ind_s1
                        & excesses$s2 == ind_s2
                        & excesses$tau == t] <- Tobs

        excesses$kij[excesses$s1 == ind_s1
                      & excesses$s2 == ind_s2
                      & excesses$tau == t] <- joint_excesses
      }
    }
  } else {
    print("The variable 'type' is not valid. It has to be 'rpareto' 
          or 'brownresnick'.")
  }
  return(excesses)
}

# THEORETICAL CHI --------------------------------------------------------------


#' Compute the theoretical chi dataframe.
#'
#' This function calculates the theoretical chi dataframe, based on the given
#' variogram parameters.
#'
#' @param params A vector of variogram parameters.
#' @param df_lags A dataframe with spatial and temporal lag values.
#' @param latlon A boolean value to indicate if the locations are in latitude
#'             and longitude. Default is FALSE.
#' @param directional A boolean value to indicate if the variogram is
#'                  directional. Default is FALSE.
#' @param convert_in_hours A boolean value to indicate if the time lag should be
#'                       converted to hours at the end. Default is TRUE.
#'
#' @return The theoretical chi
#'
#' @export
theoretical_chi <- function(params, df_lags, latlon = FALSE,
                            directional = FALSE) {
  beta1 <- params[1]
  beta2 <- params[2]
  alpha1 <- params[3]
  alpha2 <- params[4]
  adv <- params[5:6]
  chi_df <- df_lags[, c("s1", "s2", "tau", "s1x", "s1y", "s2x", "s2y", "hnorm")]
  # chi_df$tau <- chi_df$tau
  if (latlon) {
    # Compute spatial distance in km (already in km in df_lags)
    haversine_df <- haversine_distance_with_advection(chi_df$s1y, chi_df$s1x,
                                      chi_df$s2y, chi_df$s2x, adv,
                                      (chi_df$tau)) # seconds
    chi_df$hnormV <- haversine_df$distance # km
    # chi_df$theta <- haversine_df$theta
    # chi_df$hx <- chi_df$hnormV * cos(chi_df$theta)
    # chi_df$hy <- chi_df$hnormV * sin(chi_df$theta)

    vector_dist <- vector_distance_with_advection(chi_df$s1y, chi_df$s1x,
                                      chi_df$s2y, chi_df$s2x, adv,
                                      (chi_df$tau)) # seconds
    chi_df$hx <- vector_dist$hx
    chi_df$hy <- vector_dist$hy
  } else {
   # Cartesian coordinates case
   chi_df$s1xv <- chi_df$s1x
   chi_df$s1yv <- chi_df$s1y
   # - adv because tau = t - t0 and 
   # s_shifted = s - adv * t, s0_shifted = s0 - adv * t0, h = s - s0
   chi_df$s2xv <- chi_df$s2x - adv[1] * chi_df$tau # m/s * s = m
   chi_df$s2yv <- chi_df$s2y - adv[2] * chi_df$tau
   
   chi_df$hx <- chi_df$s2xv - chi_df$s1xv # m
   chi_df$hy <- chi_df$s2yv - chi_df$s1yv # m
   # Recompute distance and angle regardless of advection
   chi_df$hnormV <- sqrt((chi_df$s2xv - chi_df$s1xv)^2 +
                        (chi_df$s2yv - chi_df$s1yv)^2) # m
   # Convert to km
  #  chi_df$hx <- chi_df$hx / 1000 # km
  #  chi_df$hy <- chi_df$hy / 1000 # km
  #  chi_df$hnormV <- chi_df$hnormV / 1000 # km
  #  chi_df$theta <- atan2(chi_df$s2yv - chi_df$s1yv, chi_df$s2xv - chi_df$s1xv)

  }

  # if (convert_in_hours && any(chi_df$tau > 3600)) {
  #   # Convert tau to hours
  #   chi_df$tau <- chi_df$tau / 3600 # convert to hours
  # }

  # if (convert_in_km && any(chi_df$hnorm > 1000)) {
  #   # Convert hnormV to km
  #   chi_df$hnormV <- chi_df$hnormV / 1000 # convert to km
  #   chi_df$hx <- chi_df$hx / 1000 # km
  #   chi_df$hy <- chi_df$hy / 1000 # km
  # }


  if (directional) {
    # Apply directional distance adjustment
    chi_df$vario <- (2 * beta1) * abs(chi_df$hx)^alpha1 +
                    (2 * beta1) * abs(chi_df$hy)^alpha1 +
                    (2 * beta2) * abs(chi_df$tau)^alpha2 # same units as V
  } else {
    # Only distance lag
    chi_df$hlag <-  chi_df$hnormV
    chi_df$vario <- (2 * beta1) * abs(chi_df$hlag)^alpha1 +
                    (2 * beta2) * abs(chi_df$tau)^alpha2 # same units as V
  }

  # Compute chi from variogram
  chi_df$chi <- 2 * (1 - pnorm(sqrt(0.5 * chi_df$vario)))
  return(chi_df)
}


#' get_chi_vect function
#'
#' This function calculates the theoretical spatio-temporal extremogram vector
#' based on the given theoretical spatio-temporal extremogram matrix.
#'
#' @param chi_mat The theoretical spatio-temporal extremogram matrix.
#' @param h_vect The spatial lags vector.
#' @param tau The temporal lag vector.
#' @param df_dist The distances long dataframe.
#'
#' @return The theoretical spatio-temporal extremogram vector.
#'
#' @export
get_chi_vect <- function(chi_mat, h_vect, tau, df_dist) {
  df_dist$h <- ifelse(df_dist$value %in% h_vect, df_dist$value, NA)
  # create excesses matrix
  npairs <- nrow(df_dist) # number of pairs
  nconfig <- npairs * length(tau) # number of configurations
  h_vect_all <- rep(df_dist$h, times = length(tau))
  tau_all <- rep(tau, each = npairs) # repeat tau for each pair
  chi_vect <- numeric(length = nconfig) # init chi vector
  for (p in 1:nconfig) {
    h <- h_vect_all[p] # get h for each configuration
    if (is.na(h)) { # if NA then chi = NA
      chi_vect[p] <- NA
    } else { # else get chi value
      t <- tau_all[p] # get tau for each configuration
      h_ind <- which(h_vect == h) # get index of h
      chi <- chi_mat[t, h_ind] # get chi value from matrix
      # chi inside ]0, 1] to avoid future error
      chi_vect[p] <- ifelse(chi <= 0, 0.000001, chi) # get chi value in vector
    }
  }
  return(chi_vect)
}

# NEGATIVE LOG-LIKELIHOOD ------------------------------------------------------

#' neg_ll function
#'
#' Calculate the negative log-likelihood for a given set of variogram
#' parameters by considering excess indicators following a binomial distribution
#' with a probability parameter equals to the theoretical spatio-temporal chi
#' value.
#'
#' @param params Vector of variogram parameters (beta1, beta2, alpha1, alpha2).
#' @param df_lags The dataframe with spatial and temporal lag values.
#' @param excesses The excesses dataframe with the number of excesses kij and
#'                 the number of possible excesses Tobs.
#' @param hmax The maximum spatial lag value. Default is NA.
#' @param rpar A boolean value to indicate if the data is considered as a
#'            r-Pareto process. Default is TRUE.
#' @param threshold A boolean value to indicate if the quantile variable is a
#'                threshold value and not a uniform quantile. Default is FALSE.
#' @param latlon A boolean value to indicate if the locations are in latitude
#'             and longitude. Default is FALSE.
#' @param quantile The quantile value. Default is NA.
#' @param data The data dataframe (brown-resnick case). Default is NA.
#' @param directional A boolean value to indicate if the variogram is
#'                   directional. Default is TRUE.
#' @param convert_in_hours A boolean value to indicate if the time lag should be
#'                       converted to hours at the end. Default is TRUE.
#' @param convert_in_km A boolean value to indicate if the spatial lag should be
#'                    converted to kilometers at the end. Default is TRUE.
#'
#' @return The negative log-likelihood value.
#'
#' @export
neg_ll <- function(params, df_lags, excesses,
                      hmax = NA,  rpar = TRUE, threshold = FALSE,
                      latlon = FALSE, quantile = NA, data = NA,
                      directional = TRUE) {
  if (rpar) { # if we have a conditioning location
    p <- 1 # sure excess for r-Pareto process in (s0,t0)
  } else {
    if (is.na(data)) {
      stop("The data must be provided for the max-stable composite likelihood.")
    } else if (is.na(quantile)) {
      stop("The quantile or threshold must be provided for the
            max-stable composite likelihood.")
    }
    Tmax <- nrow(data) # number of total observations
    # number of marginal excesses
    nmarg <- get_marginal_excess(data, quantile, threshold)
    p <- nmarg / Tmax # probability of marginal excesses
  }
  
  # compute theoretical chi values
  # df_lags$hnorm <- df_lags$hnorm / 1000 # convert to km
  if (!is.na(hmax)) {
    excesses <- excesses[df_lags$hnorm <= hmax, ]
    df_lags <- df_lags[df_lags$hnorm <= hmax, ]
  }
  chi <- theoretical_chi(params, df_lags, latlon, directional)
  ll_df <- df_lags # copy the dataframe
  ll_df$kij <- excesses$kij # number of excesses
  ll_df$Tobs <- excesses$Tobs
  # ll_df$hnormV <- chi$hnormV
  ll_df$chi <- chi$chi
  # ll_df$chi <- pmax(pmin(ll_df$chi, 1 - 1e-8), 1e-8)
  ll_df$chi <- ifelse(ll_df$chi <= 0, 1e-10, ll_df$chi)
  ll_df$pchi <- 1 - p * ll_df$chi
  # ll_df$pchi <- pmax(pmin(ll_df$pchi,  1 - 1e-8), 1e-8)
  ll_df$pchi <- ifelse(ll_df$pchi <= 0, 1e-10, ll_df$pchi)

  # number of non-excesses
  ll_df$non_excesses <- ll_df$Tobs - ll_df$kij
  ll_df$ll <- ll_df$kij * log(ll_df$chi) +
              ll_df$non_excesses * log(ll_df$pchi)
  # if (!is.na(hmax)) {
  #   ll_df <- ll_df[ll_df$hnorm <= hmax, ]
  # }

  nll <- -sum(ll_df$ll, na.rm = TRUE)
  return(nll)
}



#' neg_ll_composite function
#'
#' Calculate the negative log-likelihood for a list of r-Pareto processes
#' with wind covariates if wanted.
#'
#' @param params Vector of variogram parameters (beta1, beta2, alpha1, alpha2
#'              eta1, eta2).
#' @param list_lags A list of dataframes with spatial and temporal lag values.
#' @param list_excesses A list of excesses dataframes.
#' @param wind_df A dataframe with wind values. Default is NA.
#' @param hmax The maximum spatial lag value. Default is NA.
#' @param latlon A boolean value to indicate if the locations are in latitude
#'             and longitude. Default is TRUE.
#' @param directional A boolean value to indicate if the variogram is
#'                    directional. Default is TRUE.
#' @param threshold A boolean value to indicate if the quantile variable is a
#'               threshold value and not a uniform quantile. Default is FALSE.
#' @param rpar A boolean value to indicate if the data is considered as a
#'          r-Pareto process. Default is TRUE.
#' @param fixed_eta1 A boolean value to indicate if eta1 is fixed. Default is
#'                   FALSE.
#' @param fixed_eta2 A boolean value to indicate if eta2 is fixed. Default is
#'                  FALSE.
#' @param convert_in_hours A boolean value to indicate if the time lag should be
#'                      converted to hours at the end. Default is FALSE.
#' @param convert_in_km A boolean value to indicate if the spatial lag should be
#'                     converted to kilometers at the end. Default is FALSE.
#'
#' @return The negative log-likelihood value.
#'
#' @export
neg_ll_composite <- function(params, list_episodes, list_excesses,
                             list_lags, wind_df = NA,
                             hmax = NA, latlon = TRUE,
                             directional = TRUE, threshold = FALSE,
                             rpar = TRUE, fixed_eta1 = FALSE,
                             fixed_eta2 = FALSE) {
  # Add default values for advection parameters
  print(params)
  if (length(params) == 4) {
    params <- c(params, 0, 0) # No advection
  }
  if(rpar) {
    data <- NA
    quantile <- NA
  }

  if (all(is.na(wind_df))) {
    adv <- params[5:6]
    adv_df <- adv
  } else {
    eta1 <- params[5]
    eta2 <- params[6]
    if (fixed_eta1) {
      eta1 <- 1
    }
    if (fixed_eta2) {
      eta2 <- 1
    }
    adv_x <- (abs(wind_df$vx)^eta2) * sign(wind_df$vx) * eta1
    adv_y <- (abs(wind_df$vy)^eta2) * sign(wind_df$vy) * eta1
    # adv_x <- eta1 * wind_df$vx*(abs(wind_df$vx)^(eta2 - 1))
    # adv_y <- eta1 *  wind_df$vy*(abs(wind_df$vy)^(eta2 - 1))
    adv_df <- cbind(adv_x, adv_y)
    if (length(adv_df) == 2) {
      adv <- as.vector(adv_df)
    }
  }


  m <- length(list_episodes) # Number of r-Pareto processes ie episodes
  nll_composite <- 0
  for (i in 1:m) {
    # extract episode and excesses from i-th r-pareto process from data
    excesses <- list_excesses[[i]]
    lags <- list_lags[[i]]

    if (!all(is.na(wind_df)) && length(adv_df) != 2) {
      adv <- as.vector(adv_df[i,])
    }

    params_adv <- c(params[1:4], adv) # Add advection parameters
    if (!rpar) {
      data <- list_episodes[[i]]
      quantile <- quantile
    }

    nll_i <- neg_ll(params = params_adv,
                     df_lags = lags,
                     hmax = hmax, excesses = excesses,
                     latlon = latlon, directional = directional,
                     threshold = threshold, rpar = rpar,
                     data = data, quantile = quantile)

    nll_composite <- nll_composite + nll_i
  }

  if (!is.finite(nll_composite)) {
    cat("Warning: Non-finite value detected for parameters:", params, "\n")
  }

  return(nll_composite)
}

#' neg_ll_composite_mem function
#'
#' Calculate the negative log-likelihood for a list of r-Pareto processes
#' with wind covariates if wanted.
#'
#' @param params Vector of variogram parameters (beta1, beta2, alpha1, alpha2
#'              eta1, eta2).
#' @param list_episodes A list of episodes (r-Pareto processes).
#' @param s0_list A list of s0 coordinates (latitude and longitude).
#' @param sites_coords A dataframe with the coordinates of the sites.
#' @param wind_df A dataframe with the wind data (vx, vy).
#' @param hmax The maximum spatial lag value. Default is NA.
#' @param latlon A boolean value to indicate if the locations are in latitude
#'            and longitude. Default is FALSE.
#' @param directional A boolean value to indicate if the variogram is
#'                 directional. Default is FALSE.
#' @param threshold A boolean value to indicate if the quantile variable is a
#'               threshold value and not a uniform quantile. Default is FALSE.
#' @param rpar A boolean value to indicate if the data is considered as a
#'           r-Pareto process. Default is TRUE.
#' @param t0 The starting time. Default is 0.
#' @param tau_vect The vector of temporal lags. Default is 0:10.
#' @param fixed_eta1 A boolean value to indicate if eta1 is fixed. Default is
#'                   FALSE.
#' @param fixed_eta2 A boolean value to indicate if eta2 is fixed. Default is
#'                   FALSE.
#'
#' @return The negative log-likelihood value.
#'
#' @export
neg_ll_composite_mem <- function(params, list_episodes, s0_list, sites_coords,
                             wind_df = NA,
                             hmax = NA, latlon = TRUE,
                             directional = TRUE, threshold = FALSE,
                             rpar = TRUE, t0 = 0, tau_vect = 0:10,
                             fixed_eta1 = FALSE,
                             fixed_eta2 = FALSE) {
  # Add default values for advection parameters
  # print(params)
  if (length(params) == 4) {
    params <- c(params, 0, 0) # No advection
  }
  if(rpar) {
    data <- NA
    quantile <- NA
  }
  if (all(is.na(wind_df))) {
    adv <- params[5:6]
    adv_df <- adv
  } else {
    eta1 <- params[5]
    eta2 <- params[6]
    if (fixed_eta1) {
      eta1 <- 1
    }
    if (fixed_eta2) {
      eta2 <- 1
    }
    adv_x <- (abs(wind_df$vx)^eta1) * sign(wind_df$vx) * eta2
    adv_y <- (abs(wind_df$vy)^eta1) * sign(wind_df$vy) * eta2
    adv_df <- cbind(adv_x, adv_y)
    if (length(adv_df) == 2) {
      adv <- as.vector(adv_df)
    }
  }

  m <- length(list_episodes) # Number of r-Pareto processes ie episodes
  nll_composite <- 0
  for (i in 1:m) {
    # extract episode and excesses from i-th r-pareto process from data
    # excesses <- list_excesses[[i]]
    # lags <- list_lags[[i]]
    s0_x <- s0_list[[i]]$x
    s0_y <- s0_list[[i]]$y
    s0_coords <- sites_coords[sites_coords$Longitude == s0_x &
                              sites_coords$Latitude == s0_y, ]
    episode <- list_episodes[[i]]
    lags <- get_conditional_lag_vectors(sites_coords, s0_coords,
                                        t0, tau_vect, latlon = FALSE)
    excesses <- empirical_excesses_rpar(episode, quantile = 1, lags,
                                     t0 = t0, threshold = TRUE)

    if (!all(is.na(wind_df)) && length(adv_df) != 2) {
      adv <- as.vector(adv_df[i,])
    }

    params_adv <- c(params[1:4], adv) # Add advection parameters
    if (!rpar) {
      data <- list_episodes[[i]]
      quantile <- quantile
    }

    nll_i <- neg_ll(params = params_adv,
                     df_lags = lags,
                     hmax = hmax, excesses = excesses,
                     latlon = latlon, directional = directional,
                     threshold = threshold, rpar = rpar,
                     data = data, quantile = quantile)

    nll_composite <- nll_composite + nll_i
  }

  if (!is.finite(nll_composite)) {
    cat("Warning: Non-finite value detected for parameters:", params, "\n")
  }

  return(nll_composite)
}



# WIND COVARIATES --------------------------------------------------------------

#' convert_to_cardinal function
#'
#' Convert wind direction in degrees to cardinal direction between N, NE, E, SE,
#' S, SW, W, NW.
#'
#' @param degrees The wind direction in degrees.
#'
#' @return The cardinal direction.
#'
#' @export
convert_to_cardinal <- function(degrees, nb_cardinal = 8) {
  if (is.na(degrees)) {
    return(NA)  # Return NA if the input is NA
  }
  if (nb_cardinal < 4 || nb_cardinal > 8) {
    stop("nb_cardinal must be between 4 and 8.")
  }

  if (nb_cardinal == 4) {
    # For 4 cardinal directions: N, E, S, W
    directions <- c("N", "E", "S", "W", "N")
    breaks <- c(0, 90, 180, 270, 360)
  } else if (nb_cardinal == 8) {
    # For 8 cardinal directions: N, NE, E, SE, S, SW, W, NW
    directions <- c("N", "NE", "E", "SE", "S", "SW", "W", "NW", "N")
    breaks <- c(0, 45, 90, 135, 180, 225, 270, 315, 360)
  }

  return(directions[findInterval(degrees, breaks, rightmost.closed = TRUE)])
}

#' cardinal_to_degree function
#'
#' Convert cardinal direction between N, NE, E, SE, S, SW, W, NW to degrees.
#'
#' @param cardinal The cardinal direction.
#'
#' @return The wind direction in degrees.
#'
#' @export
cardinal_to_degree <- function(cardinal, nb_cardinal = 8) {
  if (is.na(cardinal)) {
    return(NA)  # Return NA if the input is NA
  }

  if (nb_cardinal < 4 || nb_cardinal > 8) {
    stop("nb_cardinal must be between 4 and 8.")
  }

  if (nb_cardinal == 4) {
    # For 4 cardinal directions: N, E, S, W
    directions <- c("N", "E", "S", "W")
    degrees <- c(0, 90, 180, 270)
  } else if (nb_cardinal == 8) {
    # For 8 cardinal directions: N, NE, E, SE, S, SW, W, NW
    directions <- c("N", "NE", "E", "SE", "S", "SW", "W", "NW")
    degrees <- c(0, 45, 90, 135, 180, 225, 270, 315)
  }

  return(degrees[match(cardinal, directions)])
}

#' get_mode_dir function
#'
#' Get the most frequent wind direction in a vector of wind directions.
#'
#' @param x The wind data
#'
#' @return The most frequent wind direction.
#'
#' @export
get_mode_dir <- function(x) {
  x <- na.omit(x)  # Remove NA values
  uniq_values <- unique(x)  # Get unique values
  freq_table <- tabulate(match(x, uniq_values))  # Count occurrences
  mode_value <- uniq_values[which.max(freq_table)]  # Most frequent value
  return(mode_value)
}

#' compute_wind_episode function
#'
#' Compute the wind speed and direction for a given extreme episode.
#'
#' @param episode The extreme episode dataframe.
#' @param s0 The conditioning location name.
#' @param u The quantile threshold.
#' @param wind_df The wind dataframe.
#' @param speed_time The time index for wind speed.
#'
#' @return The wind episode dataframe.
#'
#' @export
compute_wind_episode <- function(episode, s0, u, wind_df, speed_time) {
  timestamps <- as.POSIXct(rownames(episode),
                            format = "%Y-%m-%d %H:%M:%S", tz = "UTC")

  # Subset the wind data for the episode
  wind_subset <- wind_df %>%
    dplyr::filter(datetime %in% timestamps)

  # get timestamp for which there is an excess above u in s0 in episode
  t_excesses <- which(episode[, s0] > u)
  
  wind_subset_excess <- wind_df %>%
    dplyr::filter(datetime %in% t_excesses)

  # Compute the mean wind speed and direction for the episode if there is data
  if (nrow(wind_subset) > 0) {
    FF <- wind_subset$FF[speed_time + 1]
    DD_t0 <- wind_subset$DD[speed_time + 1]
    DD_mean <- mean(wind_subset$DD, na.rm = TRUE)
    FF_mean <- mean(wind_subset$FF, na.rm = TRUE)
    cardDir_t0 <- wind_subset$cardDir[speed_time + 1]
  } else {
    FF <- NA
    DD_t0 <- NA
    cardDir_t0 <- NA
    DD_mean <- NA
    FF_mean <- NA
  }

  return(data.frame(FF = FF,
                    DD_t0 = DD_t0, cardDir_t0 = cardDir_t0,
                    FF_mean = FF_mean, DD_mean = DD_mean))
}



# Fonction pour mode robuste
get_mode_or_t0 <- function(x, fallback) {
  x <- na.omit(x)
  if (length(x) == 0) return(fallback)
  tab <- table(x)
  modes <- names(tab[tab == max(tab)])
  if (length(modes) == 1) {
    return(modes)
  } else {
    return(fallback)
  }
}


compute_wind_episode <- function(episode, s0, u, wind_df, speed_time,
                                 nb_cardinal = 8) {
  timestamps <- as.POSIXct(rownames(episode),
                            format = "%Y-%m-%d %H:%M:%S", tz = "UTC")

  # Subset du vent sur toute la durée de l’épisode
  wind_subset <- wind_df %>%
    dplyr::filter(datetime %in% timestamps)

  # Timestamps des excès
  t_excesses <- which(episode[, s0] > u)
  timestamps_excess <- timestamps[t_excesses]

  # Vent pendant les excès
  wind_subset_excess <- wind_df %>%
    dplyr::filter(datetime %in% timestamps_excess)

  # valeur à t0
  FF <- DD_t0 <- cardDir_t0 <- FF_mean <- DD_mean  <- DD_excess <- NA
  FF_mean_excess <- cardDir_mode_excess <- NA

  if (nrow(wind_subset) > 0) {
    FF <- wind_subset$FF[speed_time + 1]
    DD_t0 <- wind_subset$DD[speed_time + 1]
    cardDir_t0 <- wind_subset$cardDir[speed_time + 1]

    FF_mean <- mean(wind_subset$FF, na.rm = TRUE)
    DD_mean <- mean(wind_subset$DD, na.rm = TRUE)
  }



  if (nrow(wind_subset_excess) > 0) {
    FF_mean_excess <- mean(wind_subset_excess$FF, na.rm = TRUE)
    cardDir_mode_excess <- get_mode_or_t0(wind_subset_excess$cardDir, fallback = cardDir_t0)
    DD_excess <- cardinal_to_degree(cardDir_mode_excess, nb_cardinal = nb_cardinal)
  }

  return(data.frame(
    FF = FF,
    DD_t0 = DD_t0,
    cardDir_t0 = cardDir_t0,
    FF_mean = FF_mean,
    DD_mean = DD_mean,
    FF_mean_excess = FF_mean_excess,
    DD_excess = DD_excess,
    cardDir_mode_excess = cardDir_mode_excess
  ))
}



compute_wind_episode <- function(episode, s0, u, wind_df, speed_time,
                                 nb_cardinal = 8) {
  timestamps <- as.POSIXct(rownames(episode),
                            format = "%Y-%m-%d %H:%M:%S", tz = "UTC")

  # Subset du vent sur toute la durée de l’épisode
  wind_subset <- wind_df %>%
    dplyr::filter(datetime %in% timestamps)

  # Timestamps des excès
  t_excesses <- which(episode[, s0] > u)
  timestamps_excess <- timestamps[t_excesses]

  # Get pondering mean for wind direction DD with the force FF
  FF <- DD <- NA

  return(data.frame(
    FF = FF,
    DD = DD,
    cardDir_t0 = cardDir_t0,
  ))
}




# SIMU -------------------------------------------------------------------------

#' simulate_excess_ind function
#'
#' Simulate individual excess following a binomial distribution with a given
#' probability parameter equals to the chi value.
#'
#' @param Tmax The maximum time observation.
#' @param chi_h_t The spatio-temporal extremogram value at a given spatial lag h
#' and a given temporal lag t.
#' @param p_marg The probability of marginal excesses.
#'
#' @return The simulated individual excess.
#'
#' @import stats
#'
#' @export
simulate_excess_ind <- function(Tmax, chi_h_t, p_marg) {
  # simulate excesses following a binomial distribution
  vect_E <- rbinom(Tmax, 1, p_marg * chi_h_t)
  return(vect_E)
}

# VALIDATION -------------------------------------------------------------------


#' get_criterion function
#'
#' Calculates the mean, RMSE and MAE criterion values for the given dataframe
#' result and the true parameter values.
#'
#' @param df_result The data frame containing the result values
#' @param true_param The true variogram parameter (beta1, beta2, alpha1, alpha2)
#' @return The calculated criterion value in a dataframe.
#'
#' @export
get_criterion <- function(df_result, true_param) {
  # remove NA values
  df_result <- na.omit(df_result)
  # get the mean, RMSE and MAE values for each parameter
  mean_beta1 <- mean(df_result$beta1)
  mean_alpha1 <- mean(df_result$alpha1)
  mean_beta2 <- mean(df_result$beta2)
  mean_alpha2 <- mean(df_result$alpha2)

  rmse_beta1 <- sqrt(mean((true_param[1] - df_result$beta1)^2))
  rmse_alpha1 <- sqrt(mean((true_param[3] - df_result$alpha1)^2))
  rmse_beta2 <- sqrt(mean((true_param[2] - df_result$beta2)^2))
  rmse_alpha2 <- sqrt(mean((true_param[4] - df_result$alpha2)^2))

  mae_beta1 <- mean(abs(true_param[1] - df_result$beta1))
  mae_alpha1 <- mean(abs(true_param[3] - df_result$alpha1))
  mae_beta2 <- mean(abs(true_param[2] - df_result$beta2))
  mae_alpha2 <- mean(abs(true_param[4] - df_result$alpha2))

  if (length(true_param) == 6) {
    mean_adv1 <- mean(df_result$adv1)
    mean_adv2 <- mean(df_result$adv2)

    rmse_adv1 <- sqrt(mean((true_param[5] - df_result$adv1)^2))
    rmse_adv2 <- sqrt(mean((true_param[6] - df_result$adv2)^2))

    mae_adv1 <- mean(abs(true_param[5] - df_result$adv1))
    mae_adv2 <- mean(abs(true_param[6] - df_result$adv2))

    df_valid <- data.frame(mean = c(mean_beta1, mean_beta2, mean_alpha1,
                                mean_alpha2, mean_adv1, mean_adv2),
                        rmse = c(rmse_beta1, rmse_beta2, rmse_alpha1,
                                rmse_alpha2, rmse_adv1, rmse_adv2),
                        mae = c(mae_beta1, mae_beta2, mae_alpha1, mae_alpha2,
                                mae_adv1, mae_adv2),
                        row.names = c("beta1", "beta2", "alpha1", "alpha2",
                                      "adv1", "adv2"))
  } else {
    df_valid <- data.frame(mean = c(mean_beta1, mean_beta2, mean_alpha1,
                                mean_alpha2),
                        rmse = c(rmse_beta1, rmse_beta2, rmse_alpha1,
                                rmse_alpha2),
                        mae = c(mae_beta1, mae_beta2, mae_alpha1, mae_alpha2),
                        row.names = c("beta1", "beta2", "alpha1", "alpha2"))
  }
  return(df_valid)
}


#' save_results_optim function
#'
#' Save the results of the optimization process in a CSV file.
#'
#' @param result The result of the optimization process.
#' @param true_param The true variogram parameter (beta1, beta2, alpha1, alpha2)
#' @param filename The name of the file to save the results.
#' @param data_folder The folder where the results will be saved. Default is NA.
#'
#' @import utils
#'
#' @export
save_results_optim <- function(result, true_param, filename, data_folder = NA) {
  if (is.na(data_folder)) {
    data_folder <- ""
  } else {
    data_folder <- paste0(data_folder, "/simulations/simulation_rpar/")
  }

  if (result$convergence == 0) {
    rmse <- sqrt((result$par - true_param)^2)
    df_rmse <- data.frame(estim = result$par, rmse = rmse)
    rownames(df_rmse) <- c("beta1", "beta2", "alpha1", "alpha2", "Vx", "Vy")
    # save the results
    utils::write.csv(t(df_rmse), file = paste0(data_folder,
                                                    filename, ".csv"))
  } else {
    print("No convergence")
  }
}

#' get_results_optim function
#'
#' Get the results of the optimization process from a CSV file.
#'
#' @param filename The name of the file to get the results.
#' @param data_folder The folder where the results are saved. Default is NA.
#'
#' @import utils
#' 
#' @export
get_results_optim <- function(filename, data_folder = NA) {
  if (is.na(data_folder)) {
    data_folder <- ""
  } else {
    data_folder <- paste0(data_folder, "/simulations/simulation_rpar/")
  }
  # read the results
  df_rmse <- utils::read.csv(paste0(data_folder, filename,
                      ".csv"))
  return(df_rmse)
}


#' process_simulation function
#'
#' This function processes the simulation data and optimizes the parameters
#' using the negative log-likelihood function.
#'
#' @param i The index of the simulation.
#' @param m The number of simulations to process.
#' @param list_episodes A list of simulated data.
#' @param u The quantile threshold.
#' @param s0_list A list of conditioning locations.
#' @param init_params The initial parameters for optimization.
#' @param sites_coords The coordinates of the sites.
#' @param hmax The maximum spatial lag value. Default is NA.
#' @param wind_df A dataframe with wind values. Default is NA.
#' @param directional A boolean value to indicate if the variogram is
#'                   directional. Default is FALSE.
#' @param t0 The starting time. Default is 0.
#' @param tau_vect The vector of temporal lags. Default is 0:10.
#' @param fixed_eta1 A boolean value to indicate if eta1 is fixed. Default is
#'                   FALSE.
#' @param fixed_eta2 A boolean value to indicate if eta2 is fixed. Default is
#'                   FALSE.
#'
#' @return The optimized parameters.
#'
#' @export
process_simulation_mem <- function(list_episodes, u, s0_list,
                               init_params, sites_coords, hmax = NA,
                               wind_df = NA, directional = FALSE,
                               t0 = 0, tau_vect = 0:10,
                               fixed_eta1 = FALSE,
                               fixed_eta2 = FALSE) {
  # Bounds
  lower_bounds <- c(1e-8, 1e-8, 1e-8, 1e-8, -Inf, -Inf)
  upper_bounds <- c(Inf, Inf, 1.999, 1.999, Inf, Inf)

  # If there is wind data, set bounds for eta1 and eta2
  if (all(!is.na(wind_df))) {
    lower_bounds[5:6] <- c(1e-8, 1e-8)
  }

  # Optimize the negative composite log-likelihood function
  result <- optim(
    par = init_params,
    fn = neg_ll_composite_mem,
    list_episodes = list_episodes,
    s0_list = s0_list,
    hmax = hmax,
    sites_coords = sites_coords,
    wind_df = wind_df,
    threshold = TRUE,
    latlon = FALSE,
    directional = directional,
    t0 = t0,
    tau_vect = tau_vect,
    fixed_eta1 = fixed_eta1,
    fixed_eta2 = fixed_eta2,
    method = "L-BFGS-B",
    lower = lower_bounds,
    upper = upper_bounds,
    control = list(maxit = 10000)
  )

  # Check for convergence
  estimates <- if (result$convergence == 0) result$par else rep(NA, length(init_params))

  return(estimates)
}



#' process_simulation function
#'
#' This function processes the simulation data and optimizes the parameters
#' using the negative log-likelihood function.
#'
#' @param i The index of the simulation.
#' @param m The number of simulations to process.
#' @param list_simu The list of simulations.
#' @param u The quantile threshold.
#' @param list_lags The list of lags.
#' @param list_excesses The list of excesses.
#' @param init_params The initial parameters for optimization.
#' @param hmax The maximum spatial lag value. Default is NA.
#' @param wind_df The wind dataframe. Default is NA.
#' @param directional Directional variogram or not. Default is FALSE.
#'
#' @return The optimized parameters.
#'
#' @export
process_simulation <- function(i, m, list_simu, u, list_lags,
                               list_excesses,
                               init_params, hmax = NA, wind_df = NA,
                               directional = FALSE,
                               fixed_eta1 = FALSE,
                               fixed_eta2 = FALSE) {
  # Bounds
  lower_bounds <- c(1e-8, 1e-8, 1e-8, 1e-8, -Inf, -Inf)
  upper_bounds <- c(Inf, Inf, 1.999, 1.999, Inf, Inf)
  if (all(!is.na(wind_df))) {
    lower_bounds[5:6] <- c(1e-8, 1e-8)
  }
  # Get the m corresponding simulations from list_simu inside a list
  list_episodes <- list_simu[((i - 1) * m + 1):(i * m)]
  lags_episodes <- list_lags[((i - 1) * m + 1):(i * m)]
  excesses_episodes <- list_excesses[((i - 1) * m + 1):(i * m)]

  # Optimize
  result <- optim(
    par = init_params,
    fn = neg_ll_composite,
    list_episodes = list_episodes,
    list_lags = lags_episodes,
    list_excesses = excesses_episodes,
    hmax = hmax,
    wind_df = wind_df,
    threshold = TRUE,
    latlon = FALSE,
    directional = directional,
    fixed_eta1 = fixed_eta1,
    fixed_eta2 = fixed_eta2,
    method = "L-BFGS-B",
    lower = lower_bounds,
    upper = upper_bounds,
    control = list(maxit = 10000),
    hessian = TRUE
  )

  estimates <- result$par
  conv <- (result$convergence == 0)
  if (!conv) {
    estimates <- rep(NA, length(init_params))
    ci_lower <- ci_upper <- rep(NA, length(init_params))
    results_df <- data.frame(
      beta1 = estimates[1],
      beta2 = estimates[2],
      alpha1 = estimates[3],
      alpha2 = estimates[4],
      adv1 = estimates[5],
      adv2 = estimates[6],
      beta1_lower = NA,
      beta1_upper = NA,
      beta2_lower = NA,
      beta2_upper = NA,
      alpha1_lower = NA,
      alpha1_upper = NA,
      alpha2_lower = NA,
      alpha2_upper = NA,
      adv1_lower = NA,
      adv1_upper = NA,
      adv2_lower = NA,
      adv2_upper = NA
    )
  } else {
    # If convergence, extract the Hessian matrix and compute CIs
    hessian <- result$hessian
    # Check if hessian is valid
    if (is.null(hessian) || any(!is.finite(hessian)) || inherits(try(solve(hessian), silent = TRUE), "try-error")) {
      estimates <- rep(NA, length(init_params))
      ci_lower <- ci_upper <- rep(NA, length(init_params))
    } else {
      # Compute standard errors from inverse Hessian
      cov_matrix <- solve(hessian)
      se <- sqrt(diag(cov_matrix))
      
      # 95% confidence intervals
      z_val <- qnorm(0.975)  # ≈1.96
      ci_lower <- estimates - z_val * se
      ci_upper <- estimates + z_val * se
      
      # Optionally, combine into a data frame
      results_df <- data.frame(
        beta1 = estimates[1],
        beta1_lower = ci_lower[1],
        beta1_upper = ci_upper[1],
        beta2 = estimates[2],
        beta2_lower = ci_lower[2],
        beta2_upper = ci_upper[2],
        alpha1 = estimates[3],
        alpha1_lower = ci_lower[3],
        alpha1_upper = ci_upper[3],
        alpha2 = estimates[4],
        alpha2_lower = ci_lower[4],
        alpha2_upper = ci_upper[4],
        adv1 = estimates[5],
        adv1_lower = ci_lower[5],
        adv1_upper = ci_upper[5],
        adv2 = estimates[6],
        adv2_lower = ci_lower[6],
        adv2_upper = ci_upper[6]
      )
    }
  }
  if (!is.na(wind_df)) {
    # change adv1, adv2 and their CIs to eta1, eta2 and their CIs
    results_df$eta1 <- results_df$adv1
    results_df$eta2 <- results_df$adv2
    results_df$eta1_lower <- results_df$adv1_lower
    results_df$eta1_upper <- results_df$adv1_upper
    results_df$eta2_lower <- results_df$adv2_lower
    results_df$eta2_upper <- results_df$adv2_upper
    results_df <- results_df[, -c("adv1", "adv2", "adv1_lower", "adv1_upper",
                                  "adv2_lower", "adv2_upper")]
  }
  return(results_df)
}
