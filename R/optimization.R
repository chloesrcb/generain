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
#' @param remove_zeros A boolean indicating whether to remove zeros when computing quantiles.
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
        if (length(col) < 30) {
          return(NA)
        }
        stats::quantile(col, probs = quantile, na.rm = TRUE)
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
#' @param beta The temporal window extension. Default is 0.
#' @return A data.table with selected pairs of sites and times (s0, t0) and their
#'         corresponding excess values u_s0.
#' 
#' @import data.table
#' @import geosphere
#' @export
get_s0t0_pairs <- function(sites_coords, data, min_spatial_dist,
                           episode_size, set_st_excess,
                           n_max_episodes = 10000, latlon = FALSE,
                           beta = 0) {

  data <- as.matrix(data)
  site_names <- colnames(data)

  if (is.null(site_names)) stop("ERROR: 'data' must have column names.")
  coord_names <- rownames(sites_coords)
  if (is.null(coord_names)) stop("ERROR: 'sites_coords' must have row names.")

  missing_in_coords <- setdiff(site_names, coord_names)
  if (length(missing_in_coords) > 0){
    stop(paste("Missing sites in coords:", paste(missing_in_coords, collapse = ", ")))
  }
  # Reorder coords
  sites_coords <- sites_coords[site_names, , drop = FALSE]

  # Distance matrix
  if (latlon) {
    dist_matrix <- get_dist_mat(sites_coords, latlon = TRUE) # meters
    # convert in km
    dist_matrix <- dist_matrix / 1000
  } else {
    dist_matrix <- get_dist_mat(sites_coords, latlon = FALSE)
    # dist_matrix <- ceiling(round(dist_matrix, 1) * 1000) / 1000
  }
  rownames(dist_matrix) <- site_names
  colnames(dist_matrix) <- site_names

  # Extract lists
  list_s <- set_st_excess$list_s
  list_t <- set_st_excess$list_t
  list_u <- set_st_excess$list_u

  # Initialize
  selected_points <- data.table(s0 = character(), t0 = integer(),
                                t0_date = as.character(NA), u_s0 = numeric())
  nb_episode <- 0

  # Combine and order excesses
  excess_dt <- data.table(
    s = unlist(list_s),
    t = unlist(list_t),
    t_dates = names(list_t),
    u = unlist(list_u)
  )
  if (is.null(excess_dt$t_dates)) excess_dt[, t_dates := NA_character_]
  setorder(excess_dt, t)

  # Loop through events
  for (i in seq_len(nrow(excess_dt))) {
    s0 <- excess_dt$s[i]
    t0 <- excess_dt$t[i]
    t0_date <- excess_dt$t_dates[i]
    u0 <- excess_dt$u[i]

    is_valid <- TRUE
    for (j in seq_len(nrow(selected_points))) {
      s_prev <- selected_points$s0[j]
      t_prev <- selected_points$t0[j]
      if (is.na(s_prev) || is.na(t_prev)) next
      spatial_dist <- dist_matrix[s0, s_prev]
      temporal_dist <- abs(t0 - t_prev)
      if (spatial_dist < min_spatial_dist && temporal_dist < (episode_size + 2 * beta)) {
        is_valid <- FALSE
        break
      }
    }

    if (is_valid) {
      new_row <- data.table(s0 = s0, t0 = t0,
                            t0_date = ifelse(is.na(t0_date), NA_character_, t0_date),
                            u_s0 = u0)
      selected_points <- rbind(selected_points, new_row, fill = TRUE)
      nb_episode <- nb_episode + 1
      if (nb_episode >= n_max_episodes) break
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
#' @param unif A boolean value to indicate if we want to get the uniformized
#'            data episodes. Default is FALSE.
#'
#' @return The list of extreme episodes of size 2 * delta each.
#'
#' @export
get_extreme_episodes <- function(selected_points, data, episode_size,
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

    t_inf <- max(1, t0)
    t_sup <- min(nrow(data), t0 + episode_size - 1)
    # print(i)
    # Check that the episode is the correct size (episode_size)
    episode_size_test <- t_sup - t_inf + 1
    if (episode_size_test == episode_size) {
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
empirical_excesses_rpar <- function(data_rain, df_lags, threshold, t0 = 0) {
  # df_lags must contain columns: s1, s2, tau
  excesses <- as.data.table(df_lags[, c("s1", "s2", "tau")])
  excesses[, kij := 0L]
  excesses[, Tobs := 1L]
  s1_name <- df_lags$s1[1]
  ind_s1 <- which(colnames(data_rain) == s1_name)

  X_s1 <- data_rain[, ind_s1]

  for (i in seq_len(nrow(excesses))) {
    # s2_name <- excesses$s2[i]
    tau     <- excesses$tau[i]
    s2_name <- df_lags$s2[i]
    ind_s2 <- which(colnames(data_rain) == s2_name)
    X_s2   <- data_rain[, ind_s2]

    t_shift <- t0 + 1 + tau
    # if (t_shift > nrow(data_rain) || t_shift < 1) {
    #   excesses$kij[i] <- 0
    #   next
    # }

    # Missing rainfall: ignore
    if (is.na(X_s2[t_shift])) {
      excesses$kij[i] <- 0
      next
    }

    # Valid value: evaluate exceedance
    is_excess <- (X_s2[t_shift] >= threshold)
    excesses$kij[i] <- sum(is_excess)
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
#' @param distance The type of spatial norm, "euclidean" or "lalpha".
#'                Default is "euclidean".
#' @param normalize A boolean value to indicate if we want to normalize the
#'                 spatial and temporal lags. Default is FALSE.
#'
#' @return The theoretical chi
#'
#' @export
theoretical_chi <- function(params, df_lags, latlon = FALSE,
                            distance = "euclidean", normalize = FALSE) {
  beta1 <- params[1]
  beta2 <- params[2]
  alpha1 <- params[3]
  alpha2 <- params[4]
  adv <- params[5:6]
  
  chi_df <- df_lags[, c("s1", "s2", "tau",
                    "s1x", "s1y", "s2x", "s2y", "hnorm")]

  # advected coordinates
  chi_df$s1xv <- chi_df$s1x
  chi_df$s1yv <- chi_df$s1y
  # - adv because tau = t - t0 and 
  # s_shifted = s - adv * t, s0_shifted = s0 - adv * t0, h = s - s0
  chi_df$s2xv <- chi_df$s2x - adv[1] * chi_df$tau
  chi_df$s2yv <- chi_df$s2y - adv[2] * chi_df$tau
  chi_df$hx <- chi_df$s2xv - chi_df$s1xv
  chi_df$hy <- chi_df$s2yv - chi_df$s1yv
  chi_df$hnormV <- sqrt(chi_df$hx^2 + chi_df$hy^2)

  if (!normalize) {
    chi_df$h_scaled <- chi_df$hnormV
    chi_df$t_scaled <- abs(chi_df$tau)
    chi_df$hx_scaled <- chi_df$hx
    chi_df$hy_scaled <- chi_df$hy
    # effective parameters
    beta1_eff <- beta1
    beta2_eff <- beta2
  } else {
    # normalization constants
    eps <- 1e-12
    c_h <- median(chi_df$hnormV[chi_df$hnormV > 0], na.rm = TRUE)
    c_tau <- median(abs(chi_df$tau)[chi_df$tau != 0], na.rm = TRUE)
    
    chi_df$h_scaled <- chi_df$hnormV / (c_h + eps)
    chi_df$t_scaled <- abs(chi_df$tau) / (c_tau + eps)

    chi_df$hx_scaled <- chi_df$hx / (c_h + eps)
    chi_df$hy_scaled <- chi_df$hy / (c_h + eps)

    beta1_eff <- beta1 * (c_h^alpha1)
    beta2_eff <- beta2 * (c_tau^alpha2)
  }

  if (distance == "lalpha") {
    # Apply lalpha norm adjustment
    chi_df$vario <- (2 * beta1_eff) * abs(chi_df$hx_scaled)^alpha1 +
                    (2 * beta1_eff) * abs(chi_df$hy_scaled)^alpha1 +
                    (2 * beta2_eff) * abs(chi_df$t_scaled)^alpha2 # same units as V
  } else {
    # Only euclidean distance lag
    chi_df$hlag <-  chi_df$h_scaled
    chi_df$vario <- (2 * beta1_eff) * abs(chi_df$hlag)^alpha1 +
                    (2 * beta2_eff) * abs(chi_df$t_scaled)^alpha2 # same units as V
  }
  # Chi
  chi_df$chi <- 2 * (1 - pnorm(sqrt(0.5 * chi_df$vario)))
  chi_df$chi <- pmin(pmax(chi_df$chi, 1e-12), 1 - 1e-12) # numerical stability

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
#' @param distance The type of spatial norm, "euclidean" or "lalpha".
#'                Default is "euclidean".
#' @param normalize A boolean value to indicate if the lags should be
#'             normalized. Default is FALSE.
#'
#' @return The negative log-likelihood value.
#'
#' @export
neg_ll <- function(params, df_lags, excesses,
                   hmax = NA,  rpar = TRUE, threshold = FALSE,
                   latlon = FALSE, quantile = NA, data = NA,
                   distance = "euclidean", normalize = FALSE) {
  
  if (rpar) {
    p <- 1
  } else {
    if (is.na(data)) stop("Data must be provided for max-stable composite likelihood.")
    if (is.na(quantile)) stop("Quantile or threshold must be provided.")
    Tmax <- nrow(data)
    nmarg <- get_marginal_excess(data, quantile, threshold)
    p <- nmarg / Tmax
  }

  if (!is.na(hmax)) { # filter by hmax
    ind_inf_hmax <- which(df_lags$hnorm < hmax)
    excesses <- excesses[ind_inf_hmax, ]
    df_lags <- df_lags[ind_inf_hmax, ]
  }

  # Get chi values
  chi <- theoretical_chi(params, df_lags, latlon, distance, normalize)
  
  ll_df <- df_lags
  ll_df$kij <- excesses$kij
  ll_df$Tobs <- excesses$Tobs

  # Chi already clipped inside theoretical_chi()
  ll_df$chi <- chi$chi
  
  # # Clip p*chi to avoid log(0)
  eps <- 1e-12
  ll_df$pchi <- pmax(pmin(1 - p * ll_df$chi, 1 - eps), eps)

  ll_df$non_excesses <- ll_df$Tobs - ll_df$kij
  ll_df$ll <- ll_df$kij * log(ll_df$chi) +
              ll_df$non_excesses * log(ll_df$pchi)
  # lambda <- ll_df$chi  # since Tobs = 1
  # ll_df$ll <- ll_df$kij * log(ll_df$chi ) - ll_df$chi  - lfactorial(ll_df$kij)

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
#' @param list_episodes A list of episode data.
#' @param list_lags A list of dataframes with spatial and temporal lag values.
#' @param list_excesses A list of excesses dataframes.
#' @param wind_df A dataframe with wind values. Default is NA.
#' @param hmax The maximum spatial lag value. Default is NA.
#' @param latlon A boolean value to indicate if the locations are in latitude
#'             and longitude. Default is TRUE.
#' @param distance The type of spatial norm, "euclidean" or "lalpha".
#'                  Default is "euclidean".
#' @param threshold A boolean value to indicate if the quantile variable is a
#'               threshold value and not a uniform quantile. Default is FALSE.
#' @param rpar A boolean value to indicate if the data is considered as a
#'          r-Pareto process. Default is TRUE.
#' @param fixed_eta1 A value to fix eta1. Default is NA.
#' @param fixed_eta2 A value to fix eta2. Default is NA.

#' @return The negative log-likelihood value.
#'
#' @export
neg_ll_composite <- function(params, list_episodes, list_excesses,
                             list_lags, wind_df = NA,
                             hmax = NA, latlon = TRUE,
                             distance = "euclidean", threshold = FALSE,
                             rpar = TRUE, normalize = FALSE) {

  if (length(params) == 4) {
    params <- c(params, 0, 0)
  }

  data <- NA
  quantile <- NA

  if (!(distance %in% c("lalpha", "euclidean"))) {
    stop("Parameter 'distance' should be 'lalpha' or 'euclidean'")
  }

  if (all(is.na(wind_df))) {
    adv <- params[5:6]
    adv_df <- matrix(adv, nrow = 1)
  } else {
    eta1 <- params[5]
    eta2 <- params[6]

    vx <- wind_df$vx
    vy <- wind_df$vy

    r <- sqrt(vx^2 + vy^2)
    eps <- 1e-12

    scale <- eta1 * (pmax(r, eps)^(eta2 - 1))  # = eta1*r^(eta2)/r
    adv_x <- scale * vx
    adv_y <- scale * vy

    adv_df <- cbind(adv_x, adv_y)
    if (nrow(adv_df) == 1) adv <- as.vector(adv_df)
  }

  m <- length(list_episodes)
  nll_composite <- 0
  # print(params)

  for (i in seq_len(m)) {
    excesses <- list_excesses[[i]]
    lags <- list_lags[[i]]
    if (!all(is.na(wind_df)) && nrow(adv_df) > 1) {
      adv <- as.vector(adv_df[i, ])
    }

    params_adv <- c(params[1:4], adv)

    if (!rpar) {
      data <- list_episodes[[i]]
      quantile <- quantile
    }

    nll_i <- neg_ll(
      params = params_adv,
      df_lags = lags,
      excesses = excesses,
      hmax = hmax,
      latlon = latlon,
      distance = distance,
      threshold = threshold,
      rpar = rpar,
      data = data,
      quantile = quantile,
      normalize = normalize
    )

    nll_composite <- nll_composite + nll_i
  }

  if (!is.finite(nll_composite)) {
    cat("Non-finite value detected for parameters:", params, "\n")
  }

  return(nll_composite)
}

#' neg_ll_composite_fixed_eta function
#'  
#' Calculate the negative log-likelihood for a list of r-Pareto processes
#' with fixed advection parameters.
#'
#' @param params Vector of variogram parameters (beta1, beta2, alpha1
#'             alpha2).
#' @param list_episodes A list of episode data.
#' @param list_lags A list of dataframes with spatial and temporal lag values.
#' @param list_excesses A list of excesses dataframes.
#' @param wind_df A dataframe with wind values. Default is NA.
#' @param hmax The maximum spatial lag value. Default is NA.
#' @param latlon A boolean value to indicate if the locations are in latitude
#'             and longitude. Default is TRUE.
#' @param distance The type of spatial norm, "euclidean" or "lalpha".
#'                Default is "euclidean".
#' @param threshold A boolean value to indicate if the quantile variable is a
#'               threshold value and not a uniform quantile. Default is FALSE.
#' @param rpar A boolean value to indicate if the data is considered as a
#'          r-Pareto process. Default is TRUE.
#' @param fixed_eta1 A value to fix eta1.
#' @param fixed_eta2 A value to fix eta2.
#' @return The negative log-likelihood value.
#' 
#' @export
neg_ll_composite_fixed_eta <- function(params, list_episodes, list_excesses,
                                       list_lags, wind_df = NA,
                                       hmax = NA, latlon = TRUE,
                                       distance = "euclidean", threshold = FALSE,
                                       rpar = TRUE,
                                       fixed_eta1 = NA, fixed_eta2 = NA) {

  if (!is.na(fixed_eta1) && is.na(fixed_eta2)) {
    # eta1 fixed, eta2 free
    full_params <- c(params[1:4], fixed_eta1, params[5])
  } else if (is.na(fixed_eta1) && !is.na(fixed_eta2)) {
    # eta2 fixed, eta1 free
    full_params <- c(params[1:4], params[5], fixed_eta2)
  } else if (!is.na(fixed_eta1) && !is.na(fixed_eta2)) {
    # both fixed
    full_params <- c(params[1:4], fixed_eta1, fixed_eta2)
  } else {
    # none fixed
    full_params <- params
  }
  print(full_params)
  neg_ll <- neg_ll_composite(
    params = full_params,
    list_lags = list_lags,
    list_episodes = list_episodes,
    list_excesses = list_excesses,
    hmax = hmax,
    wind_df = wind_df,
    latlon = latlon,
    distance = distance,
    threshold = threshold,
    rpar = rpar
  )
  
  print(neg_ll)

  return(neg_ll)
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
#' @param list_simu The list of simulations.
#' @param u The quantile threshold.
#' @param list_lags The list of lags.
#' @param list_excesses The list of excesses.
#' @param init_params The initial parameters for optimization.
#' @param hmax The maximum spatial lag value. Default is NA.
#' @param wind_df The wind dataframe. Default is NA.
#' @param distance The type of spatial norm, "euclidean" or "lalpha". 
#'                 Default is "euclidean".
#' @param hessian A boolean value to indicate if the Hessian matrix should be
#'               computed. Default is FALSE.
#'
#' @return The optimized parameters.
#'
#' @export
process_simulation <- function(i, m, list_simu, u, list_lags,
                               list_excesses,
                               init_params, hmax = NA, wind_df = NA,
                               distance = "euclidean",
                               normalize = FALSE, fixed_eta1 = NA, fixed_eta2 = NA) {
  
  # Output names: always the same
  if (all(!is.na(wind_df))) {
    colnames(wind_df) <- c("vx","vy")
    if (nrow(wind_df) != m) stop("wind_df must have exactly m rows (one per episode).")
    out_names <- c("beta1","beta2","alpha1","alpha2","eta1","eta2")
  } else {
    out_names <- c("beta1","beta2","alpha1","alpha2","adv1","adv2")
  }

  # Bounds
  lower_bounds <- c(1e-8, 1e-8, 1e-8, 1e-8, -Inf, -Inf)
  upper_bounds <- c(Inf, Inf, 1.999, 1.999, Inf, Inf)
  if (all(!is.na(wind_df))) {
    lower_bounds[5:6] <- c(1e-8, 1e-8)
    upper_bounds[5:6] <- c(10, 10)
  }

  # Extract episodes
  list_episodes_m <- list_simu[((i - 1) * m + 1):(i * m)]
  list_lags_m <- list_lags[((i - 1) * m + 1):(i * m)]
  list_excesses_m <- list_excesses[((i - 1) * m + 1):(i * m)]

  if (!is.na(fixed_eta1) && !is.na(fixed_eta2)) {
    lower4 <- lower_bounds[1:4]
    upper4 <- upper_bounds[1:4]

    opt_res <- tryCatch({
      optim(
        par = init_params[1:4],
        fn = neg_ll_composite_fixed_eta,
        list_episodes = list_episodes_m,
        list_lags = list_lags_m,
        list_excesses = list_excesses_m,
        hmax = hmax,
        wind_df = wind_df,          # m×2
        threshold = TRUE,
        latlon = FALSE,
        distance = distance,
        fixed_eta1 = fixed_eta1,
        fixed_eta2 = fixed_eta2,
        method = "L-BFGS-B",
        lower = lower4,
        upper = upper4,
        control = list(maxit = 10000)
      )
    }, error = function(e) return(NULL))

    # opt_res$par est longueur 4 : on reconstruit longueur 6
    if (!is.null(opt_res) && opt_res$convergence == 0) {
      estimates <- c(opt_res$par, fixed_eta1, fixed_eta2)
      names(estimates) <- out_names
      return(as.data.frame(as.list(estimates)))
    } else {
      return(as.data.frame(as.list(setNames(rep(NA,6), out_names))))
    }
  } else {
    opt_res <- tryCatch({
      optim(
        par = init_params,
        fn = neg_ll_composite,
        list_episodes = list_episodes_m,
        list_lags = list_lags_m,
        list_excesses = list_excesses_m,
        hmax = hmax,
        wind_df = wind_df,
        threshold = TRUE,
        latlon = FALSE,
        distance = distance,
        normalize = normalize,
        method = "L-BFGS-B",
        lower = lower_bounds,
        upper = upper_bounds,
        control = list(maxit = 10000)
      )
    }, error = function(e) return(NULL))

      # If failed or no convergence -> return NA row
      if (is.null(opt_res) || opt_res$convergence != 0) {
        return(as.data.frame(as.list(setNames(rep(NA,6), out_names))))
      }
  }
  # Extract estimates
  estimates <- opt_res$par

  # Rename adv->eta if needed
  if (all(!is.na(wind_df))) {
    estimates <- c(estimates[1:4], estimates[5], estimates[6])
    names(estimates) <- out_names
  } else {
    names(estimates) <- out_names
  }

  return(as.data.frame(as.list(estimates)))
}

#' compute_gamma_grid function
#'
#' Computes the theoretical variogram (gamma) grid for given spatial lags, temporal lags,
#' direction, and variogram/advection parameters.
#'
#' @param h_vals Vector of spatial lag values.
#' @param tau_vals Vector of temporal lag values.
#' @param direction Unit vector indicating the direction (length 2).
#' @param theta_hat Vector of variogram parameters: beta1, beta2, alpha1, alpha2.
#' @param eta1 Advection parameter for x direction (default 1).
#' @param eta2 Advection parameter exponent (default 1).
#' @param fictive_v Vector of fictive advection velocities (vx, vy).
#'
#' @return Data frame with columns: h, h_norm, h_normV, tau, tau_min, gamma.
#' @export
compute_gamma_grid <- function(h_vals, tau_vals, direction,
                               theta_hat, eta1 = 1, eta2 = 1, fictive_v = c(1, 1)) {
  
  df_out <- data.frame()
  
  # Advection velocity components
  vx <- eta1 * abs(fictive_v[1])^eta2 * sign(fictive_v[1])
  vy <- eta1 * abs(fictive_v[2])^eta2 * sign(fictive_v[2])
  
  for (tau in tau_vals) {
    for (h in h_vals) {
      
      beta1  <- theta_hat[1]
      beta2  <- theta_hat[2]
      alpha1 <- theta_hat[3]
      alpha2 <- theta_hat[4]
      
      if (is.null(direction) || identical(direction, "none")) {
        
        h_norm <- h
        hx <- h        # artificial components for consistency
        hy <- 0
        
      } else {
        # Standard: projection of h along direction
        hx <- h * direction[1]
        hy <- h * direction[2]
        h_norm <- sqrt(hx^2 + hy^2)
      }
      
      term1 <- beta1 * abs(hx - vx * tau)^alpha1 +
              beta1 * abs(hy - vy * tau)^alpha1
      
      term2 <- beta2 * tau^alpha2
      
      gamma_hat <- term1 + term2
      
      df_out <- rbind(df_out, data.frame(
        h = h,
        h_norm = h_norm,
        h_normV = sqrt((hx - vx * tau)^2 + (hy - vy * tau)^2),
        tau = tau,
        tau_min = tau * 60,
        gamma = gamma_hat
      ))
    }
  }
  return(df_out)
}

#' convert_params function
#'
#' Convert variogram parameters based on scaling constants for space and time.
#'
#' @param beta1 The spatial scale parameter.
#' @param beta2 The temporal scale parameter.
#' @param alpha1 The spatial smoothness parameter.
#' @param alpha2 The temporal smoothness parameter.
#' @param c_x The spatial scaling constant. Default is 1.
#' @param c_t The temporal scaling constant. Default is 1.
#' @return A list containing the converted parameters beta1 and beta2.
#' @export
convert_params <- function(beta1, beta2, alpha1, alpha2, c_x = 1, c_t = 1) {
  beta1_new <- beta1 / (c_x^alpha1)
  beta2_new <- beta2 / (c_t^alpha2)
  list(beta1 = beta1_new, beta2 = beta2_new)
}