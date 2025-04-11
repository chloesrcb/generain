# EXTREME EPISODES -------------------------------------------------------------


#' select_extreme_episodes function
#'
#' This function selects extreme episodes based on the given data and a 
#' quantile. A spatio-temporal specific neighborhood is considered to avoid
#' selecting nearby episodes. Each episode is selected according to a 
#' r-Pareto process with r(X) = X(s0, t0). Each episode will have a size of
#' episode_size.
#'
#' @param sites_coords The coordinates of the sites.
#' @param data The rainfall dataframe.
#' @param quantile The quantile value.
#' @param min_spatial_dist The minimum spatial distance between two episodes.
#' @param episode_size The temporal window size, i.e. the size of the episode.
#' @param n_max_episodes The maximum number of episodes to select. Default is 10000.
#' @param time_ext The temporal window extension. Default is 0.
#' @param latlon A boolean value to indicate if the locations are in latitude and longitude. Default is TRUE.
#' @return A list of selected points (s0, t0) and the corresponding threshold value u_s0 for each episode.
#'
#' @import data.table
#' @import geosphere
#'
#' @export
select_extreme_episodes <- function(sites_coords, data, quantile,
                                    min_spatial_dist, episode_size,
                                    n_max_episodes = 10000, time_ext = 0,
                                    latlon = TRUE) {
  # Convert data to a matrix for fast access
  data <- as.matrix(data)
  site_names <- colnames(data)

  # Compute threshold for each site
  thresholds_by_site <- apply(data, 2, function(col) quantile(col, probs = quantile, na.rm = TRUE))

  # Compute distance matrix (ensure named rows/columns)
  if (latlon){
    dist_matrix <- as.matrix(distm(sites_coords[site_names, ], fun = distHaversine)) / 1000
  } else {
    dist_matrix <- as.matrix(dist(sites_coords[site_names, ])) / 1000
  }
  rownames(dist_matrix) <- site_names
  colnames(dist_matrix) <- site_names

  # Initialize
  selected_points <- data.table(s0 = character(), t0 = integer(), u_s0 = numeric())
  invalid_time_mask <- matrix(FALSE, nrow = nrow(data), ncol = ncol(data))
  nb_episode <- 0

  while (nb_episode < n_max_episodes) {
    # Identify all (s0, t0) where data exceeds threshold
    exceed_mask <- sweep(data, 2, thresholds_by_site, FUN = ">")
    exceed_mask[invalid_time_mask] <- FALSE
    exceed_indices <- which(exceed_mask, arr.ind = TRUE)

    if (nrow(exceed_indices) == 0) break

    # Sort by time to get first exceeding time points
    exceed_indices <- exceed_indices[order(exceed_indices[, 1]), , drop = FALSE]

    best_candidate <- NULL

    for (idx in seq_len(nrow(exceed_indices))) {
      t0_candidate <- exceed_indices[idx, 1]
      s0_col_index <- exceed_indices[idx, 2]
      s0_candidate <- site_names[s0_col_index]
      u_s0 <- thresholds_by_site[s0_candidate]

      if (nrow(selected_points) == 0) {
        best_candidate <- data.table(s0 = s0_candidate, t0 = t0_candidate, u_s0 = u_s0)
        nb_episode <- nb_episode + 1
        break
      } else {
        selected_sites <- selected_points$s0
        time_differences <- abs(selected_points$t0 - t0_candidate)
        distances <- dist_matrix[selected_sites, s0_candidate, drop = FALSE]
        valid_pairs <- (distances >= min_spatial_dist) | (time_differences > episode_size + 1)

        if (all(valid_pairs)) {
          best_candidate <- data.table(s0 = s0_candidate, t0 = t0_candidate, u_s0 = u_s0)
          nb_episode <- nb_episode + 1
          break
        } else {
          invalid_time_mask[t0_candidate, s0_col_index] <- TRUE
        }
      }
    }

    if (!is.null(best_candidate)) {
      # Store the selected episode
      selected_points <- rbindlist(list(selected_points, best_candidate))
      t_inf <- max(1, best_candidate$t0 - time_ext)
      t_sup <- min(nrow(data), best_candidate$t0 + episode_size + time_ext)
      invalid_time_mask[t_inf:t_sup, which(site_names == best_candidate$s0)] <- TRUE
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
    t_inf <- t0 - beta
    t_sup <- t0 + episode_size + beta

    # Check that the episode is the correct size (episode_size)
    episode_size_test <- t_sup - t_inf + 1
    if (episode_size_test == episode_size + 1 + 2 * beta) {
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

#' empirical_excesses_rpar function
#'
#' This function calculates the empirical excesses based on indicators above a
#' quantile threshold for the r-Pareto process.
#'
#' @param data_rain The rainfall data.
#' @param quantile The quantile value.
#' @param df_lags The dataframe with spatial and temporal lag values.
#' @param threshold A boolean value to indicate if the quantile variable is a
#'                threshold value and not a uniform quantile. Default is FALSE.
#' @param t0 The starting time. Default is 1.
#'
#' @return The empirical excesses dataframe with the number of excesses kij and
#' the number of possible excesses Tobs, with the lag values.
#'
#' @export
empirical_excesses_rpar <- function(data_rain, quantile, df_lags,
                                    threshold = FALSE, t0 = 0) {
  excesses <- df_lags # copy the dataframe
  unique_tau <- unique(df_lags$tau) # unique temporal lags
  ind_s1 <- df_lags$s1[1] # s0

  for (tau in unique_tau) { # loop over temporal lags
    df_h_t <- df_lags[df_lags$tau == tau, ] # get the dataframe for each tau lag

    for (i in seq_len(nrow(df_h_t))) { # loop over each pair of sites
      # get the indices of the sites
      ind_s2 <- as.numeric(as.character(df_h_t$s2[i]))

      # get the data for the second site
      X_s2 <- data_rain[, c(ind_s2), drop = FALSE]
      X_s2 <- as.vector(na.omit(X_s2))

      # if (!threshold) {
      #   X_s2 <- rank(X_s2) / (length(X_s2) + 1) # uniformize the data
      # }

      # shifted data
      X_s_t <- X_s2[t0 + 1 + tau] # X_{s,t0 + tau}
      nmargin <- sum(X_s_t > quantile) # 0 or 1

      # store the number of excesses and T - tau
      excesses$Tobs[excesses$s1 == ind_s1
                      & excesses$s2 == ind_s2
                      & excesses$tau == tau] <- 1

      excesses$kij[excesses$s1 == ind_s1
                    & excesses$s2 == ind_s2
                    & excesses$tau == tau] <- nmargin
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
#'             and longitude.
#' @param directional A boolean value to indicate if the variogram is
#'                  directional.
#'
#' @return The theoretical chi
#'
#' @export
theoretical_chi <- function(params, df_lags, latlon,
                             directional) {
  beta1 <- params[1]
  beta2 <- params[2]
  alpha1 <- params[3]
  alpha2 <- params[4]
  adv <- params[5:6]
  chi_df <- df_lags[c("s1", "s2", "tau", "s1x", "s1y", "s2x", "s2y", "hnorm")]

  if (latlon) {
    # Compute spatial distance in km (already in km in df_lags)
    haversine_df <- haversine_distance_with_advection(chi_df$s1y, chi_df$s1x,
                                      chi_df$s2y, chi_df$s2x, adv,
                                      (chi_df$tau)* 3600) # seconds
    chi_df$hnormV <- haversine_df$distance # km
    chi_df$theta <- haversine_df$theta

  } else {
   # Cartesian coordinates case
   chi_df$s1xv <- chi_df$s1x
   chi_df$s1yv <- chi_df$s1y
   # tau = t - t0 then t = tau + t0
   # + adv because tau = t0 - t (in simulation)
   chi_df$s2xv <- chi_df$s2x + adv[1] * chi_df$tau # m/s * s = m
   chi_df$s2yv <- chi_df$s2y + adv[2] * chi_df$tau

   # Recompute distance and angle regardless of advection
   chi_df$hnormV <- sqrt((chi_df$s2xv - chi_df$s1xv)^2 +
                        (chi_df$s2yv - chi_df$s1yv)^2) # m
   # Convert to km
   chi_df$hnormV <- chi_df$hnormV # km

   chi_df$theta <- atan2(chi_df$s2yv - chi_df$s1yv, chi_df$s2xv - chi_df$s1xv)

  }

  if (directional) {
      # Apply directional distance adjustment
      chi_df$x_polar <- chi_df$hnormV * cos(chi_df$theta)
      chi_df$y_polar <- chi_df$hnormV * sin(chi_df$theta)
      chi_df$hlag <- chi_df$x_polar + chi_df$y_polar
  } else {
      # Only distance lag
      chi_df$hlag <-  chi_df$hnormV
  }

  # Compute variogram and chi
  chi_df$vario <- (2 * beta1) * abs(chi_df$hlag)^alpha1 +
                  (2 * beta2) * abs(chi_df$tau)^alpha2 # km and hours

  chi_df$chi <- 2 * (1 - pnorm(sqrt(0.5 * chi_df$vario)))
  # print only chi > 0
  # chi_no0 <- chi_df[chi_df$chi > 0 & chi_df$chi < 1, ]
  # print(chi_no0$chi)
  return(chi_df)
}

theoretical_chi <- function(params, df_lags, latlon=F, directional=F) {
  beta1 <- params[1]
  beta2 <- params[2]
  alpha1 <- params[3]
  alpha2 <- params[4]
  adv <- params[5:6]
  chi_df <- df_lags[c("s1", "s2", "tau", "s1x", "s1y", "s2x", "s2y")]

  # Advected coordinates
  chi_df$s1x_shifted <- chi_df$s1x
  chi_df$s1y_shifted <- chi_df$s1y
  chi_df$s2x_shifted <- chi_df$s2x + adv[1] * chi_df$tau
  chi_df$s2y_shifted <- chi_df$s2y + adv[2] * chi_df$tau

  # Euclidean distance
  chi_df$hnormV <- sqrt((chi_df$s2x_shifted - chi_df$s1x_shifted)^2 +
                        (chi_df$s2y_shifted - chi_df$s1y_shifted)^2)

  # Compute variogram and chi
  chi_df$vario <- (2 * beta1) * abs(chi_df$hnormV)^alpha1 +
                  (2 * beta2) * abs(chi_df$tau)^alpha2

  chi_df$chi <- 2 * (1 - pnorm(sqrt(0.5 * chi_df$vario)))

  return(chi_df)
}


neg_ll <- function(params, data, df_lags, quantile, excesses, hmax = NA,
                  s0 = NA, t0 = NA, threshold = FALSE) {
  Tmax <- nrow(data) # number of total observations
  # print(params)
  if (all(!is.na(s0))) { # if we have a conditioning location
    p <- 1 # sure excess for r-Pareto process in (s0,t0)
  } else {
    # number of marginal excesses
    nmarg <- get_marginal_excess(data, quantile, threshold)
    p <- nmarg / Tmax # probability of marginal excesses
  }

  # # Bounds for the parameters
  # lower.bound <- c(1e-6, 1e-6, 1e-6, 1e-6)
  # upper.bound <- c(Inf, Inf, 1.999, 1.999)
  # if (length(params) == 6) {
  #   lower.bound <- c(lower.bound, -Inf, -Inf)
  #   upper.bound <- c(upper.bound, Inf, Inf)
  # }

  # # Check if the parameters are in the bounds
  # if (any(params < lower.bound) || any(params > upper.bound)) {
  #   return(1e50)
  # }

  chi <- theoretical_chi(params, df_lags) # get chi matrix
  ll_df <- df_lags # copy the dataframe
  ll_df$kij <- excesses$kij # number of excesses
  ll_df$Tobs <- excesses$Tobs
  ll_df$hnormV <- chi$hnormV
  ll_df$chi <- chi$chi
  ll_df$chi <- ifelse(ll_df$chi <= 0, 1e-10, ll_df$chi)
  ll_df$pchi <- 1 - p * ll_df$chi
  ll_df$pchi <- ifelse(ll_df$pchi <= 0, 1e-10, ll_df$pchi)

  # number of non-excesses
  ll_df$non_excesses <- ll_df$Tobs - ll_df$kij
  ll_df$ll <- ll_df$kij * log(ll_df$chi) +
              ll_df$non_excesses * log(ll_df$pchi)
  if (!is.na(hmax)) {
    ll_df <- ll_df[ll_df$hnorm <= hmax, ]
  }

  nll <- -sum(ll_df$ll, na.rm = TRUE)
  return(nll)
}


neg_ll_composite_simu <- function(params, list_simu, df_lags, quantile,
                    list_excesses, hmax = NA, s0 = NA,
                    t0 = NA, threshold = FALSE) {
  print(params)
  # Bounds for the parameters
  lower.bound <- c(1e-6, 1e-6, 1e-6, 1e-6)
  upper.bound <- c(Inf, Inf, 1.999, 1.999)
  if (length(params) == 6) {
    lower.bound <- c(lower.bound, -Inf, -Inf)
    upper.bound <- c(upper.bound, Inf, Inf)
  }

  # Check if the parameters are in the bounds
  if (any(params < lower.bound) || any(params > upper.bound)) {
    return(1e50)
  }

  m <- length(list_simu) # number of replicates
  nll_composite <- 0 # composite negative log-likelihood
  for (i in 1:m) {
    # extract simulation data from i-th simulation
    simu <- list_simu[[i]]
    excesses <- list_excesses[[i]]
    nll_i <- neg_ll(params, simu, df_lags, quantile, hmax = hmax,
                    excesses = excesses, s0 = s0, t0 = t0,
                    threshold = threshold)
    nll_composite <- nll_composite + nll_i
  }
  return(nll_composite)
}



# theoretical_chi <- function(params, df_lags, latlon,
#                                             directional) {
#   beta1 <- params[1]
#   beta2 <- params[2]
#   alpha1 <- params[3]
#   alpha2 <- params[4]
#   adv <- params[5:6]
#   chi_df <- df_lags[c("s1", "s2", "tau", "s1x", "s1y", "s2x", "s2y", "hnorm")]

#   # Convert tau from hours to seconds (assuming tau is in hours)
#   # chi_df$tau_seconds <- chi_df$tau * 3600  # Convert tau from hours to seconds

#   if (latlon) {
#     # Compute spatial distance in km (already in km in df_lags)
#     haversine_df <- haversine_distance_with_advection(chi_df$s1y, chi_df$s1x,
#                                       chi_df$s2y, chi_df$s2x, adv,
#                                       chi_df$tau)  # Use tau in seconds
#     chi_df$hnormV <- haversine_df$distance

#     # Directional adjustment: compute the angle (theta) for the directional
#     # variogram
#     chi_df$theta <- haversine_df$theta

#   } else {
#     # Cartesian coordinates case
#     chi_df$s1xv <- chi_df$s1x
#     chi_df$s1yv <- chi_df$s1y
#     chi_df$s2xv <- chi_df$s2x + adv[1] * chi_df$tau
#     chi_df$s2yv <- chi_df$s2y + adv[2] * chi_df$tau
#     if (all(adv == 0)) {
#       chi_df$hnormV <- chi_df$hnorm
#     } else {
#       chi_df$hnormV <- sqrt((chi_df$s2xv - chi_df$s1xv)^2 +
#                             (chi_df$s2yv - chi_df$s1yv)^2)
#     }
#     # Directional adjustment: compute the angle (theta) for the directional
#     # variogram
#     chi_df$theta <- atan2(chi_df$s2yv - chi_df$s1yv, chi_df$s2xv - chi_df$s1xv)
#   }

#   if (directional) {
#       # Apply directional distance adjustment
#       chi_df$x_polar <- chi_df$hnormV * cos(chi_df$theta)
#       chi_df$y_polar <- chi_df$hnormV * sin(chi_df$theta)
#       chi_df$hlag <- chi_df$x_polar + chi_df$y_polar
#   } else {
#       # Only distance lag
#       chi_df$hlag <-  chi_df$hnormV
#   }

#   # Compute variogram and chi
#   chi_df$vario <- (2 * beta1) * abs(chi_df$hlag)^alpha1 +
#                   (2 * beta2) * abs(chi_df$tau)^alpha2 # km and hours

#   chi_df$chi <- 2 * (1 - pnorm(sqrt(0.5 * chi_df$vario)))

#   return(chi_df)
# }



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
#'
#' @return The negative log-likelihood value.
#'
#' @export
neg_ll <- function(params, df_lags, excesses,
                      hmax = NA,  rpar = TRUE, threshold = FALSE,
                      latlon = FALSE, quantile = NA, data = NA,
                      directional = TRUE) {
  # if (rpar) { # if we have a conditioning location
    p <- 1 # sure excess for r-Pareto process in (s0,t0)
  # } else {
  #   if (is.na(data)) {
  #     stop("The data must be provided for the max-stable composite likelihood.")
  #   } else if (is.na(quantile)) {
  #     stop("The quantile or threshold must be provided for the
  #           max-stable composite likelihood.")
  #   }
  #   Tmax <- nrow(data) # number of total observations
  #   # number of marginal excesses
  #   nmarg <- get_marginal_excess(data, quantile, threshold)
  #   p <- nmarg / Tmax # probability of marginal excesses
  # }

  # compute theoretical chi values
  chi <- theoretical_chi(params, df_lags, latlon, directional)
  ll_df <- df_lags # copy the dataframe
  ll_df$kij <- excesses$kij # number of excesses
  ll_df$Tobs <- excesses$Tobs
  ll_df$hnormV <- chi$hnormV
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
  if (!is.na(hmax)) {
    ll_df <- ll_df[ll_df$hnorm <= hmax, ]
  }

  nll <- -sum(ll_df$ll, na.rm = TRUE)
  return(nll)
}

#' neg_ll_composite_simu function
#'
#' Calculate the negative log-likelihood for a list of simulations.
#'
#' @param params Vector of variogram parameters (beta1, beta2, alpha1, alpha2).
#' @param list_simu A list of simulated data.
#' @param df_lags The dataframe with spatial and temporal lag values.
#' @param quantile The quantile value.
#' @param list_excesses A list of excesses dataframes.
#' @param hmax The maximum spatial lag value.
#' @param threshold A boolean value to indicate if the quantile variable is a
#'               threshold value and not a uniform quantile. Default is FALSE.
#' @param directional A boolean value to indicate if the variogram is
#'                   directional. Default is TRUE.
#' @param rpar A boolean value to indicate if the data is considered as an
#'           r-Pareto process. Default is TRUE.
#'
#' @return The negative log-likelihood value.
#'
#' @import stats
#'
#' @export
neg_ll_composite_simu <- function(params, list_simu, df_lags, quantile,
                    list_excesses, adv=NA, hmax = NA, threshold = FALSE,
                    directional = FALSE, rpar = TRUE) {
  m <- length(list_simu) # number of replicates
  nll_composite <- 0 # composite negative log-likelihood
  print(params)
  # if (length(adv) == 1 && !is.na(adv)) {
  #   eta1 <- params[5]
  #   eta2 <- params[6]
  #   adv_x <- (abs(adv[1])^eta1) * sign(adv[1]) * eta2
  #   adv_y <- (abs(adv[2])^eta1) * sign(adv[2]) * eta2
  #   adv <- cbind(adv_x, adv_y)
  #   # params <- c(params, adv) # add advection parameters
  # }
  
  for (i in 1:m) {
    # extract simulation data from i-th simulation
    simu <- list_simu[[i]]
    excesses <- list_excesses[[i]]
    # nll_i <- neg_ll(params, data = simu, df_lags = df_lags,
    #                 quantile = quantile, hmax = hmax,
    #                 excesses = excesses, threshold = threshold,
    #                 directional = directional,
    #                 rpar = rpar)
    nll_i <- neg_ll(params = params,
                     df_lags = df_lags,
                     hmax = hmax, excesses = excesses,
                     directional = directional,
                     threshold = threshold, rpar = rpar,
                     data = data, quantile = quantile)
    nll_composite <- nll_composite + nll_i
  }
  return(nll_composite)
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
#'
#' @return The negative log-likelihood value.
#'
#' @export
neg_ll_composite <- function(params, list_episodes, list_excesses,
                             list_lags, wind_df = NA,
                             hmax = NA, latlon = TRUE,
                             directional = TRUE, threshold = FALSE,
                             rpar = TRUE) {
  # params[1:4] <-  c(0.01, 0.2, 1.5, 1)
  # Add default values for advection parameters
  if (length(params) == 4) {
    params <- c(params, 0, 0) # No advection
  }
  if(rpar) {
    data <- NA
    quantile <- NA
  }
  if (all(is.na(wind_df))) {
    adv <- params[5:6]
  } else if (length(wind_df) == 2) {
    eta1 <- params[5]
    eta2 <- params[6]
    adv_x <- (abs(wind_df[1])^eta1) * sign(wind_df[1]) * eta2
    adv_y <- (abs(wind_df[2])^eta1) * sign(wind_df[2]) * eta2
    adv <- cbind(adv_x, adv_y)
  } else {
    eta1 <- params[5]
    eta2 <- params[6]
    adv_x <- (abs(wind_df$vx)^eta1) * sign(wind_df$vx) * eta2
    adv_y <- (abs(wind_df$vy)^eta1) * sign(wind_df$vy) * eta2
    adv_df <- cbind(adv_x, adv_y)
  }
  
  print(params)
  # # Ensure list_lags has the same length as list_episodes
  # if (is.data.frame(list_lags)) {
  #   list_lags <- list(list_lags)
  #   list_lags <- rep(list_lags, length(list_episodes))
  # }
  m <- length(list_episodes) # Number of r-Pareto processes
  nll_composite <- 0
  for (i in 1:m) {
    # extract episode and excesses from i-th r-pareto process from data
    excesses <- list_excesses[[i]]
    lags <- list_lags[[i]]

    if (!all(is.na(wind_df)) && length(wind_df) != 2) {
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
convert_to_cardinal <- function(degrees) {
  if (is.na(degrees)) {
    return(NA)  # Return NA if the input is NA
  }

  directions <- c("N", "NE", "E", "SE", "S", "SW", "W", "NW", "N")
  breaks <- c(0, 22.5, 67.5, 112.5, 157.5, 202.5, 247.5, 292.5, 337.5, 360)

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
cardinal_to_degree <- function(cardinal) {
  if (is.na(cardinal)) {
    return(NA)  # Return NA if the input is NA
  }

  directions <- c("N", "NE", "E", "SE", "S", "SW", "W", "NW")
  degrees <- c(0, 45, 90, 135, 180, 225, 270, 315)

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

  # Time when there is excess above quantile in s0
  # s0_excess_time <- which(episode[, s0] > u)
  excess_idx <- which(episode[[s0]] > u)
  excess_times <- timestamps[excess_idx]
  wind_excess_subset <- wind_df %>%
    dplyr::filter(datetime %in% excess_times)
  if (nrow(wind_excess_subset) > 0) {
    max_idx <- which.max(wind_excess_subset$FF)
    DD_excess <- wind_excess_subset$DD[max_idx]
    cardDir_excess <- wind_excess_subset$cardDir[max_idx]
  } else {
    DD_excess <- NA
    cardDir_excess <- NA
  }
  # Compute the mean wind speed and direction for the episode if there is data
  if (nrow(wind_subset) > 0) {
    FF <- wind_subset$FF[speed_time + 1]
    cardDir <- get_mode_dir(wind_subset$cardDir)
    DD <- mean(wind_subset$DD[wind_subset$cardDir == cardDir])
    DD_deg <- cardinal_to_degree(cardDir)
    # cardDir_excess <- get_mode_dir(wind_subset_excess$cardDir)
    # DD_excess <- mean(wind_subset_excess$DD[wind_subset_excess$cardDir ==
    #                                                         cardDir_excess])
    # DD_excess <- mean(wind_subset_excess$DD)
    DD_t0 <- wind_subset$DD[speed_time + 1]
    cardDir_t0 <- wind_subset$cardDir[speed_time + 1]
  } else {
    FF <- NA
    DD <- NA
    DD_deg <- NA
    cardDir <- NA
    # cardDir_excess <- NA
    # DD_excess <- NA
    DD_t0 <- NA
    cardDir_t0 <- NA
  }

  return(data.frame(FF = FF, DD = DD, DD_deg = DD_deg, cardDir = cardDir,
                    DD_t0 = DD_t0, cardDir_t0 = cardDir_t0,
                    DD_excess = DD_excess, cardDir_excess = cardDir_excess))
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

#' evaluate_optim function
#'
#' This function evaluates the optimization process on simulated data (for
#' example, following a Brown-Resnick process) using the negative
#' log-likelihood function.
#'
#' @param list_simu A list of simulated data
#' @param quantile The quantile value
#' @param true_param The true variogram parameter (beta1, beta2, alpha1, alpha2)
#' @param df_lags The dataframe with spatial and temporal lag values
#' @param locations The locations dataframe
#' @param parscale The scaling parameter for the optimization process. Default
#'                 is c(1, 1, 1, 1).
#' @param latlon A boolean value to indicate if the locations are in latitude
#'              and longitude. Default is FALSE.
#'
#' @return The result of the optimization process as a dataframe.
#'
#' @import spam
#' @import stats
#' @import doParallel
#' @import foreach
#'
#' @export
evaluate_optim <- function(list_simu, quantile, true_param, df_lags,
                           locations, parscale = c(1, 1, 1, 1),
                           latlon = FALSE) {

  if (length(true_param) == 6 && length(parscale) == 4) {
    parscale <- c(parscale, 1, 1)
  }

  n_res <- length(list_simu)
  df_result <- data.frame(beta1 = rep(NA, n_res), beta2 = rep(NA, n_res),
                          alpha1 = rep(NA, n_res), alpha2 = rep(NA, n_res))

  if (length(true_param) == 6) {
    df_result$adv1 <- rep(NA, n_res)
    df_result$adv2 <- rep(NA, n_res)
  }

  count_cv <- 0

  cl <- parallel::makeCluster(parallel::detectCores() - 1)
  # clusterEvalQ(cl, library(generain))
  parallel::clusterExport(cl, c("list_simu", "neg_ll", "true_param",
                      "quantile", "df_lags", "locations", "nmin",
                      "parscale", "latlon", "empirical_excesses",
                      "optim"), envir = NULL)

  results <- parallel::parLapply(cl, 1:n_res, function(n) {
    simu_df <- as.data.frame(list_simu[[n]])

    if (length(true_param) == 6) {
      excesses <- NULL
    } else {
      excesses <- empirical_excesses(simu_df, quantile, df_lags)
    }

    result <- tryCatch({
      optim(par = true_param, fn = neg_ll, excesses = excesses,
            quantile = quantile,
            df_lags = df_lags, locations = locations,
            simu = simu_df,
            method = "CG", control = list(parscale = parscale, maxit = 10000))
    }, error = function(e) {
      NULL
    })

    if (!is.null(result) && result$convergence == 0) {
      params <- result$par
      return(c(params, TRUE))
    } else {
      return(rep(NA, length(true_param) + 1))
    }
  })

  parallel::stopCluster(cl)

  for (n in 1:n_res) {
    res <- results[[n]]
    if (!is.na(res[length(res)])) {
      count_cv <- count_cv + 1
      df_result[n, 1:length(true_param)] <- res[1:length(true_param)]
    }
  }

  print(paste0("Number of convergences: ", count_cv))
  return(df_result)
}

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
#' @param M The number of simulations.
#' @param m The number of simulations to process.
#' @param list_simuM The list of simulations.
#' @param u The quantile value.
#' @param df_lags The dataframe with spatial and temporal lag values.
#' @param t0 The starting time.
#' @param true_param The true variogram parameter (beta1, beta2, alpha1, alpha2)
#' @param hmax The maximum spatial lag value. Default is NA.
#'
#' @return The optimized parameters.
#'
#' @export
process_simulation <- function(i, M, m, list_simuM, u, df_lags, t0,
                              true_param, hmax = NA, wind_df = NA) {
  # Get the m corresponding simulations from list_simu inside a list
  list_episodes <- list_simuM[((i - 1) * m + 1):(i * m)]
  # Compute excesses
  list_excesses <- lapply(list_episodes, function(episode) {
    empirical_excesses_rpar(episode, u, df_lags, threshold = TRUE, t0 = t0)
  })

  list_lags <- lapply(1:m, function(i) {
   df_lags
})

  # Optimize
  result <- optim(
    par = true_param,
    fn = neg_ll_composite,
    list_episodes = list_episodes,
    list_lags = list_lags,
    list_excesses = list_excesses,
    hmax = hmax,
    wind_df = wind_df,
    threshold = TRUE,
    latlon = FALSE,
    directional = TRUE,
    method = "L-BFGS-B",
    lower = c(1e-8, 1e-8, 1e-8, 1e-8, -Inf, -Inf),
    upper = c(Inf, Inf, 1.999, 1.999, Inf, Inf),
    control = list(maxit = 10000)
  )

  estimates <- result$par
  if (result$convergence != 0) {
    estimates <- rep(NA, length(true_param))
  }
  
  return(estimates)
}
