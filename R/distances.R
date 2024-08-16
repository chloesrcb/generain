#' This function reshapes a distance matrix into a long dataframe.
#'
#' @param locations A matrix or data frame containing the coordinates of the
#' locations.
#' @param dmax The maximum distance threshold if specified.
#' @param latlon Logical indicating whether the coordinates are in latitude and
#' longitude format. Default is TRUE.
#' @param adv A vector of advection values. Default is c(0, 0).
#' @param tau A vector of temporal lags. Default is 1:10.
#' @param alpha_spa The spatial lag exponent. Default is 1.5.
#' @return A distance matrix containing the reshaped distances.
#' 
#' @import geodist
#' 
#' @export
get_dist_mat <- function(locations, dmax = NA, latlon = TRUE,
                                   adv = c(0, 0), tau = 1:10, alpha_spa = 1.5) {
  # get longitude and latitude in a dataframe to get distance between points
  loc <- data.frame(lat = locations$Latitude, lon = locations$Longitude)

  # without advection
  if (latlon) {
    dist_mat <- geodist(loc, measure = "haversine") # meters by default
  } else {
    dist_mat <- as.matrix(dist(loc))
  }

  if (!is.na(dmax)) { # threshold
    dist_mat[dist_mat > dmax] <- 0
  }

  if (all(adv == c(0, 0))) {
    return(dist_mat)
  } else { # with advection
    dist_mats_adv <- list()
    for (t in tau) {
      adv_t <- adv * t
      dist_mat_adv <- matrix(0, nrow = nrow(dist_mat), ncol = ncol(dist_mat))

      for (i in 1:nrow(dist_mat)) {
        for (j in 1:ncol(dist_mat)) {
          h <- c(locations$Longitude[j] - locations$Longitude[i],
                 locations$Latitude[j] - locations$Latitude[i])
          dist_mat_adv[i, j] <- sqrt(sum((h - adv_t)^2))
          # dist_mat_adv[i, j] <- norm_Lp(h[1] - adv_t[1], h[2] - adv_t[2], alpha_spa)
        }
      }

      dist_mats_adv[[paste0("t", t)]] <- dist_mat_adv
    }

    return(dist_mats_adv)
  }
}

#' This function reshapes a distance matrix into a long dataframe.
#'
#' @param dist_mat The input distance matrix.
#' @return A dataframe containing the reshaped distances.
#'
#' @import reshape2
#'
#' @export
reshape_distances <- function(dist_mat) {
  # convert distance matrix into a dataframe
  df_dist <- as.data.frame(dist_mat)
  n <- nrow(df_dist)
  colnames(df_dist) <- c(1:n)
  rownames(df_dist) <- c(1:n)

  # Make a triangle
  df_dist[lower.tri(df_dist)] <- NA

  # Convert to a data frame, and add tenure labels
  df_dist <- as.data.frame(df_dist)
  df_dist$Y <- 1:n

  # Reshape to suit ggplot, remove NAs, and sort the labels
  df_dist <- na.omit(reshape2::melt(df_dist, "Y", variable_name = "X"))
  colnames(df_dist) <- c("Y", "X", "value")
  df_dist$X <- factor(df_dist$X, levels = rev(levels(df_dist$X)))

  return(df_dist)
}


#' get_h_vect function
#'
#' This function calculates the spatial lag vector based on the given distance
#' dataframe and maximum spacial lag value.
#'
#' @param df_dist The distance dataframe or matrix for intervals.
#' @param hmax The maximum spacial lag value. Default is NA.
#' @param intervals A logical value indicating whether we work with intervals.
#' @return The spatial distance lags vector or list of lag vectors.
#'
#' @export
get_h_vect <- function(dist, hmax = NA, intervals = FALSE) {
  # get unique distances
  if (intervals) {
    dist_int <- na.omit(as.vector(dist))
    h_vect <- unique(dist_int)
  } else {
    h_vect <- sort(dist$value)
    h_vect <- unique(h_vect)
    if (!is.na(hmax)) {
        h_vect <- h_vect[h_vect <= hmax]
    }
  }
  h_vect <- h_vect[h_vect > 0]
  return(h_vect)
}


#' norm_Lp function
#'
#' This function calculates the Lp norm of two vectors.
#'
#' @param x The first vector.
#' @param y The second vector.
#' @param p The exponent value.
#' @return The Lp norm of the two vectors.
#'
#' @export
norm_Lp <- function(x, y, p) {
  return((abs(x)^p + abs(y)^p)^(1 / p))
}

#' get_lag_vectors function
#'
#' This function calculates the lag vectors between pairs of points.
#'
#' @param df_coords A dataframe containing the coordinates of the points.
#' @param params A vector of parameters.
#' @param hmax The maximum distance threshold. Default is NA.
#' @param tau_vect A vector of temporal lags. Default is 1:10.
#'
#' @return A dataframe containing the lag vectors.
#'
#' @export
get_lag_vectors <- function(df_coords, params, hmax = NA, tau_vect = 1:10) {
  alpha_spa <- params[3]

  if (length(params) != 6) {
    adv <- c(0, 0)
  } else {
    adv <- params[5:6]
  }

  n <- nrow(df_coords)
  tau_len <- length(tau_vect)
  total_combinations <- (n * (n - 1) / 2) * tau_len

  # Preallocate the dataframe
  lags <- data.frame(
    s1 = integer(total_combinations),
    s2 = integer(total_combinations),
    h1 = numeric(total_combinations),
    h2 = numeric(total_combinations),
    tau = integer(total_combinations),
    hnorm = numeric(total_combinations)
  )

  idx <- 1
  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      lag_latitude <- df_coords$Latitude[j] - df_coords$Latitude[i]
      lag_longitude <- df_coords$Longitude[j] - df_coords$Longitude[i]

      # hnorm <- norm_Lp(lag_latitude, lag_longitude, alpha_spa)
      hnorm <- sqrt(lag_latitude^2 + lag_longitude^2)

      if (is.na(hmax) || hnorm <= hmax) {
        for (tau in tau_vect) {
          lags$s1[idx] <- i
          lags$s2[idx] <- j
          lags$h1[idx] <- lag_latitude
          lags$h2[idx] <- lag_longitude
          lags$tau[idx] <- tau
          lags$hnorm[idx] <- hnorm
          idx <- idx + 1
        }
      }
    }
  }

  # Remove the unused preallocated space
  lags <- lags[1:(idx - 1), ]

  return(lags)
}



#' Calculate the Euclidean distance between two points.
#'
#' @param point1 The first point.
#' @param point2 The second point.
#' @return The Euclidean distance between the two points.
#'
#' @examples
#' euclidean_distance(c(0, 0), c(3, 4))
#'
#' @export
get_euclidean_distance <- function(point1, point2) {
  sqrt(sum((point1 - point2)^2))
}

#' distances_regular_grid calculates the distances between sites on a regular
#' grid.
#'
#' @param nsites the number of sites on the grid ie nb of pixels squared
#' @return a matrix of distances between sites
#'
#' @import spam
#' @import tidyr
#' @import reshape2
#'
#' @examples
#' distances_regular_grid(25)
#'
#' @export
distances_regular_grid <- function(nsites, adv = c(0, 0), tau = 1:10) {
  # distances
  grid_size <- sqrt(nsites)
  # Calculate the spatial coordinates of the regular grid
  x_coords <- rep(1:grid_size, each = grid_size)
  y_coords <- rep(1:grid_size, times = grid_size)
  grid_points <- cbind(x_coords, y_coords)
  sites_coords <- data.frame(grid_points)
  colnames(sites_coords) <- c("Latitude", "Longitude")
  if (all(adv == c(0, 0))) {
    distances <- matrix(0, nrow = grid_size^2, ncol = grid_size^2)
    for (i in 1:(grid_size^2)) {
      for (j in 1:(grid_size^2)) {
        distances[i, j] <- get_euclidean_distance(grid_points[i, ],
                                                  grid_points[j, ])
      }
    }
  } else {
    distances <- array(0, dim = c(grid_size^2, grid_size^2, length(tau)))
    for (i in 1:(grid_size^2)) {
      for (j in 1:(grid_size^2)) {
        for (k in seq_along(tau)) {
          x1 <- grid_points[i, 1]  - adv[1] * tau[k]
          y1 <- grid_points[i, 2] - adv[2] * tau[k]
          x2 <- grid_points[j, 1] - adv[1] * tau[k]
          y2 <- grid_points[j, 2] - adv[2] * tau[k]
          distances[i, j, k] <- get_euclidean_distance(c(x1, y1), c(x2, y2))
        }
      }
    }
  }
  # Convert array to long dataframe
  df_dist_long <- melt(distances)
  if (all(adv == c(0, 0))) {
    colnames(df_dist_long) <- c("Y", "X", "value")
  } else {
    colnames(df_dist_long) <- c("Y", "X", "tau", "value")
  }
  # df_dist <- reshape_distances(dist_matrix)
  return(df_dist_long)
}


#' generate_grid_coords generates the spatial coordinates of a regular grid.
#'
#' @param grid_size the size of the grid (number of rows/columns)
#' @return a data frame containing the grid coordinates
#'
#' @examples
#' generate_grid_coords(5)
#'
#' @export
generate_grid_coords <- function(grid_size) {
  # Generate x and y coordinates for the grid points
  sites_coords <- data.frame(
    Latitude = rep(1:grid_size, each = grid_size),   # Y
    Longitude = rep(1:grid_size, times = grid_size)  # X
  )
  return(sites_coords)
}
