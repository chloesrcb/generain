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
#' @param dist The distance dataframe or matrix for intervals.
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

#' Calculate the Euclidean distance between two points.
#'
#' @param point1 The first point.
#' @param point2 The second point.
#' @return The Euclidean distance between the two points.
#'
#' @examples
#' get_euclidean_distance(c(0, 0), c(3, 4))
#'
#' @export
get_euclidean_distance <- function(point1, point2) {
  sqrt(sum((point1 - point2)^2))
}

#' get_lag_vectors function
#'
#' This function calculates the lag vectors between pairs of points.
#'
#' @param df_coords A dataframe containing the coordinates of the points.
#' @param params A vector of parameters.
#' @param tau_vect A vector of temporal lags. Default is 0:10.
#' @param hmax The maximum distance threshold. Default is NA.
#' @param norm The norm to use. Can be "euclidean" or "Lp". Default
#'             is "euclidean".
#'
#' @import utils
#'
#' @return A dataframe containing the lag vectors.
#'
#' @export
get_lag_vectors <- function(df_coords, params, tau_vect = 0:10, hmax = NA,
                            norm = "euclidean") {
  # Advection vector
  adv <- if (length(params) == 6) params[5:6] else c(0, 0)

  # Dimensions
  n <- nrow(df_coords)
  tau_len <- length(tau_vect)

  # Create index combinations
  pairs_comb <- utils::combn(n, 2)
  # indices <- pairs_comb

  # Add n,n pairs
  diag_pairs <- matrix(rep(1:n, each = 2), nrow = 2)
  indices <- cbind(pairs_comb, diag_pairs)

  # Get indices
  i_vals <- indices[1, ]
  j_vals <- indices[2, ]
  # Number of pairs
  num_pairs <- length(i_vals)
  # Number of spatio-temporal lags
  num_st_lags <- num_pairs * tau_len

  # init dataframe
  lags <- data.frame(s1 = integer(num_st_lags), s2 = integer(num_st_lags),
                     s1x = integer(num_st_lags), s1y = integer(num_st_lags),
                     s2x = integer(num_st_lags), s2y = integer(num_st_lags),
                     hx = numeric(num_st_lags), hy = numeric(num_st_lags),
                     tau = integer(num_st_lags), hnorm = numeric(num_st_lags))
  lags$s1 <- rep(i_vals, each = tau_len)
  lags$s2 <- rep(j_vals, each = tau_len)
  lags$tau <- rep(tau_vect, times = num_pairs)

  # Get coordinates
  lags$s1x <- df_coords$Longitude[lags$s1]
  lags$s1y <- df_coords$Latitude[lags$s1]
  lags$s2x <- df_coords$Longitude[lags$s2]
  lags$s2y <- df_coords$Latitude[lags$s2]

  # Vector coordinates between two sites
  lags$hx <- lags$s1x - lags$s2x
  lags$hy <- lags$s1y - lags$s2y

  # Add advection
  lags$hx <- lags$hx - adv[1] * lags$tau
  lags$hy <- lags$hy - adv[2] * lags$tau

  # Calculate the norm
  if (norm == "euclidean") {
    lags$hnorm <- sqrt(lags$hx^2 +  lags$hy^2)
  } else if (norm == "Lp") {
    lags$hnorm <- norm_Lp(lags$hx, lags$hy, params[3])
  }

  # Filter based on hmax
  if (!is.na(hmax)) {
    lags <- lags[lags$hnorm <= hmax, ]
  }

  return(lags)
}


#' get_conditional_lag_vectors function
#'
#' This function calculates the lag vectors between pairs of points.
#'
#' @param df_coords A dataframe containing the coordinates of the points.
#' @param params A vector of parameters.
#' @param s0 Conditional spatial point coordinates.
#' @param t0 Conditional temporal point.
#' @param tau_vect A vector of temporal lags. Default is 0:10.
#' @param hmax The maximum distance threshold. Default is NA.
#' @param norm The norm to use. Can be "euclidean" or "Lp". Default
#'            is "euclidean".
#'
#' @import utils
#'
#' @return A dataframe containing the lag vectors.
#'
#' @export
get_conditional_lag_vectors <- function(df_coords, params, s0 = c(1, 1),
                                        t0 = 1, tau_vect = 0:10, hmax = NA,
                                        norm = "euclidean") {
  # Advection
  adv <- if (length(params) == 6) params[5:6] else c(0, 0)

  # Conditional point index in df_coords
  ind_s0 <- which(df_coords$Latitude == s0[1] & df_coords$Longitude == s0[2])

  if (length(ind_s0) == 0) {
    stop("The conditional site (s0) is not found in df_coords")
  }

  # Dimensions
  n <- nrow(df_coords)
  tau_lag <- sort(unique(abs(tau_vect - t0))) # tau - t0
  tau_len <- length(tau_lag)

  # Index of pairs (s0, si)
  j_vals <- 1:n

  # Number of pairs
  num_pairs <- n
  # Number of spatio-temporal lags
  num_st_lags <- n * tau_len

  # init dataframe
  lags <- data.frame(s1 = integer(num_st_lags), s2 = integer(num_st_lags),
                     s1x = integer(num_st_lags), s1y = integer(num_st_lags),
                     s2x = integer(num_st_lags), s2y = integer(num_st_lags),
                     hx = numeric(num_st_lags), hy = numeric(num_st_lags),
                     tau = integer(num_st_lags), hnorm = numeric(num_st_lags))
  lags$s1 <- rep(ind_s0, times = num_st_lags)
  lags$s2 <- rep(j_vals, each = tau_len)
  lags$tau <- rep(tau_lag, times = num_pairs)

  # Get coordinates
  lags$s1x <- df_coords$Longitude[lags$s1]
  lags$s1y <- df_coords$Latitude[lags$s1]
  lags$s2x <- df_coords$Longitude[lags$s2]
  lags$s2y <- df_coords$Latitude[lags$s2]

  # Vector coordinates between two sites
  lags$hx <- lags$s1x - lags$s2x
  lags$hy <- lags$s1y - lags$s2y

  # Add advection
  lags$hx <- lags$hx - adv[1] * lags$tau
  lags$hy <- lags$hy - adv[2] * lags$tau

  # Calculate the norm
  if (norm == "euclidean") {
    lags$hnorm <- sqrt(lags$hx^2 +  lags$hy^2)
  } else if (norm == "Lp") {
    lags$hnorm <- norm_Lp(lags$hx, lags$hy, params[3])
  }

  # Filter < hmax
  if (!is.na(hmax)) {
    lags <- lags[lags$hnorm <= hmax, ]
  }

  return(lags)
}



#' distances_regular_grid calculates the distances between sites on a regular
#' grid.
#'
#' @param nsites the number of sites on the grid ie nb of pixels squared
#' @param adv a vector of advection values. Default is c(0, 0).
#' @param tau a vector of temporal lags. Default is 1:10.
#'
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
