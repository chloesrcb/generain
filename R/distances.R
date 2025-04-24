#' This function reshapes a distance matrix into a long dataframe.
#'
#' @param locations A matrix or data frame containing the coordinates of the
#' locations.
#' @param dmax The maximum distance threshold if specified.
#' @param latlon Logical indicating whether the coordinates are in latitude and
#' longitude format. Default is TRUE.
#'
#' @return A distance matrix containing the reshaped distances.
#'
#' @import geodist
#'
#' @export
get_dist_mat <- function(locations, dmax = NA, latlon = TRUE) {
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

  return(dist_mat)
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
#' @param tau_vect A vector of temporal lags.
#' @param latlon Logical indicating whether the coordinates are in latitude and
#'               longitude format. Default is FALSE.
#'
#' @import utils
#'
#' @return A dataframe containing the lag vectors.
#'
#' @export
get_lag_vectors <- function(df_coords, tau_vect, latlon = FALSE) {
  # Dimensions
  n <- nrow(df_coords)
  tau_len <- length(tau_vect)

  # Create index combinations
  pairs_comb <- utils::combn(n, 2)

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
                     tau = integer(num_st_lags))
  lags$s1 <- rep(i_vals, each = tau_len)
  lags$s2 <- rep(j_vals, each = tau_len)
  lags$tau <- rep(tau_vect, times = num_pairs)

  # Get coordinates
  lags$s1x <- df_coords$Longitude[lags$s1]
  lags$s1y <- df_coords$Latitude[lags$s1]
  lags$s2x <- df_coords$Longitude[lags$s2]
  lags$s2y <- df_coords$Latitude[lags$s2]

  # Convert to meters using Haversine distance for hx, hy
  if (latlon) {
    lags$hnorm <- haversine_distance_with_advection(lags$s1y, lags$s1x,
                                              lags$s2y, lags$s2x,
                                              c(0, 0), lags$tau)$distance

  } else {
    lags$hx <- lags$s2x - lags$s1x # s - s0
    lags$hy <- lags$s2y - lags$s1y
    lags$hnorm <- sqrt(lags$hx^2 + lags$hy^2) # Euclidean distance
  }

  return(lags)
}

#' get_conditional_lag_vectors function
#'
#' This function calculates the lag vectors between pairs of points.
#'
#' @param df_coords A dataframe containing the coordinates of the points.
#' @param s0 Conditional spatial point coordinates.
#' @param t0 Conditional temporal point.
#' @param tau_vect A vector of temporal lags. Default is 0:10.
#' @param latlon Logical indicating whether the coordinates are in latitude and
#'
#' @import utils
#' @import geosphere
#'
#' @return A dataframe containing the lag vectors.
#'
#' @export
get_conditional_lag_vectors <- function(df_coords, s0 = c(1, 1),
                                    t0 = 0, tau_vect = 0:10, latlon = FALSE) {
  # Conditional point index in df_coords
  if (typeof(s0) == "list") { # TODO t0 inutile
    ind_s0 <- which(df_coords$Latitude == s0$Latitude &
                    df_coords$Longitude == s0$Longitude)
  } else {
    ind_s0 <- which(df_coords$Latitude == s0[2] & df_coords$Longitude == s0[1])
  }

  if (length(ind_s0) == 0) {
    stop("The conditional site (s0) is not found in df_coords")
  }

  # Dimensions
  n <- nrow(df_coords)
  tau_lag <- tau_vect
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

  # s1_coords <- unlist(df_coords[ind_s0, c("Longitude", "Latitude")])
  # s2_coords <- df_coords[unique(lags$s2), c("Longitude", "Latitude")]

  # Vector coordinates between two sites
  # Convert to meters using Haversine distance for hx, hy
  if (latlon) {
    # Geodesic distance
    for (i in 1:nrow(lags)) {
      lags$hnorm <- haversine_distance_with_advection(lags$s1y, lags$s1x,
                                              lags$s2y, lags$s2x,
                                              c(0, 0), lags$tau)$distance

    }
  } else {
    # lags$hx <- s2_coords$Longitude - s1_coords[1]
    # lags$hy <- s2_coords$Latitude - s1_coords[2]
    lags$hx <- lags$s2x  - lags$s1x
    lags$hy <- lags$s2y - lags$s1y
    lags$hnorm <- sqrt(lags$hx^2 + lags$hy^2)
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
    Longitude = rep(1:grid_size, times = grid_size), # X
    Latitude = rep(1:grid_size, each = grid_size)  # Y
  )
  rownames(sites_coords) <- paste0("S", 1:(grid_size^2))
  return(sites_coords)
}


#' generate_realistic_latlon_grid generates a data frame with site IDs and
#' random latitude and longitude values within given ranges.
#'
#' @param n_sites the number of sites to generate
#' @param lat_range the range of latitude values. Default is c(40, 50).
#' @param lon_range the range of longitude values. Default is c(-5, 5).
#' @return a data frame with site IDs and coordinates
#'
#' @export
generate_realistic_latlon_grid <- function(n_sites, lat_range = c(40, 50),
                                                        lon_range = c(-5, 5)) {
  # Generate random lat/lon values within given ranges
  latitudes <- runif(n_sites, min = lat_range[1], max = lat_range[2])
  longitudes <- runif(n_sites, min = lon_range[1], max = lon_range[2])

  # Create a data frame with site IDs
  sites_coords <- data.frame(
    site_id = 1:n_sites,
    Latitude = latitudes,
    Longitude = longitudes
  )

  return(sites_coords)
}


#' haversine_distance_with_advection function
#'
#' This function calculates the Haversine distance between two points
#' with advection applied and gives the direction of the vector in radians
#' with meteorological orientation (0 degrees is North).
#'
#' @param lat1 Latitude of the first point.
#' @param lon1 Longitude of the first point.
#' @param lat2 Latitude of the second point.
#' @param lon2 Longitude of the second point.
#' @param adv A vector of advection values (in m/s).
#' @param tau Temporal lag (in seconds).
#' @return The Haversine distance between the two points with advection applied
#' and the direction of the vector in radians.
#'
#' @export
haversine_distance_with_advection <- function(lat1, lon1, lat2, lon2, adv, tau) {
  # Convert degrees to radians
  to_rad <- function(deg) deg * pi / 180

  # Earth's radius in meters
  R <- 6371000

  # Ensure inputs are numeric vectors of the same length
  if (!all(length(lat1) == length(lat2), length(lat1) == length(tau))) {
    stop("lat/lon and tau inputs must be of the same length.")
  }

  # Initial haversine calculation
  deltaLat <- to_rad(lat2 - lat1)
  deltaLon <- to_rad(lon2 - lon1)

  a <- sin(deltaLat / 2)^2 +
       cos(to_rad(lat1)) * cos(to_rad(lat2)) * sin(deltaLon / 2)^2
  c <- 2 * atan2(sqrt(a), sqrt(1 - a))
  # distance <- R * c

  # Add advection to the coordinates
  adv_x <- adv[1] * tau # m/s * s = m or km/h * h = km
  adv_y <- adv[2] * tau

  # Convert displacement to degrees
  lat2_adj <- lat2 + (adv_y / 111132.92)   # 1° latitude ≈ 111132 m
  lon2_adj <- lon2 + (adv_x / (111412.84 * cos(to_rad(lat2))))  # lon varies by latitude

  # Haversine distance after advection
  deltaLat_adj <- to_rad(lat2_adj - lat1)
  deltaLon_adj <- to_rad(lon2_adj - lon1)

  a_adj <- sin(deltaLat_adj / 2)^2 +
           cos(to_rad(lat1)) * cos(to_rad(lat2_adj)) * sin(deltaLon_adj / 2)^2
  c_adj <- 2 * atan2(sqrt(a_adj), sqrt(1 - a_adj))
  distance_adj <- R * c_adj

  # Direction calculation
  theta <- atan2(
    sin(deltaLon_adj) * cos(to_rad(lat2_adj)),
    cos(to_rad(lat1)) * sin(to_rad(lat2_adj)) -
      sin(to_rad(lat1)) * cos(to_rad(lat2_adj)) * cos(deltaLon_adj)
  )

  # Meteorological convention (from North, clockwise)
  theta_meteo <- (pi / 2 - theta) %% (2 * pi)

  # In kilometers
  distance_km <- distance_adj / 1000
  distance_km[distance_km < 1e-6] <- 0

  return(list(distance = distance_km, theta = theta, theta_meteo = theta_meteo))
}


# haversine_distance_with_advection <- function(lat1, lon1, lat2, lon2, adv,
#                                               tau) {
 
#   deltaLat <- (lat2 - lat1) * pi / 180
#   deltaLon <- (lon2 - lon1) * pi / 180
#   theta <- atan2(deltaLat, deltaLon)  # in radians
#   a <- sin(deltaLat / 2)^2 + cos(lat1 * pi / 180) *
#                             cos(lat2 * pi / 180) * sin(deltaLon / 2)^2
#   c <- 2 * atan2(sqrt(a), sqrt(1 - a))
#   distance <- 6371000 * c  # Radius of the Earth in meters

#   # With advection
#   if (length(adv) == 2) {
#     adv_x <- adv[1] * tau  * 3600
#     adv_y <- adv[2] * tau  * 3600

#     # Adjust the coordinates based on advection
#     lon2_adj <- lon2 + adv_x / (111319 * cos(lat2 * pi / 180))
#     lat2_adj <- lat2 + adv_y / 111319

#     # Recalculate the distance after adjustment
#     deltaLat_adj <- (lat2_adj - lat1) * pi / 180
#     deltaLon_adj <- (lon2_adj - lon1) * pi / 180
#     a_adj <- sin(deltaLat_adj / 2)^2 + cos(lat1 * pi / 180) *
#                         cos(lat2_adj * pi / 180) * sin(deltaLon_adj / 2)^2
#     c_adj <- 2 * atan2(sqrt(a_adj), sqrt(1 - a_adj))
#     distance <- 6371000 * c_adj  # Radius of the Earth in meters
#     theta <- atan2(deltaLat_adj, deltaLon_adj)  # in radians
#   }

#   theta_meteo <- (pi / 2 - theta) %% (2 * pi) # Meteorological direction
#   # If the distance is less than a very small value, set to zero for
#   # each element
#   distance <- ifelse(distance < 1e-06, 0, distance)

#   # Return the norm and angle (theta)
#   return(list(distance = distance / 1000, theta = theta,
#                                           theta_meteo = theta_meteo))
# }

