
#' Compute the distance matrix between locations
#'
#' This function computes the distance matrix between a set of locations.
#'
#' @param locations A matrix or data frame containing the coordinates of the
#' locations.
#' @param dmax The maximum distance threshold if specified.
#' @param latlon Logical indicating whether the coordinates are in latitude and
#' longitude format.
#' @return A distance matrix where each element represents the distance between
#' two locations.
#'
#' @import geodist
#' @importFrom stats dist
#'
#' @export
get_dist_mat <- function(locations, dmax = NA, latlon = TRUE) {
  # get longitude and latitude in a dataframe to get distance between points
  loc <- data.frame(lat = locations$Latitude, lon = locations$Longitude)
  # get distance matrix
  if(latlon) {
    dist_mat <- geodist(loc, measure = "haversine") # meters by default
  } else {
    dist_mat <- as.matrix(dist(loc))
  }
  if (!is.na(dmax)) {
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
#' @import spam
#' @import terra
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
  df_dist <- na.omit(melt(df_dist, "Y", variable_name = "X"))
  df_dist$X <- factor(df_dist$X, levels = rev(levels(df_dist$X)))

  return(df_dist)
}


#' get_h_vect function
#'
#' This function calculates the spatial lag vector based on the given distance 
#' dataframe and maximum spacial lag value.
#'
#' @param df_dist The distance dataframe.
#' @param hmax The maximum spacial lag value.
#' @return The h vector.
#'
#' @import terra
#'
#' @examples
#' hmax <- 2
#' get_h_vect(df_dist, hmax)
#'
#' @export
get_h_vect <- function(df_dist, hmax) {
  # get unique distances
  h_vect <- sort(df_dist$value)
  h_vect <- unique(h_vect)
  if (!is.na(hmax)) {
      h_vect <- h_vect[h_vect <= hmax]
  }
  return(h_vect)
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
#' @import geodist
#' @import stats
#' @import terra
#' @import reshape2
#' @import spam
#' @import tidyr
#'
#' @examples
#' distances_regular_grid(10)
#'
#' @export
distances_regular_grid <- function(nsites) {
  # distances
  grid_size <- sqrt(nsites)
  # Calculate the spatial coordinates of the regular grid
  x_coords <- rep(1:grid_size, each = grid_size)
  y_coords <- rep(1:grid_size, times = grid_size)
  grid_points <- cbind(x_coords, y_coords)
  sites_coords <- data.frame(grid_points)
  colnames(sites_coords) <- c("Latitude", "Longitude")
  distances <- matrix(0, nrow = grid_size^2, ncol = grid_size^2)

  for (i in 1:(grid_size^2)) {
    for (j in 1:(grid_size^2)) {
      distances[i, j] <- get_euclidean_distance(grid_points[i, ],
                                                grid_points[j, ])
    }
  }

  dist_matrix <- as.matrix(distances)
  rownames(dist_matrix) <- c(1:nsites)
  colnames(dist_matrix) <- c(1:nsites)
  df_dist <- reshape_distances(dist_matrix)
  return(df_dist)
}