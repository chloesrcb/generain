#' vario function
#'
#' This function computes the variogram value for a given distance matrix and
#' variogram parameters.
#'
#' @param h The spatial lags.
#' @param tau The temporal lags.
#' @param beta1 The value of beta1.
#' @param beta2 The value of beta2.
#' @param alpha1 The value of alpha1.
#' @param alpha2 The value of alpha2.
#'
#' @return The variogram value.
#'
#' @export
gamma_theta <- function(h, tau, beta1, beta2, alpha1, alpha2) {
  gamma <- (h * beta1)^alpha1 + (abs(tau) * beta2)^alpha2
  return(gamma)
}


#' dist_adv function
#'
#' This function computes the distance between two points with advection.
#'
#' @param s1 The first point coordinates.
#' @param s2 The second point coordinates.
#' @param t1 The first time.
#' @param t2 The second time.
#' @param adv The advection coordinates vector.
#'
#' @return The distance between the two points with advection.
#'
#' @export
dist_adv <- function(s1, s2, t1, t2, adv) {
  h <- s1 - s2
  tau <- abs(t1 - t2)
  hnorm_adv <- sqrt((h[1] - adv[1] * tau)^2 + (h[2] - adv[2] * tau)^2)
  return(hnorm_adv)
}

#' conditional_variogram function
#'
#' This function computes the spatio-temporal conditional variogram.
#'
#' @param x Vector for the first dimension (spatial x in our case).
#' @param y Vector for the second dimension (spatial y in our case)
#' @param t Vector for the third dimension (time in our case).
#' @param s0 Vector of dimension 2 for the spatial conditioning point.
#' @param t0 Value of the temporal conditioning point.
#' @param grid The grid matrix.
#' @param model The model to use for the variogram.
#' @param adv The advection coordinates vector. Default is c(0, 0).
#'
#' @return The spatio-temporal conditional variogram.
#'
#' @export
conditional_variogram <- function(x, y, t, s0, t0, grid, model, adv = c(0, 0)) {
  lx <- length(x)
  ly <- length(y)
  lt <- length(t)

  # Spatial conditioning point
  s0_x <- s0[2]
  s0_y <- s0[1]

  gamma_s0_t0 <- array(0, dim = c(lx, ly, lt))
  for (i in seq_len(lx)) {
      for (j in seq_len(ly)) {
          for (k in seq_len(lt)) {
            gamma_s0_t0[i, j, k] <- RandomFields::RFvariogram( # i,j or j,i??
                model,
                x = x[i] - s0_x - adv[1] * (t[k] - t0),
                y = y[j] - s0_y - adv[2] * (t[k] - t0),
                z = t[k] - t0
            )
          }
      }
  }
  return(gamma_s0_t0)
}


#' advected_variogram function
#'
#' This function computes the spatio-temporal advected variogram.
#'
#' @param x Vector for the first dimension (spatial x in our case).
#' @param y Vector for the second dimension (spatial y in our case)
#' @param t Vector for the third dimension (time in our case).
#' @param grid The grid matrix.
#' @param model The model to use for the variogram.
#' @param adv The advection coordinates vector. Default is c(0, 0).
#'
#' @return An array for the spatio-temporal advected variogram.
#'
#'
#' @export
advected_variogram <- function(x, y, z, grid, model, adv) {
  lx <- length(x)
  ly <- length(y)
  lz <- length(z)
  N <- length(grid[, 1]) # spatio-temporal dimension

  gamma <- array(0, dim = c(lx, ly, lz, N)) # variogram

  for (n in seq_len(N)) {
      for (i in seq_len(lx)) {
          for (j in seq_len(ly)) {
              for (k in seq_len(lz)) {
                  gamma[j, i, k, n] <- RandomFields::RFvariogram(
                      model,
                      x = x[i] - grid[n, 1] - adv[1] * (z[k] - grid[n, 3]),
                      y = y[j] - grid[n, 2] - adv[2] * (z[k] - grid[n, 3]),
                      z = z[k] - grid[n, 3]
                  )
              }
          }
      }
  }
  return(gamma)
}

#' compute_gamma_point function
#'
#' This function computes the spatio-temporal variogram for a given point.
#'
#' @param grid The grid matrix.
#' @param gamma_space The spatial variogram.
#' @param gamma_temp The temporal variogram.
#' @param adv The advection coordinates vector.
#' @param s The spatial coordinates vector. Default is NA.
#' @param t The temporal time. Default is NA.
#'
#' @return The spatio-temporal variogram for the given point.
#'
#' @export
compute_gamma_point <- function(grid, gamma_space, gamma_temp, adv, s=NA, t=NA) {
  # Get length of unique values
  lx <- length(unique(grid$x))
  ly <- length(unique(grid$y))
  lt <- length(unique(grid$t))

  # Get spatial variogram for s
  if(all(!(is.na(s)) & !(is.na(t)))) {
    index_grid <- which(grid$shifted_x == s[1] & grid$shifted_y == s[2] &
                          grid$t == t)
    gamma_s <- gamma_space[[index_grid]]
    gamma_t <- gamma_temp[[t + 1]]
  } else {
    gamma_s <- gamma_space
    gamma_t <- gamma_temp
  }

  coords <- cbind(grid$shifted_x, grid$shifted_y)
  if (all(adv == 0)) {
    duplicates <- duplicated(coords)
    filtered_coords <- coords[!duplicates, ]
    coords <- filtered_coords
  }
  nsites <- nrow(coords)

  # t_index <- grid$t + 1 # index starts at 1
  t_levels <- sort(unique(grid$t))
  t_index  <- match(grid$t, t_levels)  # 1..lt


  # compute gamma space-time
  gamma <- array(NA, dim = c(lx, ly, lt))
  for (i in seq_along(t_index)) {
    if (all(adv == 0) & i > nsites) {
        if (i %% nsites == 0) { # get index when no advection
          ind_g_s <- i - nsites * (i %/% nsites - 1)
        } else {
          ind_g_s <- i - nsites * (i %/% nsites)
        }
    } else {
      ind_g_s <- i # get index when advection
    }
    vario <- gamma_s[ind_g_s] + gamma_t[t_index[i]]
    gamma[grid$x[i], grid$y[i], t_index[i]] <- vario
  }

  return(gamma)
}


#' compute_st_variogram function
#'
#' Computes the spatio-temporal variogram over a grid.
#'
#' @param grid The grid data.frame with columns x, y, t, shifted_x, shifted_y.
#' @param gamma_space The isotropic spatial variogram (optional).
#' @param gamma_space_x The spatial variogram in the x direction (optional).
#' @param gamma_space_y The spatial variogram in the y direction (optional).
#' @param gamma_temp The temporal variogram.
#' @param adv The advection vector (length 2 numeric).
#'
#' @return A 3D array with the spatio-temporal variogram values over the grid.
#'
#' @export
compute_st_variogram <- function(grid,
                                 gamma_space = NULL,
                                 gamma_space_x = NULL,
                                 gamma_space_y = NULL,
                                 gamma_temp,
                                 adv) {
  lx <- length(unique(grid$x))
  ly <- length(unique(grid$y))
  lt <- length(unique(grid$t))

  # Determine distance type
  is_lalpha <- !is.null(gamma_space_x) && !is.null(gamma_space_y)
  is_euclidean <- !is.null(gamma_space)

  if (!is_euclidean && !is_lalpha) {
    stop("Provide either gamma_space for euclidean or gamma_space_x and gamma_space_y for lalpha case.")
  }

  # Determine gamma components
  if (!any(is.na(adv)) && all(adv == 0)) {
    coords <- cbind(grid$shifted_x, grid$shifted_y)
    coords <- coords[!duplicated(coords), ]
  } else {
    coords <- cbind(grid$shifted_x, grid$shifted_y)
  }
  if (is.vector(coords)) {
    coords <- matrix(coords, nrow = 1)
  }
  nsites <- nrow(coords)
  # t_index <- grid$t + 1
  t_levels <- sort(unique(grid$t))
  t_index  <- match(grid$t, t_levels)  # 1..lt

  gamma <- array(NA, dim = c(lx, ly, lt))

  for (i in seq_along(t_index)) {
    if (!any(is.na(adv)) && all(adv == 0) && i > nsites) {
      ind_g_s <- if (i %% nsites == 0) {
        i - nsites * (i %/% nsites - 1)
      } else {
        i - nsites * (i %/% nsites)
      }
    } else {
      ind_g_s <- i
    }

    if (is_lalpha) {
      vario <- gamma_space_x[[ind_g_s]] + gamma_space_y[[ind_g_s]] +
                                                      gamma_temp[[t_index[i]]]
    } else {
      vario <- gamma_space[[ind_g_s]] + gamma_temp[[t_index[i]]]
    }

    gamma[grid$x[i], grid$y[i], t_index[i]] <- vario
  }

  return(gamma)
}


#' compute_W_s_t function
#'
#' This function computes the spatio-temporal gaussian random field given the
#' spatial and temporal gaussian random fields and the advection.
#'
#' @param grid The grid matrix.
#' @param W_s The spatial gaussian random field.
#' @param W_t The temporal gaussian random field.
#' @param adv The advection coordinates vector.
#'
#' @return The spatio-temporal gaussian random field.
#'
#' @export
compute_W_s_t <- function(grid, W_s, W_t, adv) {
  lx <- length(unique(grid$x))
  ly <- length(unique(grid$y))
  lt <- length(unique(grid$t))
  coords <- cbind(grid$shifted_x, grid$shifted_y)
  if (all(adv == 0)) {
    duplicates <- duplicated(coords)
    filtered_coords <- coords[!duplicates, ]
    coords <- filtered_coords
  }
  nsites <- nrow(coords)
  W_s_t <- array(NA, dim = c(lx, ly, lt))
  # t_index <- grid$t + 1 # index starts at 1
  t_levels <- sort(unique(grid$t))
  t_index  <- match(grid$t, t_levels)  # 1..lt

  for (i in seq_len(nrow(grid))) {
      s_x <- grid$x[i]
      s_y <- grid$y[i]
      s <- c(s_x, s_y)
      if (all(adv == 0) & i > nsites) {
        if (i %% nsites == 0) {
          ind_W_s <- i - nsites * (i %/% nsites - 1)
        } else {
          ind_W_s <- i - nsites * (i %/% nsites)
        }
      } else {
        ind_W_s <- i
      }
      W_s_t[s[1], s[2], t_index[i]] <- W_s[ind_W_s] + W_t[t_index[i]]
  }
  return(W_s_t)
}


#' compute_st_gaussian_process function
#'
#' This function computes the spatio-temporal Gaussian random field given the
#' spatial and temporal Gaussian random fields and the advection.
#'
#' @param grid The grid matrix.
#' @param W_s The spatial Gaussian random field (for euclidean case).
#' @param W_s_x The spatial Gaussian random field in the x direction (for lalpha case).
#' @param W_s_y The spatial Gaussian random field in the y direction (for lalpha case).
#' @param W_t The temporal Gaussian random field.
#' @param adv The advection coordinates vector.
#'
#' @return The spatio-temporal Gaussian random field.
#'
#' @export
compute_st_gaussian_process <- function(grid, W_s = NULL,
                                        W_s_x = NULL, W_s_y = NULL, W_t,
                                        adv) {
  # Dimensions of the original (unshifted) grid
  lx <- length(unique(grid$x))
  ly <- length(unique(grid$y))
  lt <- length(unique(grid$t))

  # Shifted coordinates (advection already applied upstream)
  coords <- cbind(grid$shifted_x, grid$shifted_y)

  # Determine distance type
  is_lalpha <- !is.null(W_s_x) && !is.null(W_s_y)
  is_euclidean <- !is.null(W_s)

  if (!is_euclidean && !is_lalpha) {
    stop("Provide either W_s (euclidean) or W_s_x and W_s_y (lalpha).")
  }

  key <- paste(coords[,1], coords[,2], sep = "_")
  key_u <- unique(key)
  map_idx <- match(key, key_u)

  n_unique <- length(key_u)

  if (is_euclidean) {
    W_s_vec <- as.numeric(W_s)
    if (length(W_s_vec) != n_unique && length(W_s_vec) != nrow(grid)) {
      stop("Length of W_s does not match n_unique(coords) nor nrow(grid).")
    }
  } else {
    W_sx_vec <- as.numeric(W_s_x)
    W_sy_vec <- as.numeric(W_s_y)
    if ((length(W_sx_vec) != n_unique && length(W_sx_vec) != nrow(grid)) ||
        (length(W_sy_vec) != n_unique && length(W_sy_vec) != nrow(grid))) {
      stop("Length of W_s_x / W_s_y does not match n_unique(coords) nor nrow(grid).")
    }
  }

  # If grid$t are not 0:(lt-1), match them to positions in sorted unique times
  t_levels <- sort(unique(grid$t))
  t_index <- match(grid$t, t_levels)  # 1..lt

  W_t_vec <- as.numeric(W_t)
  if (length(W_t_vec) != lt) {
    stop("Length of W_t must equal the number of unique time points in grid$t.")
  }

  W_s_t <- array(NA_real_, dim = c(lx, ly, lt))

  for (i in seq_len(nrow(grid))) {
    sx <- grid$x[i]
    sy <- grid$y[i]
    ti <- t_index[i]

    if (is_euclidean) {
      ind_ws <- if (length(W_s_vec) == nrow(grid)) i else map_idx[i]
      W_s_t[sx, sy, ti] <- W_s_vec[ind_ws] + W_t_vec[ti]
    } else {
      ind_ws <- if (length(W_sx_vec) == nrow(grid)) i else map_idx[i]
      W_s_t[sx, sy, ti] <- W_sx_vec[ind_ws] + W_sy_vec[ind_ws] + W_t_vec[ti]
    }
  }

  return(W_s_t)
}


#' sim_BR function
#'
#' This function performs a simulation of a spatio-temporal Brown-Resnick
#' process using a fractionnal Brownian motion model and based on the
#' David Leber code with advection.
#'
#' @param beta1 The value of beta1.
#' @param beta2 The value of beta2.
#' @param alpha1 The value of alpha1.
#' @param alpha2 The value of alpha2.
#' @param x Vector for the first dimension (spatial x in our case).
#' @param y Vector for the second dimension (spatial y in our case)
#' @param t Vector for the third dimension (time in our case).
#' @param adv The advection coordinates vector. Default is c(0, 0).
#' @param nres The number of simulations to perform. Default is 1.
#'
#' @return The result of the simulation.
#'
#' @import stats
#' @import RandomFields
#' @import RandomFieldsUtils
#'
#' @export
sim_BR <- function(beta1, beta2, alpha1, alpha2, x, y, t, adv = c(0, 0), 
                   nres = 1) {
  ## Setup
  options(str = list())
  # RandomFields::RFoptions(spConform = FALSE)
  RandomFields::RFoptions(spConform = FALSE, allow_duplicated_locations = T)

  lx <- length(sx <- seq_along(x))
  ly <- length(sy <- seq_along(y))
  lt <- length(st <- seq_along(t))

  ## Model-Variogram BuhlCklu
  modelSpace <- RandomFields::RMfbm(alpha = alpha1, var = 2*beta1)
  modelTime <- RandomFields::RMfbm(alpha = alpha2, var = 2*beta2)

  ## Construct grid
  grid_with_advection <- expand.grid(
    x = seq_len(lx),
    y = seq_len(ly),
    t = t
  )

  grid_with_advection$shifted_x <- grid_with_advection$x +
                                    grid_with_advection$t * adv[1]
  grid_with_advection$shifted_y <- grid_with_advection$y +
                                    grid_with_advection$t * adv[2]

  grid <- grid_with_advection

  coords <- grid[, 4:5]

  if (all(adv == 0)) {
    duplicates <- duplicated(coords)
    filtered_coords <- coords[!duplicates, ]
    coords <- filtered_coords
  }

  # Variogram
  N <- nrow(grid) # number of points
  # x_shifted <- rep(x, times = ly * lt)
  # y_shifted <- rep(rep(y, each = lx), times = lt)
  gamma_space <- lapply(seq_len(N), function(n) 
          RandomFields::RFvariogram(modelSpace,
            x = x_shifted - grid$shifted_x[n],
            y = y_shifted - grid$shifted_y[n])) # for all s,t

  gamma_temp <- lapply(seq_len(lt), function(n) # for all t
          RandomFields::RFvariogram(modelTime,
            x = t - t[n]))

  ## Main
  Z <- array(, dim = c(lx, ly, lt, nres)) # 4d array
  E <- matrix(rexp(nres * N), nrow = nres, ncol = N)
  for (i in seq_len(nres)) { ## n=1
    V <- 1 / E[i, 1]
 
    # Spatial gaussian random field on shifted coords
    W_s <- RandomFields::RFsimulate(modelSpace, coords[, 1], coords[, 2],
                                      grid = FALSE)
    # Temporal gaussian random field
    W_t <- RandomFields::RFsimulate(modelTime, t, n = 1, grid = TRUE)
    # Spatio-temporal random field
    W <- compute_W_s_t(grid, W_s, W_t, adv)

    Y <- exp(W - W[1] - gamma[, , , 1])
    Z[, , , i] <- V * Y
    for (n in 2:N) {
      Exp <- E[i, n]
      V <- 1 / Exp 
      while (V > Z[N * (i - 1) + n]) {
        # Spatial gaussian random field on shifted coords
        W_s <- RandomFields::RFsimulate(modelSpace, coords[, 1], coords[, 2],
                                          grid = FALSE)
        # Temporal gaussian random field
        W_t <- RandomFields::RFsimulate(modelTime, t, n = 1, grid = TRUE)
        # Spatio-temporal gaussian random field
        W <- compute_W_s_t(grid, W_s, W_t, adv)
        s <- c(grid$shifted_x[n], grid$shifted_y[n])
        time <- grid$t[n]
        gamma <- compute_gamma_point(grid, gamma_space, gamma_temp, adv, s,
                                     time)
        Y <- exp(W - W[n] - gamma[, , , n])
        if(all(V * Y[seq_len(n-1)] < Z[(N*(i-1)+1):(N*(i-1)+(n-1))]))
          Z[, , , i] <- pmax(V * Y, Z[, , , i])
          Exp <- Exp + rexp(1)
          V <- 1 / Exp
      }
    }
  }
  Z
}

#' sim_rpareto function (extended)
#'
#' This function simulates a spatio-temporal r-Pareto process using either 
#' a fractional Brownian motion model.
#' It supports conditioning on one or multiple spatial points s0.
#'
#' @param beta1, beta2 Variogram scale parameters for space and time.
#' @param alpha1, alpha2 Variogram smoothness parameters for space and time.
#' @param x, y, t Vectors representing spatial (x, y) and temporal (t) grids.
#' @param adv Advection vector (default = c(0, 0)).
#' @param t0 Time point at which the process is conditioned.
#' @param nres Number of simulations to perform.
#' @param random_s0 If TRUE, selects conditioning point(s) randomly within radius.
#' @param s0 Conditioning spatial location(s). 
#'        Either a vector c(x,y) for one point, or a matrix/data.frame with 2 columns for multiple points.
#' @param s0_radius Radius used if random_s0 is TRUE.
#' @param distance Type of distance metric to use ("euclidean" or "lalpha").
#'
#' @return A list with simulated field Z (5D array: x × y × t × nres × n_s0)
#'         and list of conditioning points s0_used.
#' @export
sim_rpareto <- function(beta1, beta2, alpha1, alpha2, x, y, t,
                        adv = c(0, 0), t0 = 0, nres = 1,
                        random_s0 = FALSE, s0 = c(1, 1),
                        s0_radius = Inf,
                        distance = "euclidean", seed = NULL) {
  if (!is.null(seed)) set.seed(seed)

  RandomFields::RFoptions(spConform = FALSE, allow_duplicated_locations = TRUE,
                          install = "no")

  if (!(distance %in% c("euclidean", "lalpha"))) {
    stop('Invalid distance type. Choose either "euclidean" or "lalpha".')
  }

  lx <- length(x); ly <- length(y); lt <- length(t)
  site_names <- paste0("S", seq_len(lx * ly))

  modelSpace <- RandomFields::RMfbm(alpha = alpha1, var = 2 * beta1)
  modelTime  <- RandomFields::RMfbm(alpha = alpha2, var = 2 * beta2)

  # adv handling
  if (is.null(dim(adv))) {
    # vector length 2
    if (length(adv) != 2) stop("adv must be length-2 vector or nres x 2 matrix/data.frame.")
    adv_mat <- matrix(rep(adv, nres), ncol = 2, byrow = TRUE)
  } else {
    adv_mat <- as.matrix(adv)
    if (ncol(adv_mat) != 2) stop("adv must have 2 columns (adv1, adv2).")
    if (nrow(adv_mat) != nres) stop("If adv is a matrix/data.frame, it must have nres rows.")
  }

  # Normalize s0 input
  if (is.null(dim(s0))) {
    s0_df <- data.frame(x = s0[1], y = s0[2])
  } else {
    s0_df <- as.data.frame(s0); names(s0_df) <- c("x", "y")
  }
  n_s0 <- nrow(s0_df)

  # If random s0, candidates
  if (random_s0) {
    s0_center <- as.numeric(s0_df[1, ])
    grid_points <- expand.grid(x = seq_len(lx), y = seq_len(ly))
    distances <- sqrt((grid_points$x - s0_center[1])^2 +
                      (grid_points$y - s0_center[2])^2)
    candidate_points <- grid_points[distances <= s0_radius, ]
    if (nrow(candidate_points) == 0) stop("No grid points found within specified radius.")
  }

  s0_list <- vector("list", nres)
  Z <- array(NA_real_, dim = c(lx, ly, lt, nres, n_s0))

  # Base grid
  base_grid <- expand.grid(x = seq_len(lx), y = seq_len(ly), t = t)
  base_grid$site <- rep(site_names, times = lt)

  for (i in seq_len(nres)) {
    s0_list[[i]] <- list()

    adv_i <- adv_mat[i, ]  # (adv1, adv2)

    # grid for one episode
    grid <- base_grid
    grid$shifted_x <- grid$x - grid$t * adv_i[1]
    grid$shifted_y <- grid$y - grid$t * adv_i[2]
    coords <- grid[, c("shifted_x", "shifted_y")]

    if (all(adv_i == 0)) {
      coords <- coords[!duplicated(coords), ]
    }

    for (j in seq_len(n_s0)) {
      # Select s0
      s0_curr <- s0_df[j, , drop = FALSE]
      if (random_s0) {
        s0_curr <- candidate_points[sample(nrow(candidate_points), 1), , drop = FALSE]
      }
      s0_list[[i]][[j]] <- s0_curr

      ind_s0_t0 <- which(grid$x == s0_curr$x & grid$y == s0_curr$y & grid$t == t0)
      if (length(ind_s0_t0) != 1) stop("Conditioning index ind_s0_t0 not found uniquely.")
      
      # Temporal variogram at t0
      gamma_temp <- RandomFields::RFvariogram(modelTime, x = t - grid$t[ind_s0_t0])

      # Spatial variogram and Gaussian process
      if (distance == "lalpha") {
        gamma_space_x <- RandomFields::RFvariogram(
          modelSpace, x = coords$shifted_x - grid$shifted_x[ind_s0_t0])
        gamma_space_y <- RandomFields::RFvariogram(
          modelSpace, x = coords$shifted_y - grid$shifted_y[ind_s0_t0])

        gamma_0 <- compute_st_variogram(
          grid,
          gamma_space_x = gamma_space_x,
          gamma_space_y = gamma_space_y,
          gamma_temp = gamma_temp,
          adv = adv_i
        )

        W_s_x <- RandomFields::RFsimulate(modelSpace, coords[, 1], grid = FALSE)
        W_s_y <- RandomFields::RFsimulate(modelSpace, coords[, 2], grid = FALSE)
        W_t   <- RandomFields::RFsimulate(modelTime, t, n = 1, grid = TRUE)

        W <- compute_st_gaussian_process(
          grid, W_s_x = W_s_x, W_s_y = W_s_y, W_t = W_t, adv = adv_i
        )
      } else {
        gamma_space <- RandomFields::RFvariogram(
          modelSpace,
          x = coords$shifted_x - grid$shifted_x[ind_s0_t0],
          y = coords$shifted_y - grid$shifted_y[ind_s0_t0]
        )

        gamma_0 <- compute_st_variogram(
          grid,
          gamma_space = gamma_space,
          gamma_temp = gamma_temp,
          adv = adv_i
        )

        W_s <- RandomFields::RFsimulate(modelSpace, coords[, 1], coords[, 2], grid = FALSE)
        W_t <- RandomFields::RFsimulate(modelTime, t, n = 1, grid = TRUE)

        W <- compute_st_gaussian_process(
          grid, W_s = W_s, W_t = W_t, adv = adv_i
        )
      }

      Y <- exp(W - W[ind_s0_t0] - gamma_0)
      R <- 1 / runif(1) # Pareto(1)
      Z[,,, i, j] <- R * Y
    }
  }

  return(list(Z = Z, s0_used = s0_list))
}


#' sim_rpareto_coords function
#' This function simulates a spatio-temporal r-Pareto process using a
#' fractional Brownian motion model on given coordinates grid, regular or not.
#' @param coords Data frame or matrix with columns 'Longitude' and 'Latitude'.
#' @param times Vector of time points.
#' @param beta1, beta2 Variogram scale parameters for space and time.
#' @param alpha1, alpha2 Variogram smoothness parameters for space and time.
#' @param adv Advection vector (default = c(0, 0)).
#' @param threshold Threshold for the r-Pareto process (default = 1).
#' @param s0_index Index of the spatial conditioning point.
#' @param t0_index Index of the temporal conditioning point.
#' @param seed Optional random seed for reproducibility.
#' @return A list with simulated field Z, Gaussian process W, and variogram gamma0.
#' @export
sim_rpareto_coords <- function(coords, times,
                              beta1, beta2, alpha1, alpha2,
                              adv = c(0, 0), threshold = 1,
                              s0_index = 1, t0_index = 1, seed = NULL) {

  if (!is.null(seed)) set.seed(seed)

  RandomFields::RFoptions(spConform = FALSE, allow_duplicated_locations = TRUE,
                          install = "no")

  coords <- as.data.frame(coords)
  n_sites <- nrow(coords)
  lt <- length(times)

  # Models: variogram gamma(r)=var*r^alpha
  modelSpace <- RandomFields::RMfbm(alpha = alpha1, var = 2 * beta1, scale = 1)
  modelTime  <- RandomFields::RMfbm(alpha = alpha2, var = 2 * beta2, scale = 1)

  # spatio-temporal grid
  grid <- expand.grid(site = seq_len(n_sites), it = seq_len(lt))
  x  <- coords$Longitude[grid$site]
  y  <- coords$Latitude [grid$site]
  tt <- times[grid$it]

  # Shift by advection
  x_shift <- x - tt * adv[1]
  y_shift <- y - tt * adv[2]

  # Conditioning point
  ind0 <- which(grid$site == s0_index & grid$it == t0_index)

  x0s <- x_shift[ind0]
  y0s <- y_shift[ind0]
  t0  <- tt[ind0]

  # gamma0(s,t) = gamma_space(s-s0) + gamma_time(t-t0)
  gamma_space_vec <- RandomFields::RFvariogram(modelSpace,
                                               x = x_shift - x0s,
                                               y = y_shift - y0s)
  gamma_time_vec  <- RandomFields::RFvariogram(modelTime,
                                               x = tt - t0)
  gamma0_vec <- gamma_space_vec + gamma_time_vec
  gamma0 <- matrix(gamma0_vec, nrow = n_sites, ncol = lt, byrow = FALSE)

  # simulate Ws on the unique shifted coordinates
  # (avoid simulating duplicates)
  key <- paste0(sprintf("%.16g", x_shift), "_", sprintf("%.16g", y_shift))
  key_u <- unique(key)
  map_idx <- match(key, key_u)

  xy_u <- do.call(rbind, strsplit(key_u, "_", fixed = TRUE))
  x_u <- as.numeric(xy_u[, 1])
  y_u <- as.numeric(xy_u[, 2])

  W_s_u <- RandomFields::RFsimulate(modelSpace, x = x_u, y = y_u, grid = FALSE)
  W_s   <- W_s_u[map_idx]

  # simulate Wt 
  W_t <- RandomFields::RFsimulate(modelTime, times, grid = TRUE)
  # map back to full GP
  W <- W_s + W_t[grid$it]
  W_mat <- matrix(W, nrow = n_sites, ncol = lt, byrow = FALSE)

  # r-Pareto transform
  W0 <- W_mat[s0_index, t0_index]
  Y  <- exp(W_mat - W0 - gamma0)
  # R  <- 1 / runif(1) # Pareto(1)
  R = evd::rgpd(n = 1, loc = 1, scale = 1, shape = 1)
  Z <- threshold * R * Y
  rownames(Z) <- rownames(coords)

  list(Z = Z, W = W_mat, gamma0 = gamma0, R = R)
}



#' save_simulations function
#'
#' This function transforms the Brown-Resnick simulations into dataframes and
#' save them into CSV files rain_BR_i.csv.
#'
#' @param simu The BR simulationsas array.
#' @param sites_coords The coordinates of the sites.
#' @param nsimu The number of BR simulations.
#' @param folder The folder path.
#' @param file The filename without extension, default is "rainBR".
#' @param forcedind The index of the simulation to save, default is NA.
#'
#' @return None
#'
#' @export
save_simulations <- function(simu, ngrid, folder, file = "rainBR",
                             forcedind = NA) {
  # Create a data frame with coordinates
  coord_df <- expand.grid(x = 1:ngrid, y = 1:ngrid)
  coord_df$point <- paste0("S", seq_along(coord_df$y))

  simu_long <- as.data.frame.table(simu, responseName = "value")
  names(simu_long) <- c("x", "y", "t", "nsimu", "value")

  simu_long$y <- as.integer(factor(simu_long$y))
  simu_long$x <- as.integer(factor(simu_long$x))
  simu_long$t <- as.integer(factor(simu_long$t))
  simu_long$nsimu <- as.integer(factor(simu_long$nsimu))
  simu_long <- merge(simu_long, coord_df, by = c("x", "y"))

  simulations <- sort(unique(simu_long$nsimu))  # List of available simulations

  # Loop over each simulation
  for (simu in simulations) {
    # Filter simu_long for this simulation
    simu_data <- simu_long %>%
      dplyr::filter(simu_long$nsimu == simu) %>%
      dplyr::select(x, y, t, value, point)


    data_wide <- simu_data %>%
      dplyr::select(t, value, point) %>%
      tidyr::pivot_wider(names_from = point, values_from = value)


    # Add the reshaped data to the list with a unique name based on the simulation number
    df <- as.data.frame(data_wide)
    # sort by t
    df <- df[order(df$t), ]
    # Remove t column
    df <- df[, -1]

    # list_res[[simu]] <- df
    # Define filename and write
    index_to_use <- if (!is.na(forcedind)) forcedind else simu
    filename <- paste0(folder, file, "_", index_to_use, ".csv")
    write.csv(df, file = filename, row.names = FALSE)
  }
}



#' Convert simulations from array to list of data.frames
#'
#' @param simu The simulation array (x, y, t, nsimu)
#' @param ngrid Size of spatial grid (used to build coordinates)
#'
#' @return A list of data.frames, one per episode
#' @export
convert_simulations_to_list <- function(simu, sites_coords) {
  if (length(dim(simu)) == 5) {
    # remove last dimension if present
    simu <- simu[,,, ,1]
  }
  ngrid <- sqrt(nrow(sites_coords))  # n^2 rows => ngrid
  coord_df <- data.frame(
    x = sites_coords$Longitude,
    y = sites_coords$Latitude,
    point = rownames(sites_coords)
  )

  # Flatten the 4D array
  dimnames(simu) <- list(
    x = 1:ngrid,
    y = 1:ngrid,
    t = seq_len(dim(simu)[3]),
    nsimu = seq_len(dim(simu)[4])
  )

  simu_long <- as.data.frame.table(simu, responseName = "value")
  names(simu_long) <- c("x", "y", "t", "nsimu", "value")

  # Ensure factors are integers
  simu_long$y <- as.integer(factor(simu_long$y))
  simu_long$x <- as.integer(factor(simu_long$x))
  simu_long$t <- as.integer(factor(simu_long$t))
  simu_long$nsimu <- as.integer(factor(simu_long$nsimu))

  # Merge to get point names
  simu_long <- merge(simu_long, coord_df, by = c("x", "y"))

  simulations <- sort(unique(simu_long$nsimu))  # simulation indices

  list_res <- vector("list", length(simulations))

  for (i in simulations) {
    df_i <- simu_long %>%
      dplyr::filter(nsimu == i) %>%
      dplyr::select(t, value, point) %>%
      tidyr::pivot_wider(names_from = point, values_from = value) %>%
      arrange(t)

    # remove 't' column, keep only data
    df_i <- df_i[, -1]
    list_res[[i]] <- as.data.frame(df_i)
  }

  return(list_res)
}


#' Convert a single simulation (3D array) to a data.frame
#'
#' @param simu_i A 3D array of dimensions (x, y, t) corresponding to one simulation
#' @param ngrid The size of the spatial grid
#'
#' @return A data.frame with each column being a spatial point and rows as time steps
#' @export
convert_single_simulation_to_df <- function(simu_i, ngrid) {
  coord_df <- expand.grid(x = 1:ngrid, y = 1:ngrid)
  coord_df$point <- paste0("S", seq_len(nrow(coord_df)))

  # Flatten the 3D array into a long data.frame
  simu_long <- as.data.frame.table(simu_i, responseName = "value")
  names(simu_long) <- c("x", "y", "t", "value")

  # Ensure proper types
  simu_long$y <- as.integer(factor(simu_long$y))
  simu_long$x <- as.integer(factor(simu_long$x))
  simu_long$t <- as.integer(factor(simu_long$t))

  # Merge with coordinates to get point labels
  simu_long <- merge(simu_long, coord_df, by = c("x", "y"))

  # Pivot to wide format: one column per spatial point
  df <- simu_long %>%
    dplyr::select(t, value, point) %>%
    tidyr::pivot_wider(names_from = point, values_from = value) %>%
    dplyr::arrange(t)

  # Remove time column and return just the matrix of values
  df <- df[, -1]
  return(as.data.frame(df))
}


# #' sim_episode function
# #' This function simulates a single episode of spatio-temporal data
# #' based on r-Pareto process simulations.
# #'' @param params_vario A list containing variogram parameters:
# #'                     beta1, beta2, alpha1, alpha2.
# #' @param params_margins A list containing marginal parameters:
# #'                       xi, sigma, kappa (vectors per site).
# #' @param x Vector for the first spatial dimension.
# #' @param y Vector for the second spatial dimension.
# #' @param times Vector for the temporal dimension.
# #' @param adv The advection coordinates vector.
# #' @param t0 The conditioning time point.
# #' @param s0 The conditioning spatial location (index).
# #' @param u_s Vector of thresholds per site for the Pareto transformation.
# #' 
# #' @return A 3D array of simulated data (sites x times).
# #' @export
# sim_episode <- function(params_vario, params_margins, coords, times, adv, t0, s0, u_s) {
#   #---------------------------------------------------------
#   # Simulate r-Pareto process on the spatial grid
#   #---------------------------------------------------------
#   sim <- sim_rpareto(
#     beta1 = params_vario$beta1,
#     beta2 = params_vario$beta2,
#     alpha1 = params_vario$alpha1,
#     alpha2 = params_vario$alpha2,
#     adv = adv,
#     x = unique(coords$x),
#     y = unique(coords$y),
#     t = times,
#     t0 = t0,
#     s0 = s0
#   )

#   # Keep only first realization / first conditioning site
#   Z <- sim$Z[,,,1,1, drop = TRUE] * u_s[s0]

#   n_sites <- nrow(coords)
#   nx <- length(unique(coords$x))
#   ny <- length(unique(coords$y))
#   nt <- length(times)

#   U <- array(NA_real_, dim = c(nx, ny, nt))
#   X <- array(NA_real_, dim = c(nx, ny, nt))

#   #---------------------------------------------------------
#   # Transform site by site
#   #---------------------------------------------------------
#   for (k in seq_len(n_sites)) {
#     # Transform Pareto → Uniform using local threshold
#     xk <- coords$x[k]
#     yk <- coords$y[k]
#     U[xk, yk, ] <- pmax(1 - u_s[k] / Z[xk, yk, ], 0)

#     # Transform Uniform → rainfall intensity (EGPD)
#     X[xk, yk, ] <- qextgp(
#       U[xk, yk, ],
#       type  = 1,
#       xi    = params_margins$xi[k],
#       sigma = params_margins$sigma[k],
#       kappa = params_margins$kappa[k]
#     )
#   }

#   # Return rainfall simulations (matrix [site × time])
#   return(X)
# }



# sim_episode_coords <- function(params_vario, params_margins, coords, times, adv, t0, s0, u_s) {
#   #---------------------------------------------------------
#   # Simulate r-Pareto process on the spatial grid
#   #---------------------------------------------------------
#   sim <- sim_rpareto_coords(
#     beta1 = params_vario$beta1,
#     beta2 = params_vario$beta2,
#     alpha1 = params_vario$alpha1,
#     alpha2 = params_vario$alpha2,
#     adv = adv,
#     coords = coords,
#     t = times,
#     t0 = t0,
#     s0 = s0
#   )

#   # Keep only first realization / first conditioning site
#   index_s0 <- which(coords$Latitude == s0[2] & coords$Longitude == s0[1])
#   site_name <- rownames(coords)[index_s0]
#   Y <- sim$Z[,,1, drop = TRUE]
#   Z <- Y * u_s[site_name]
#   n_sites <- nrow(coords)
#   nt <- length(times)
#   U <- matrix(NA_real_, n_sites, nt, dimnames = list(rownames(coords), NULL))
#   X <- matrix(NA_real_, n_sites, nt, dimnames = list(rownames(coords), NULL))
#   #---------------------------------------------------------
#   # Transform site by site
#   #---------------------------------------------------------
#   for (k in seq_len(n_sites)) {
#     # Transform Pareto → Uniform using local threshold
#     sk <- rownames(coords)[k]
#     U[sk, ] <- pmax(1 - u_s[k] / Z[sk, ], 0.00001)

#     # Transform Uniform → rainfall intensity (EGPD)
#     X[sk, ] <- qextgp(
#       U[sk, ],
#       type  = 1,
#       xi    = params_margins$xi[k],
#       sigma = params_margins$sigma[k],
#       kappa = params_margins$kappa[k]
#     )
#   }

#   # Return rainfall simulations (matrix [site × time])
#   return(X)
# }
