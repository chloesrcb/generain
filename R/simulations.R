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

  t_index <- grid$t + 1 # index starts at 1

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


#' compute_st_variogram function (supports isotropic and anisotropic variograms)
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

  # Determine isotropic or anisotropic
  is_anisotropic <- !is.null(gamma_space_x) && !is.null(gamma_space_y)
  is_isotropic <- !is.null(gamma_space)

  if (!is_isotropic && !is_anisotropic) {
    stop("Provide either gamma_space for isotropic or gamma_space_x and gamma_space_y for anisotropic case.")
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
  t_index <- grid$t + 1
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

    if (is_anisotropic) {
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
  t_index <- grid$t + 1 # index starts at 1
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
#' @param W_s The spatial Gaussian random field (for isotropic case).
#' @param W_s_x The spatial Gaussian random field in the x direction (for anisotropic case).
#' @param W_s_y The spatial Gaussian random field in the y direction (for anisotropic case).
#' @param W_t The temporal Gaussian random field.
#' @param adv The advection coordinates vector.
#'
#' @return The spatio-temporal Gaussian random field.
#'
#' @export
compute_st_gaussian_process <- function(grid, W_s = NULL, 
                                        W_s_x = NULL, W_s_y = NULL, W_t,
                                        adv) {
  lx <- length(unique(grid$x))
  ly <- length(unique(grid$y))
  lt <- length(unique(grid$t))
  coords <- cbind(grid$shifted_x, grid$shifted_y)
  

  # Determine isotropic or anisotropic
  is_anisotropic <- !is.null(W_s_x) && !is.null(W_s_y)
  is_isotropic <- !is.null(W_s)

  if (!is_isotropic && !is_anisotropic) {
    stop("Provide either W_s for isotropic or W_s_x and W_s_y for anisotropic case.")
  }

  # Remove duplicates if no advection
  if (!any(is.na(adv)) && all(adv == 0)) {
    duplicates <- duplicated(coords)
    filtered_coords <- coords[!duplicates, ]
    coords <- filtered_coords
  }

  if (is.vector(coords)) {
    coords <- matrix(coords, nrow = 1)
  }
  nsites <- nrow(coords)

  W_s_t <- array(NA, dim = c(lx, ly, lt))  # Initialize the 3D array for results
  t_index <- grid$t + 1  # Adjust for indexing (R is 1-based)
  
  # Loop over each grid point
  for (i in seq_along(t_index)) {
    s_x <- grid$x[i]
    s_y <- grid$y[i]
    t_idx <- t_index[i]
    
    # Handle advection: adjust the index based on the grid size
    if (!any(is.na(adv)) && all(adv == 0) && i > nsites) {
      if (i %% nsites == 0) {
        ind_W_s <- i - nsites * (i %/% nsites - 1)
      } else {
        ind_W_s <- i - nsites * (i %/% nsites)
      }
    } else {
      ind_W_s <- i
    }

    # Check if anisotropic or isotropic, and handle accordingly
    if (is_anisotropic) {
      # Anisotropic case: Combine W_s_x and W_s_y
      W_s_t_point <- W_s_x[ind_W_s] + W_s_y[ind_W_s] + W_t[t_index[i]]
    } else {
      # Isotropic case: Use W_s directly
      W_s_t_point <- W_s[ind_W_s] + W_t[t_index[i]]
    }
    # Store the result in the 3D array
    W_s_t[s_x, s_y, t_idx] <- W_s_t_point
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

  ## Variogram
  N <- nrow(grid) # number of points
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
    ## n in {2,..,N}
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
  ## Return
  Z
}

#' sim_BR_aniso function
#' 
#' Simulate the spatio-temporal Brown-Resnick process with anisotropy and 
#' advection.
#' 
#' @param beta1 The value of beta1.
#' @param beta2 The value of beta2.
#' @param alpha1 The value of alpha1.
#' @param alpha2 The value of alpha2.
#' @param x Vector for the first dimension (spatial x in our case).
#' @param y Vector for the second dimension (spatial y in our case)
#' @param z Vector for the third dimension (time in our case).
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
sim_BR_aniso <- function(beta1, beta2, alpha1, alpha2, x, y, z, adv = NA,
                          nres = 1) {
  ## Setup
  RandomFields::RFoptions(spConform = FALSE)
  lx <- length(sx <- seq_along(x))
  ly <- length(sy <- seq_along(y))
  lz <- length(sz <- seq_along(z))

  ## Model-Variogram BuhlCklu
  modelBuhlCklu <- RandomFields::RMfbm(alpha = alpha1, var = 2*beta1, proj = 1) +
                   RandomFields::RMfbm(alpha = alpha1, var = 2*beta1, proj = 2) +
                   RandomFields::RMfbm(alpha = alpha2, var = 2*beta2, proj = 3)

  ## Construct grid
  Nxy <- lx * ly
  N <- Nxy * lz
  grid <- matrix(0, nrow = N, ncol = 3) # (N,3)-matrix

  for (i in sx)
    for (j in seq_len(ly * lz))
      grid[i + (j - 1) * ly, 1] <- i

  for (i in sy)
    for (j in sx)
      for (k in sz)
        grid[j + lx * (i - 1) + (k - 1) * Nxy, 2] <- i

  for (i in sz)
    for (j in seq_len(Nxy))
      grid[j + Nxy * (i - 1), 3] <- i

  # # Construct shifted grid with advected coordinates
  # grid[, 1] <- grid[, 1] - grid[, 3] * adv[1]
  # grid[, 2] <- grid[, 2] - grid[, 3] * adv[2]

  ## Construct shifted variogram
  if (all(adv == c(0, 0))) {
    gamma <- vapply(seq_len(N), function(n)
        RandomFields::RFvariogram(modelBuhlCklu,
          x = sx - grid[n, 1],
          y = sy - grid[n, 2],
          z = sz - grid[n, 3]),
          array(NA_real_, dim = c(lx, ly, lz))) ## => (lx, ly, lz, N)-array
  } else { # with advection, a lot longer
    gamma <- advected_variogram(x, y, z, grid, modelBuhlCklu, adv)
  }

  ## Main
  Z <- array(, dim = c(lx, ly, lz, nres)) # 4d array
  E <- matrix(rexp(nres * N), nrow = nres, ncol = N)
  for (i in seq_len(nres)) { ## n=1
    V <- 1 / E[i, 1]
    W <- RandomFields::RFsimulate(modelBuhlCklu, x, y, z, n = 1)
    Y <- exp(W - W[1] - gamma[, , , 1])
    Z[, , , i] <- V * Y
    ## n in {2,..,N}
    for (n in 2:N) {
      Exp <- E[i, n]
      V <- 1 / Exp 
      while (V > Z[N * (i - 1) + n]) {
        W <- RandomFields::RFsimulate(modelBuhlCklu, x, y, z)
        Y <- exp(W - W[n] - gamma[, , , n])
        if(all(V * Y[seq_len(n-1)] < Z[(N*(i-1)+1):(N*(i-1)+(n-1))]))
          Z[, , , i] <- pmax(V * Y, Z[, , , i])
          Exp <- Exp + rexp(1)
          V <- 1 / Exp
      }
    }
  }
  ## Return
  Z
}


#' sim_rpareto function
#'
#' This function performs a simulation of a spatio-temporal r-Pareto process
#' using a fractionnal Brownian motion model and based on the David Leber code.
#'
#' @param beta1 The value of beta1.
#' @param beta2 The value of beta2.
#' @param alpha1 The value of alpha1.
#' @param alpha2 The value of alpha2.
#' @param x Vector for the first dimension (spatial x in our case).
#' @param y Vector for the second dimension (spatial y in our case)
#' @param t Vector for the third dimension (time in our case).
#' @param adv The advection coordinates vector. Default is c(0, 0).
#' @param t0 Conditional temporal time. Default is 1.
#' @param nres The number of simulations to perform. Default is 1.
#' @param random_s0 Logical value indicating whether to choose a random s0.
#'                  Default is FALSE.
#' @param s0 Vector of dimension 2 for the spatial conditioning point.
#'           Default is c(1, 1).
#' @param s0_radius The radius for random s0 selection. Default is Inf.
#'
#' @return The result of the simulation.
#'
#' @import stats
#' @import RandomFields
#' @import RandomFieldsUtils
#'
#' @export
sim_rpareto <- function(beta1, beta2, alpha1, alpha2, x, y, t,
                        adv = c(0, 0), t0 = 0, nres = 1,
                        random_s0 = FALSE, s0 = c(1, 1),
                        s0_radius = Inf) {
  # beta1, beta2, alpha1, alpha2 are variogram parameters
  # x is the first dimension (spatial x in our case)
  # y is the second dimension (spatial y in our case)
  # z is the third dimension (time in our case)
  # (adv1, adv2) advection coordinates vector
  ## Setups 
  # RandomFields::RFoptions(spConform = FALSE, install = "no")
  # RandomFields::RFoptions(spConform = FALSE, allow_duplicated_locations = T)

  lx <- length(sx <- seq_along(x))  # spatial
  ly <- length(sy <- seq_along(y))  # spatial
  lt <- length(st <- seq_along(t))  # temporal
  site_names <- paste0("S", seq_len(lx * ly))

  ## Model-Variogram BuhlCklu
  modelSpace <- RandomFields::RMfbm(alpha = alpha1, var = 2 * beta1)
  modelTime <- RandomFields::RMfbm(alpha = alpha2, var = 2 * beta2)

  ## Construct grid
  grid_with_advection <- expand.grid(
    x = seq_len(lx),
    y = seq_len(ly),
    t = t
  )

  grid_with_advection$shifted_x <- grid_with_advection$x -
                                    grid_with_advection$t * adv[1]
  grid_with_advection$shifted_y <- grid_with_advection$y -
                                    grid_with_advection$t * adv[2]

  grid <- grid_with_advection

  grid$site <- rep(site_names, times = lt)
  coords <- grid[, 4:5]

  if (all(adv == 0)) {
    duplicates <- duplicated(coords)
    filtered_coords <- coords[!duplicates, ]
    coords <- filtered_coords
  }
  # Possible random s0
  s0_center <- s0
  if (random_s0) {
    grid_points <- expand.grid(x = seq_len(lx), y = seq_len(ly))
    distances <- sqrt((grid_points$x - s0_center[1])^2 + (grid_points$y -
                                                            s0_center[2])^2)
    candidate_points <- grid_points[distances <= s0_radius, ]
    if (nrow(candidate_points) == 0) {
      stop("No grid points found within specified radius of s0_center.")
    }
  }


  s0_list <- list()

  # Main
  Z <- array(, dim = c(lx, ly, lt, nres)) # 4d array
  for (i in seq_len(nres)) {
    # Choose random s from grid
    if (random_s0) {
      if (!is.null(s0_radius) && !is.null(s0_center)) {
        selected_index <- sample(nrow(candidate_points), 1)
        s0 <- as.integer(candidate_points[selected_index, ])
      } else {
        s0 <- c(sample(seq_len(lx), 1), sample(seq_len(ly), 1))
      }
    }
    s0 <- data.frame(x = s0[1], y = s0[2])
    s0_list <- c(s0_list, list(s0))
    # s0 <- c(1, 1) # default


    ## Variogram for s0, t0
    ind_s0_t0 <- which(grid$x == s0$x & grid$y == s0$y &
                              grid$t == t0)

    gamma_space <- RandomFields::RFvariogram( # for s0,t0
        modelSpace,
        x = coords$shifted_x - grid$shifted_x[ind_s0_t0], # s - s0
        y = coords$shifted_y - grid$shifted_y[ind_s0_t0],
      )

    gamma_temp <- RandomFields::RFvariogram( # for t0
        modelTime,
        x = t - grid$t[ind_s0_t0] # t-t0
      )


    # Get gamma spatio-temporal for s0, t0
    gamma_0 <- compute_st_variogram(
      grid,
      gamma_space = gamma_space,
      gamma_temp = gamma_temp,
      adv = adv
    )


    # Spatial gaussian random field on shifted coords
    W_s <- RandomFields::RFsimulate(modelSpace, coords[, 1], coords[, 2],
                                      grid = FALSE)
    # Temporal gaussian random field
    W_t <- RandomFields::RFsimulate(modelTime, t, n = 1, grid = TRUE)
    # Spatio-temporal random field
    W <- compute_st_gaussian_process(grid,
                                     W_s = W_s, W_t = W_t,
                                     adv = adv)
    Y <- exp(W - W[ind_s0_t0] - gamma_0)
    R <- evd::rgpd(n = 1, loc = 1, scale = 1, shape = 1) # simple Pareto
    Z[,,, i] <- R * Y
  }
  # Return
  # Z
  return(list(
    Z = Z,
    s0_used = s0_list
  ))
}


#' sim_rpareto function
#'
#' This function performs a simulation of a spatio-temporal r-Pareto process
#' using a fractionnal Brownian motion model and based on the David Leber code.
#'
#' @param beta1 The value of beta1.
#' @param beta2 The value of beta2.
#' @param alpha1 The value of alpha1.
#' @param alpha2 The value of alpha2.
#' @param x Vector for the first dimension (spatial x in our case).
#' @param y Vector for the second dimension (spatial y in our case)
#' @param t Vector for the third dimension (time in our case).
#' @param adv The advection coordinates vector. Default is c(0, 0).
#' @param t0 Conditional temporal time. Default is 1.
#' @param nres The number of simulations to perform. Default is 1.
#' @param random_s0 Logical value indicating whether to choose a random s0.
#'                  Default is FALSE.
#' @param s0 Vector of dimension 2 for the spatial conditioning point.
#'           Default is c(1, 1).
#' @param s0_radius The radius for random s0 selection. Default is Inf.
#'
#' @return The result of the simulation.
#'
#' @import stats
#' @import RandomFields
#' @import RandomFieldsUtils
#'
#' @export
sim_rpareto_dir <- function(beta1, beta2, alpha1, alpha2, x, y, t,
                        adv = c(0, 0), t0 = 0, nres = 1,
                        random_s0 = FALSE, s0 = c(1, 1),
                        s0_radius = Inf) {
  # beta1, beta2, alpha1, alpha2 are variogram parameters
  # x is the first dimension (spatial x in our case)
  # y is the second dimension (spatial y in our case)
  # z is the third dimension (time in our case)
  # (adv1, adv2) advection coordinates vector
  ## Setups 
  # RandomFields::RFoptions(spConform = FALSE, install = "no")
  RandomFields::RFoptions(spConform = FALSE, allow_duplicated_locations = T)

  lx <- length(sx <- seq_along(x))  # spatial
  ly <- length(sy <- seq_along(y))  # spatial
  lt <- length(st <- seq_along(t))  # temporal
  site_names <- paste0("S", seq_len(lx * ly))

  ## Model-Variogram BuhlCklu
  modelSpace <- RandomFields::RMfbm(alpha = alpha1, var = 2 * beta1)
  modelTime <- RandomFields::RMfbm(alpha = alpha2, var = 2 * beta2)

  ## Construct grid
  grid_with_advection <- expand.grid(
    x = seq_len(lx),
    y = seq_len(ly),
    t = t
  )

  grid_with_advection$shifted_x <- grid_with_advection$x -
                                    grid_with_advection$t * adv[1]
  grid_with_advection$shifted_y <- grid_with_advection$y -
                                    grid_with_advection$t * adv[2]

  grid <- grid_with_advection

  grid$site <- rep(site_names, times = lt)
  coords <- grid[, 4:5]

  if (all(adv == 0)) {
    duplicates <- duplicated(coords)
    filtered_coords <- coords[!duplicates, ]
    coords <- filtered_coords
  }
  # Possible random s0
  s0_center <- s0
  if (random_s0) {
    grid_points <- expand.grid(x = seq_len(lx), y = seq_len(ly))
    distances <- sqrt((grid_points$x - s0_center[1])^2 + (grid_points$y -
                                                            s0_center[2])^2)
    candidate_points <- grid_points[distances <= s0_radius, ]
    if (nrow(candidate_points) == 0) {
      stop("No grid points found within specified radius of s0_center.")
    }
  }


  s0_list <- list()

  # Main
  Z <- array(, dim = c(lx, ly, lt, nres)) # 4d array
  for (i in seq_len(nres)) {
    # Choose random s from grid
    if (random_s0) {
      if (!is.null(s0_radius) && !is.null(s0_center)) {
        selected_index <- sample(nrow(candidate_points), 1)
        s0 <- as.integer(candidate_points[selected_index, ])
      } else {
        s0 <- c(sample(seq_len(lx), 1), sample(seq_len(ly), 1))
      }
    }
    s0 <- data.frame(x = s0[1], y = s0[2])
    s0_list <- c(s0_list, list(s0))
    # s0 <- c(1, 1) # default


    ## Variogram for s0, t0
    ind_s0_t0 <- which(grid$x == s0$x & grid$y == s0$y & grid$t == t0)

    gamma_space_x <- RandomFields::RFvariogram( # for s0,t0
        modelSpace,
        x = coords$shifted_x - grid$shifted_x[ind_s0_t0]
      )

    gamma_space_y <- RandomFields::RFvariogram( # for s0,t0
        modelSpace,
        x = coords$shifted_y - grid$shifted_y[ind_s0_t0],
      )

    gamma_temp <- RandomFields::RFvariogram( # for t0
        modelTime,
        x = t - grid$t[ind_s0_t0] # t-t0
      )


    # Get gamma spatio-temporal for s0, t0
    gamma_0 <- compute_st_variogram(grid,
                               gamma_space_x = gamma_space_x,
                               gamma_space_y = gamma_space_y,
                               gamma_temp = gamma_temp,
                               adv = adv)

    # Spatial gaussian random field on shifted coords
    W_s_x <- RandomFields::RFsimulate(modelSpace, coords[, 1],
                                    grid = FALSE)
    W_s_y <- RandomFields::RFsimulate(modelSpace, coords[, 2],
                                    grid = FALSE)
    # Temporal gaussian random field
    W_t <- RandomFields::RFsimulate(modelTime, t, n = 1, grid = TRUE)
    # Spatio-temporal random field
    W <- compute_st_gaussian_process(grid = grid, W_s_x = W_s_x,
                                    W_s_y = W_s_y, W_t = W_t, adv = adv)
    Y <- exp(W - W[ind_s0_t0] - gamma_0)
    R <- evd::rgpd(n = 1, loc = 1, scale = 1, shape = 1) # simple Pareto
    Z[,,, i] <- R * Y
  }
  # Return
  # Z
  return(list(
    Z = Z,
    s0_used = s0_list
  ))
}



#' sim_rpareto function
#'
#' This function simulates a spatio-temporal r-Pareto process using either 
#' an isotropic or anisotropic fractional Brownian motion model. It is a unified 
#' version of sim_rpareto (isotropic) and sim_rpareto_dir (anisotropic).
#'
#' @param beta1, beta2 Variogram scale parameters for space and time.
#' @param alpha1, alpha2 Variogram smoothness parameters for space and time.
#' @param x, y, t Vectors representing spatial (x, y) and temporal (t) grids.
#' @param adv Advection vector (default = c(0, 0)).
#' @param t0 Time point at which the process is conditioned.
#' @param nres Number of simulations to perform.
#' @param random_s0 If TRUE, selects conditioning point s0 randomly within radius.
#' @param s0 Conditioning spatial location (default = c(1, 1)).
#' @param s0_radius Radius used if random_s0 is TRUE.
#' @param anisotropic If TRUE, uses directional (anisotropic) model.
#'
#' @return A list with simulated field Z and list of conditioning points s0_used.
#' @export
sim_rpareto <- function(beta1, beta2, alpha1, alpha2, x, y, t,
                                adv = c(0, 0), t0 = 0, nres = 1,
                                random_s0 = FALSE, s0 = c(1, 1),
                                s0_radius = Inf,
                                anisotropic = FALSE) {
  # Ensure RandomFields works with duplicated coordinates if needed
  RandomFields::RFoptions(spConform = FALSE, allow_duplicated_locations = TRUE)

  # Dimensions
  lx <- length(x)
  ly <- length(y)
  lt <- length(t)
  site_names <- paste0("S", seq_len(lx * ly))

  # Define fractional Brownian motion models
  modelSpace <- RandomFields::RMfbm(alpha = alpha1, var = 2 * beta1)
  modelTime  <- RandomFields::RMfbm(alpha = alpha2, var = 2 * beta2)

  # Create spatio-temporal grid with advection
  grid <- expand.grid(x = seq_len(lx), y = seq_len(ly), t = t)
  grid$shifted_x <- grid$x - grid$t * adv[1]
  grid$shifted_y <- grid$y - grid$t * adv[2]
  grid$site <- rep(site_names, times = lt)
  coords <- grid[, c("shifted_x", "shifted_y")]

  # Remove duplicate coordinates if no advection
  if (all(adv == 0)) {
    coords <- coords[!duplicated(coords), ]
  }

  # If random s0, precompute possible points within radius
  s0_center <- s0
  if (random_s0) {
    grid_points <- expand.grid(x = seq_len(lx), y = seq_len(ly))
    distances <- sqrt((grid_points$x - s0_center[1])^2 + 
                      (grid_points$y - s0_center[2])^2)
    candidate_points <- grid_points[distances <= s0_radius, ]
    if (nrow(candidate_points) == 0) {
      stop("No grid points found within specified radius of s0_center.")
    }
  }

  s0_list <- list()  # To store all s0 used
  Z <- array(NA, dim = c(lx, ly, lt, nres))  # 4D output array

  for (i in seq_len(nres)) {
    # Select s0 (random or fixed)
    if (random_s0) {
      selected_index <- sample(nrow(candidate_points), 1)
      s0 <- as.integer(candidate_points[selected_index, ])
    }
    s0 <- data.frame(x = s0[1], y = s0[2])
    s0_list[[i]] <- s0

    # Identify index in grid for conditioning point at time t0
    ind_s0_t0 <- which(grid$x == s0$x & grid$y == s0$y & grid$t == t0)

    # Temporal variogram centered at t0
    gamma_temp <- RandomFields::RFvariogram(modelTime, x = t - grid$t[ind_s0_t0])

    if (anisotropic) {
      # Anisotropic: separate space into x and y directions
      gamma_space_x <- RandomFields::RFvariogram(modelSpace,
                                x = coords$shifted_x - grid$shifted_x[ind_s0_t0])
      gamma_space_y <- RandomFields::RFvariogram(modelSpace,
                                x = coords$shifted_y - grid$shifted_y[ind_s0_t0])

      # Combine variograms
      gamma_0 <- compute_st_variogram(
        grid,
        gamma_space_x = gamma_space_x,
        gamma_space_y = gamma_space_y,
        gamma_temp = gamma_temp,
        adv = adv
      )

      # Simulate independent Gaussian fields in each spatial direction
      W_s_x <- RandomFields::RFsimulate(modelSpace, coords[, 1], grid = FALSE)
      W_s_y <- RandomFields::RFsimulate(modelSpace, coords[, 2], grid = FALSE)
      W_t   <- RandomFields::RFsimulate(modelTime, t, n = 1, grid = TRUE)

      # Combine spatial and temporal processes
      W <- compute_st_gaussian_process(
        grid, W_s_x = W_s_x, W_s_y = W_s_y, W_t = W_t, adv = adv
      )
    } else {
      # Isotropic case: compute single spatial variogram
      gamma_space <- RandomFields::RFvariogram(modelSpace,
                          x = coords$shifted_x - grid$shifted_x[ind_s0_t0],
                          y = coords$shifted_y - grid$shifted_y[ind_s0_t0])

      gamma_0 <- compute_st_variogram(
        grid,
        gamma_space = gamma_space,
        gamma_temp = gamma_temp,
        adv = adv
      )

      # Simulate isotropic spatial Gaussian field
      W_s <- RandomFields::RFsimulate(modelSpace, coords[, 1], coords[, 2], grid = FALSE)
      W_t <- RandomFields::RFsimulate(modelTime, t, n = 1, grid = TRUE)

      # Combine spatial and temporal processes
      W <- compute_st_gaussian_process(
        grid, W_s = W_s, W_t = W_t, adv = adv
      )
    }

    # Construct Pareto process
    Y <- exp(W - W[ind_s0_t0] - gamma_0)  # Normalize and shift
    R <- evd::rgpd(n = 1, loc = 1, scale = 1, shape = 1)  # Generate radial component
    Z[,,, i] <- R * Y  # Final field
  }

  return(list(Z = Z, s0_used = s0_list))
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
    colnames(df)

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
convert_simulations_to_list <- function(simu, ngrid) {
  coord_df <- expand.grid(x = 1:ngrid, y = 1:ngrid)
  coord_df$point <- paste0("S", seq_len(nrow(coord_df)))

  # Flatten the 4D array
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

