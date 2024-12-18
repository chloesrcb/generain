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
#' @import RandomFields
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
#' @import RandomFields
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
    gamma[grid$y[i], grid$x[i], t_index[i]] <- vario
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
#' @import RandomFields
#' @import stats
#'
#' @export
sim_BR <- function(beta1, beta2, alpha1, alpha2, x, y, t, adv = c(0, 0), 
                   nres = 1) {
  ## Setup
  RandomFields::RFoptions(spConform = FALSE)
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
        # Spatio-temporal random field
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
#' @import RandomFields
#' @import stats
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
#' @param s0 Conditional vector spatial point. Default is c(1, 1).
#' @param t0 Conditional temporal time. Default is 1.
#' @param nres The number of simulations to perform. Default is 1.
#'
#' @return The result of the simulation.
#'
#' @import RandomFields
#' @import stats
#'
#' @export
sim_rpareto <- function(beta1, beta2, alpha1, alpha2, x, y, t,
                        adv = c(0, 0), s0 = c(1, 1), t0 = 0, nres = 1) {
  # beta1, beta2, alpha1, alpha2 are variogram parameters
  # x is the first dimension (spatial x in our case)
  # y is the second dimension (spatial y in our case)
  # z is the third dimension (time in our case)
  # (adv1, adv2) advection coordinates vector
  ## Setups
  RandomFields::RFoptions(spConform = FALSE, install = "no")
  lx <- length(sx <- seq_along(x))  # spatial
  ly <- length(sy <- seq_along(y))  # spatial
  lt <- length(st <- seq_along(t))  # temporal

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

  ## Variogram for s0, t0
  ind_s0_t0 <- which(grid$x == s0[1] & grid$y == s0[2] &
                            grid$t == t0)

  gamma_space <- RandomFields::RFvariogram( # for s0,t0
      modelSpace,
      x = coords$shifted_x - grid$shifted_x[ind_s0_t0],
      y = coords$shifted_y - grid$shifted_y[ind_s0_t0],
    )

  gamma_temp <- RandomFields::RFvariogram( # for t0
      modelTime,
      x = grid$t[ind_s0_t0] - t # Todo: verif t0-t or t-t0 for t0 != 0
    )


  # Get gamma spatio-temporal for s0, t0
  gamma_0 <- compute_gamma_point(grid, gamma_space, gamma_temp, adv)


  # Main
  Z <- array(, dim = c(lx, ly, lt, nres)) # 4d array
  for (i in seq_len(nres)) {
    # Spatial gaussian random field on shifted coords
    W_s <- RandomFields::RFsimulate(modelSpace, coords[, 1], coords[, 2],
                                      grid = FALSE)
    # Temporal gaussian random field
    W_t <- RandomFields::RFsimulate(modelTime, t, n = 1, grid = TRUE)
    # Spatio-temporal random field
    W <- compute_W_s_t(grid, W_s, W_t, adv)
    Y <- exp(W - W[ind_s0_t0] - gamma_0)
    R <- evd::rgpd(n = 1, loc = 1, scale = 1, shape = 1) # simple Pareto
    Z[,,, i] <- R * Y
  }
  # Return
  Z
}

#' Simulate the spatio-temporal r-Pareto process with anisotropy
#' and advection.
#' 
#' @param beta1 The value of beta1.
#' @param beta2 The value of beta2.
#' @param alpha1 The value of alpha1.
#' @param alpha2 The value of alpha2.
#' @param x Vector for the first dimension (spatial x in our case).
#' @param y Vector for the second dimension (spatial y in our case)
#' @param t Vector for the third dimension (time in our case).
#' @param adv The advection coordinates vector. Default is c(0, 0).
#' 
#' @return The result of the simulation.
#' 
#' @import RandomFields
#' @import stats
#' 
#' @export
sim_rpareto_aniso <- function(beta1, beta2, alpha1, alpha2, x, y, t,
                        adv = c(0, 0), s0 = c(1, 1), t0 = 1, nres = 1) {
  # beta1, beta2, alpha1, alpha2 are variogram parameters
  # x is the first dimension (spatial x in our case)
  # y is the second dimension (spatial y in our case)
  # z is the third dimension (time in our case)
  # (adv1, adv2) advection coordinates vector
  ## Setup
  RandomFields::RFoptions(spConform = FALSE, install = "no")
  lx <- length(sx <- seq_along(x))  # spatial
  ly <- length(sy <- seq_along(y))  # spatial
  lt <- length(st <- seq_along(t))  # temporal

  ## Model-Variogram BuhlCklu
  modelBuhlCklu <- RandomFields::RMfbm(alpha = alpha1, var = 2*beta1, proj = 1) +
                   RandomFields::RMfbm(alpha = alpha1, var = 2*beta1, proj = 2) +
                   RandomFields::RMfbm(alpha = alpha2, var = 2*beta2, proj = 3)

  ## Construct grid
  Nxy <- lx * ly # spatial grid size
  N <- Nxy * lt # spatio-temporal grid size
  grid <- matrix(0, nrow=N, ncol=3) # (N,3)-matrix

  for (i in sx)
    for (j in seq_len(ly*lt))
      grid[i+(j-1)*ly, 1] <- i

  for (i in sy)
    for (j in sx)
      for(k in st)
        grid[j+lx*(i-1)+(k-1)*Nxy, 2] <- i

  for (i in st)
    for (j in seq_len(Nxy))
      grid[j+Nxy*(i-1), 3] <- i

  # Construct shifted variogram for conditional spatio-temporal point
  gamma <-  conditional_variogram(x, y, t, s0, t0, grid, modelBuhlCklu, adv)
  # params <- c(beta1, beta2, alpha1, alpha2, adv[1], adv[2])
  # gamma <- conditional_variogram_hand(x, y, t, s0, t0, params)

  # Main
  Z <- array(, dim = c(lx, ly, lt, nres)) # 3d array
  for (i in seq_len(nres)) {
    W <- RandomFields::RFsimulate(modelBuhlCklu, x, y, t) # GP
    Y <- exp(W - W[s0[2], s0[1], t0] - gamma)
    R <- evd::rgpd(n = 1, loc = 1, scale = 1, shape = 1) # simple Pareto
    Z[,,, i] <- R * Y
  }
  # Return
  Z
}


#' save_simulations function
#'
#' This function transforms the Brown-Resnick simulations into dataframes and
#' save them into CSV files rain_BR_i.csv.
#'
#' @param simu The BR simulationsas array.
#' @param ngrid The number of grid points.
#' @param nsimu The number of BR simulations.
#' @param folder The folder path.
#' @param file The filename without extension, default is "rainBR".
#' @param forcedind The index of the simulation to save, default is NA.
#'
#' @return None
#'
#' @export
save_simulations <- function(simu, ngrid, nsimu, folder, file = "rainBR",
                             forcedind = NA) {
  # Initialize the list to store the dataframes
  list_dataframes <- list()

  # Loop over the simulations
  for (i in 1:nsimu) {
    # Extract values corresponding to simulation i
    values <- simu[, , , i]

    # Initialize the dataframe with the first column
    df <- data.frame(values[1, 1, ])

    # Loop over the other columns for each site
    for (row in 1:ngrid) {
      for (col in 1:ngrid) {
        # Exclude the first column already added
        if (!(row == 1 && col == 1)) {
          v <- c(values[row, col, ])  # Values corresponding to this column
          df <- cbind(df, v)  # Add the column to the dataframe
        }
      }
    }

    colnames(df) <- paste0("S", 1:(ngrid^2))  # Rename the columns
    # Add the dataframe to the list
    list_dataframes[[i]] <- df
    if (!is.na(forcedind)) {
      i <- forcedind
    }
    # Save the dataframe in a file
    filename <- paste0(folder, file, "_", i, ".csv")
    write.csv(df, file = filename, row.names = FALSE)
  }
}
