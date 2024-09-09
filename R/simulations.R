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

#' sim_BR function
#'
#' This function performs a simulation of a spatio-temporal Brown-Resnick
#' process using a fractionnal Brownian motion model and based on the
#' David Leber code.
#'
#' @param beta1 The value of beta1.
#' @param beta2 The value of beta2.
#' @param alpha1 The value of alpha1.
#' @param alpha2 The value of alpha2.
#' @param x Vector for the first dimension (spatial x in our case).
#' @param y Vector for the second dimension (spatial y in our case)
#' @param t Vector for the third dimension (time in our case).
#' @param n.BR The number of BR simulations to perform.
#' @param adv The advection coordinates vector. Default is c(0, 0).
#'
#' @return The result of the simulation.
#'
#' @import RandomFields
#' @import stats
#'
#' @export
sim_BR <- function(beta1, beta2, alpha1, alpha2, x, y, t, n.BR, adv = c(0, 0)) {
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

  ## Model-Variogram BuhlCklu (fractional Brownian motion)
  modelBuhlCklu <- RandomFields::RMfbm(alpha = alpha1, var = beta1, proj = 1) +
                   RandomFields::RMfbm(alpha = alpha1, var = beta1, proj = 2) +
                   RandomFields::RMfbm(alpha = alpha2, var = beta2, proj = 3)

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

  if (all(adv == c(0, 0))) {
    grid_space <- grid[, 1:2]
    distmat_space <- as.matrix(dist(grid_space))
    grid_time <- grid[, 3]
    distmat_time <- as.matrix(dist(matrix(grid_time, ncol = 1)))
  } else {
    distmat_space <- matrix(0, nrow=N, ncol=N)
    distmat_time <- matrix(0, nrow=N, ncol=N)

    for (i in 1:N) {
      for (j in 1:N) {
        s1 <- c(grid[i, 1], grid[i, 2])
        s2 <- c(grid[j, 1], grid[j, 2])
        t1 <- grid[i, 3]
        t2 <- grid[j, 3]
        distmat_space[i, j] <- dist_adv(s1, s2, t1, t2, adv)
        distmat_time[i, j] <- abs(t1 - t2)
      }
    }
  }

  gamma <- gamma_theta(distmat_space, distmat_time, beta1, beta2,
                 alpha1, alpha2)

  ## Main
  Z <- array(, dim = c(lx, ly, lt, n.BR)) # 3d array
  E <- matrix(rexp(n.BR * N), nrow = n.BR, ncol = N)

  for (i in seq_len(n.BR)) {
    # for n = 1
    V <- 1 / E[i, 1] # poisson process
    W <- RandomFields::RFsimulate(modelBuhlCklu, x, y, t) # GP
    Y <- exp(W - W[1] - gamma[,1])
    Z[,,,i] <- V * Y

    # n in {2,..,N}
    for (n in 2:N) {
      Exp <- E[i, n]
      V <- 1 / Exp
      while(V > Z[N * (i - 1) + n]) {
        W <- RandomFields::RFsimulate(modelBuhlCklu, x, y, t)
        Y <- exp(W - W[n] - gamma[,n])
        if (all(V * Y[seq_len(n - 1)] < Z[(N * (i - 1) + 1):(N * (i - 1) + (n - 1))])) {
          Z[,,, i] <- pmax(V * Y, Z[,,, i])
        }
        Exp <- Exp + rexp(1)
        V <- 1 / Exp
      }
    }
  }
  ## Return
  Z
}

#' sim_BR_adv function
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
#' @param n.BR The number of BR simulations to perform.
#' @param adv The advection coordinates vector. Default is c(0, 0).
#'
#' @return The result of the simulation.
#'
#' @import RandomFields
#' @import stats
#'
#' @export
sim_BR_adv <- function(beta1, beta2, alpha1, alpha2, x, y, z, n.BR,
                       adv) {
  ## Setup
  RandomFields::RFoptions(spConform = FALSE)
  lx <- length(sx <- seq_along(x))
  ly <- length(sy <- seq_along(y))
  lz <- length(sz <- seq_along(z))

  ## Model-Variogram BuhlCklu
  modelBuhlCklu <- RandomFields::RMfbm(alpha = alpha1, var = beta1, proj = 1) +
                   RandomFields::RMfbm(alpha = alpha1, var = beta1, proj = 2) +
                   RandomFields::RMfbm(alpha = alpha2, var = beta2, proj = 3)

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

  # Construct shifted grid with advected coordinates
  grid[, 1] <- grid[, 1] - grid[, 3] * adv[1]
  grid[, 2] <- grid[, 2] - grid[, 3] * adv[2]

  ## Construct shifted variogram
  Varm1 <- vapply(seq_len(N), function(n)
      RandomFields::RFvariogram(modelBuhlCklu,
        x = sx - grid[n, 1],
        y = sy - grid[n, 2],
        z = sz - grid[n, 3]),
        array(NA_real_, dim = c(lx, ly, lz))) ## => (lx, ly, lz, N)-array

  ## Main
  set.seed(123)
  Z <- array(, dim = c(lx, ly, lz, n.BR)) # 4d array
  E <- matrix(rexp(n.BR * N), nrow = n.BR, ncol = N)
  for (i in seq_len(n.BR)) { ## n=1 
    V <- 1 / E[i, 1]
    W <- RandomFields::RFsimulate(modelBuhlCklu, x, y, z, n = 1)
    Y <- exp(W - W[1] - Varm1[, , , 1])
    Z[, , , i] <- V * Y
    ## n in {2,..,N}
    for (n in 2:N) {
      Exp <- E[i, n]
      V <- 1 / Exp 
      while(V > Z[N * (i - 1) + n]) {
        W <- RandomFields::RFsimulate(modelBuhlCklu, x, y, z)
        Y <- exp(W - W[n] - Varm1[, , , n])
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
#' @param n.res The number of simulations to perform.
#' @param adv The advection coordinates vector. Default is c(0, 0).
#'
#' @return The result of the simulation.
#'
#' @import RandomFields
#' @import stats
#'
#' @export
sim_rpareto <- function(beta1, beta2, alpha1, alpha2, x, y, t, n.res,
                        adv = c(0, 0)) {
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
  modelBuhlCklu <- RandomFields::RMfbm(alpha = alpha1, var = beta1, proj = 1) +
                   RandomFields::RMfbm(alpha = alpha1, var = beta1, proj = 2) +
                   RandomFields::RMfbm(alpha = alpha2, var = beta2, proj = 3)

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

  s0_space <- 1 # Spatial conditioning point
  s0_time <- 5  # Temporal conditioning point
  s0 <- s0_space + (s0_time - 1) * lx^2
  # grid[,s0]

  # Spatio-temporal matrix of distances

  if (all(adv == c(0, 0))) {
    grid_space <- grid[, 1:2]
    distmat_space <- as.matrix(dist(grid_space))
    grid_time <- grid[, 3]
    distmat_time <- as.matrix(dist(matrix(grid_time, ncol = 1)))
  } else {
    distmat_space <- matrix(0, nrow=N, ncol=N)
    distmat_time <- matrix(0, nrow=N, ncol=N)

    for (i in 1:N) {
      for (j in 1:N) {
        s1 <- c(grid[i, 1], grid[i, 2])
        s2 <- c(grid[j, 1], grid[j, 2])
        t1 <- grid[i, 3]
        t2 <- grid[j, 3]
        distmat_space[i, j] <- dist_adv(s1, s2, t1, t2, adv)
        distmat_time[i, j] <- abs(t1 - t2)
      }
    }
  }

  gamma <- gamma_theta(distmat_space, distmat_time, beta1, beta2, 
                 alpha1, alpha2)

  # Main
  Z <- array(, dim = c(lx, ly, lt, n.res)) # 3d array
  for (i in seq_len(n.res)) {
    W <- RandomFields::RFsimulate(modelBuhlCklu, x, y, t) # GP
    Y <- exp(W - W[s0] - gamma[, s0])
    R <- evd::rgpd(n = 1, loc = 1, scale = 1, shape = 1)
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
