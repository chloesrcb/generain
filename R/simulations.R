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

  ## Construct shifted variogram
  Varm1 <- vapply(seq_len(N), function(n) {
    dx <- sx - grid[n, 1] # spatial lags
    dy <- sy - grid[n, 2] # spatial lags
    dt <- st - grid[n, 3] # temporal lags
    # norm_h <- sqrt(combi$dx_adv^2 + combi$dy_adv^2)
    # compute variogram for each combination
    if (all(adv == c(0, 0))) { # without advection
      result <- RandomFields::RFvariogram(modelBuhlCklu,
                  x = dx,
                  y = dy,
                  z = dt
                )
    } else { # with advection
      combi <- expand.grid(dx = dx, dy = dy, dt = dt) # combinations of lags
      combi$dx_adv <- combi$dx - adv[1] * combi$dt # adding advection on x
      combi$dy_adv <- combi$dy - adv[2] * combi$dt # adding advection on y
      result <- RandomFields::RFvariogram(modelBuhlCklu,
                  x = combi$dx_adv,
                  y = combi$dy_adv,
                  z = combi$dt
                )
    }
    return(result)
    }, array(NA_real_, dim = c(lx, ly, lt)))

  ## Main
  Z <- array(, dim = c(lx, ly, lt, n.BR)) # 3d array
  E <- matrix(rexp(n.BR * N), nrow = n.BR, ncol = N)

  for (i in seq_len(n.BR)) {
    # for n = 1
    V <- 1 / E[i, 1] # poisson process
    W <- RandomFields::RFsimulate(modelBuhlCklu, x, y, t) # GP
    Y <- exp(W - W[1] - Varm1[,,, 1])
    Z[,,,i] <- V * Y

    # n in {2,..,N}
    for (n in 2:N) {
      Exp <- E[i, n]
      V <- 1 / Exp
      while(V > Z[N * (i - 1) + n]) {
        W <- RandomFields::RFsimulate(modelBuhlCklu, x, y, t)
        Y <- exp(W - W[n] - Varm1[,,, n])
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


#' sim_BR_parallel function
#'
#' This function performs parallel simulations using the BR method.
#'
#' @param params The parameters for the simulation.
#' @param n.BR The number of BR simulations to perform.
#' @param spa The spatial data for the simulation.
#' @param temp The temporal data for the simulation.
#'
#' @return The results of the parallel simulations.
#'
#' @import RandomFields
#' @import foreach
#' @import doParallel
#'
#' @export
sim_BR_parallel <- function(params, num_iterations, spa, temp, folder) {
  # Parallelize the loop
  cl <- makeCluster(detectCores())
  clusterEvalQ(cl, library(RandomFields))
  clusterEvalQ(cl, library(foreach))
  clusterEvalQ(cl, library(doParallel))
  clusterEvalQ(cl, library(parallel))

  ngrid <- length(spa)

  simulations <- foreach(i = 1:num_iterations, .combine = "c") %dopar% {
      BR <- sim_BR(params[1], params[2], params[3], params[4], spa, spa, temp,
                   1, params[5:6])
      save_simulations(BR, ngrid, 1,
                  folder = folder,
                  file = paste0("br_", ngrid^2, "s_",
                                length(temp), "t"), forcedind = i)
  }
  stopCluster(cl)

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
  grid <- expand.grid(x = x, y = y, t = t)

  ## Construct shifted variogram
  Varm1 <- vapply(seq_len(N), function(n) {
    dx <- sx - grid[n, 1] # spatial lags
    dy <- sy - grid[n, 2] # spatial lags
    dt <- st - grid[n, 3] # temporal lags
    combi <- expand.grid(dx = dx, dy = dy, dt = dt) # combinations of lags
    combi$dx_adv <- combi$dx - adv[1] * combi$dt # adding advection on x
    combi$dy_adv <- combi$dy - adv[2] * combi$dt # adding advection on y
    # norm_h <- sqrt(combi$dx_adv^2 + combi$dy_adv^2)
    # compute variogram for each combination
    result <- RandomFields::RFvariogram(modelBuhlCklu,
                x = combi$dx_adv,
                y = combi$dy_adv,
                z = combi$dt
              )
    return(result)
    }, array(NA_real_, dim = c(lx, ly, lt)))
  ## => (lx, ly, lt, N)-array
  
  # Main
  Z <- array(, dim = c(lx, ly, lt, n.res)) # 3d array
  for (i in seq_len(n.res)) {
    W <- RandomFields::RFsimulate(modelBuhlCklu, x, y, t) # GP
    Y <- exp(W - W[1] - Varm1[,,, 1])
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
#' @param BR The BR simulationsas array.
#' @param ngrid The number of grid points.
#' @param n.BR The number of BR simulations.
#' @param folder The folder path.
#' @param file The filename without extension, default is "rainBR".
#' @param forcedind The index of the simulation to save, default is NA.
#'
#' @return None
#'
#' @export
save_simulations <- function(BR, ngrid, n.BR, folder, file = "rainBR",
                             forcedind = NA) {
  # Initialize the list to store the dataframes
  list_dataframes <- list()

  # Loop over the simulations
  for (i in 1:n.BR) {
    # Extract values corresponding to simulation i
    values <- BR[, , , i]

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
