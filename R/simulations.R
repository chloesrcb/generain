#' sim_BR function
#'
#' This function performs a simulation of a spatio-temporal Brown-Resnick
#' process using a fractionnal Brownian motion model and based on 
#' David Leber code.
#'
#' @param beta1 The value of beta1.
#' @param beta2 The value of beta2.
#' @param alpha1 The value of alpha1.
#' @param alpha2 The value of alpha2.
#' @param x Vector for the first dimension (spatial x in our case).
#' @param y Vector for the second dimension (spatial y in our case)
#' @param z Vector for the third dimension (time in our case).
#' @param n.BR The number of BR simulations to perform.
#' @param adv The advection coordinates vector. Default is c(0, 0).
#'
#' @return The result of the simulation.
#'
#' @import RandomFields
#' @import stats
#'
#' @export
sim_BR <- function(beta1, beta2, alpha1, alpha2, x, y, z, n.BR, adv = c(0, 0)) {
  # beta1, beta2, alpha1, alpha2 are variogram parameters
  # x is the first dimension (spatial x in our case)
  # y is the second dimension (spatial y in our case)
  # z is the third dimension (time in our case)
  # (adv1, adv2) advection coordinates vector
  ## Setup
  RandomFields::RFoptions(spConform = FALSE, install="no")
  lx <- length(sx <- seq_along(x))  # spatial
  ly <- length(sy <- seq_along(y))  # spatial
  lz <- length(sz <- seq_along(z))  # temporal
  ## Model-Variogram BuhlCklu
  modelBuhlCklu <- RandomFields::RMfbm(alpha = alpha1, var = beta1, proj = 1) +
                   RandomFields::RMfbm(alpha = alpha1, var = beta1, proj = 2) +
                   RandomFields::RMfbm(alpha = alpha2, var = beta2, proj = 3)

  ## Construct grid
  Nxy <- lx * ly
  N <- Nxy * lz
  grid <- matrix(0, nrow=N, ncol=3) # (N,3)-matrix

  for (i in sx)
    for (j in seq_len(ly*lz))
      grid[i+(j-1)*ly, 1] <- i

  for (i in sy)
    for (j in sx)
      for(k in sz)
        grid[j+lx*(i-1)+(k-1)*Nxy, 2] <- i

  for (i in sz)
    for (j in seq_len(Nxy))
      grid[j+Nxy*(i-1), 3] <- i

  ## Construct shifted variogram
  Varm1 <- vapply(seq_len(N), function(n)
    RandomFields::RFvariogram(modelBuhlCklu,
                x = (sx - grid[n, 1]) - adv[1] * grid[n, 3],
                y = (sy - grid[n, 2]) - adv[2] * grid[n, 3],
                z = sz - grid[n, 3]),
    array(NA_real_, dim = c(lx, ly, lz)))

  ## => (lx, ly, lz, N)-array

  ## Main
  Z <- array(, dim = c(lx, ly, lz, n.BR)) # 3d array
  E <- matrix(rexp(n.BR * N), nrow=n.BR, ncol=N)

  for (i in seq_len(n.BR)) {
    ## n=1
    V <- 1 / E[i, 1] # poisson process
    W <- RandomFields::RFsimulate(modelBuhlCklu, x, y, z, n=1) # gaussian process
    Y <- exp(W - W[1] - Varm1[,,, 1])
    Z[,,, i] <- V * Y

    ## n in {2,..,N}
    for (n in 2:N) {
      Exp <- E[i, n]
      V <- 1 / Exp
      while(V > Z[N * (i - 1) + n]) {
        W <- RandomFields::RFsimulate(modelBuhlCklu, x, y, z)
        Y <- exp(W - W[n] - Varm1[,,, n])
        if (all(V * Y[seq_len(n-1)] < Z[(N * (i-1) + 1):(N * (i-1) + (n-1))])) {
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
sim_BR_parallel <- function(params, n.BR, spa, temp){
  # Parallelize the loop
  cl <- makeCluster(detectCores())
  clusterEvalQ(cl, library(RandomFields))
  clusterEvalQ(cl, library(foreach))
  clusterEvalQ(cl, library(doParallel))
  clusterEvalQ(cl, library(parallel))

  simulations <- foreach(i = 1:n.BR, .combine = "c") %dopar% {
        sim_BR(params[1], params[2], params[3], params[4],
                spa, spa, temp, 10)
  }
  stopCluster(cl)

}

#' save_simulations function
#'
#' This function transforms the Brown-Resnick simulations into dataframes and
#' save them into CSV files rain_BR_i.csv.
#'
#' @param BR The BR simulationsas array.
#' @param ngrid The number of grid points.
#' @param n.BR The number of BR simulations.
#' @param path The path to save the dataframes.
#'
#' @return None
#'
#' @export
save_simulations <- function(BR, ngrid, n.BR, path) {
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

    # Save the dataframe in a file
    filename <- paste0(path, "rainBR_", i, ".csv")
    write.csv(df, file = filename, row.names = FALSE)
  }
}