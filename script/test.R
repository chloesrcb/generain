sim_BR <- function(beta1, beta2, alpha1, alpha2, x, y, t, n.BR, adv = c(0, 0)) {
  # beta1, beta2, alpha1, alpha2 are variogram parameters
  # x is the first dimension (spatial x in our case)
  # y is the second dimension (spatial y in our case)
  # t is the third dimension (time in our case)
  # adv is the advection vector
  
  ## Setup
  RandomFields::RFoptions(spConform = FALSE, install = "no")
  lx <- length(sx <- seq_along(x))  # spatial
  ly <- length(sy <- seq_along(y))  # spatial
  lt <- length(st <- seq_along(t))  # temporal
  
  ## Model-Variogram BuhlCklu
  modelBuhlCklu <- RandomFields::RMfbm(alpha = alpha1, var = beta1, proj = 1) +
                   RandomFields::RMfbm(alpha = alpha2, var = beta2, proj = 2)

  ## Construct grid
  Nxy <- lx * ly
  N <- Nxy * lt
  grid <- matrix(0, nrow = N, ncol = 3) # (N,3)-matrix

  # fill the first column with x coordinates
  for (i in sx)
    for (j in seq_len(ly * lt))
      grid[i + (j - 1) * ly, 1] <- i

  # fill the second column with y coordinates
  for (i in sy)
    for (j in sx)
      for (k in st)
        grid[j + lx * (i - 1) + (k - 1) * Nxy, 2] <- i

  # fill the third column with temporal t coordinates
  for (i in st)
    for (j in seq_len(Nxy))
      grid[j + Nxy * (i - 1), 3] <- i

  ## Construct shifted variogram
  Varm1 <- vapply(seq_len(N), function(n) {
    dx <- sx - grid[n, 1] # spatial x lags
    dy <- sy - grid[n, 2] # spatial y lags
    dt <- st - grid[n, 3] # temporal lags
    combi <- expand.grid(dx = dx, dy = dy, dt = dt) # combinations of lags
    combi$dx_adv <- combi$dx - adv[1] * combi$dt # adding advection on x
    combi$dy_adv <- combi$dy - adv[2] * combi$dt # adding advection on y
    norm_h <- sqrt(combi$dx_adv^2 + combi$dy_adv^2)
    # compute variogram for each combination
    result <- RandomFields::RFvariogram(modelBuhlCklu,
                x = norm_h,
                y = combi$dt
              )
    return(result)
  }, array(NA_real_, dim = c(lx, ly, lt)))
  
  ## Main
  Z <- array(, dim = c(lx, ly, lt, n.BR)) # 3d array
  E <- matrix(rexp(n.BR * N), nrow = n.BR, ncol = N)
  
  for (i in seq_len(n.BR)) {
    # Adjusted coordinates for Gaussian process simulation
    x_adv <- x - adv[1] * t
    y_adv <- y - adv[2] * t

    W <- RandomFields::RFsimulate(modelBuhlCklu, x_adv, y_adv, t) # Gaussian process with advection
    
    # for n = 1
    V <- 1 / E[i, 1] # Poisson process
    Y <- exp(W - W[1] - Varm1[,,, 1])
    Z[,,,i] <- V * Y

    # n in {2,..,N}
    for (n in 2:N) {
      Exp <- E[i, n]
      V <- 1 / Exp
      while (V > Z[N * (i - 1) + n]) {
        W <- RandomFields::RFsimulate(modelBuhlCklu, x_adv, y_adv, t) # Gaussian process with advection
        Y <- exp(W - W[n] - Varm1[,,, n])
        if (all(!is.na(Y[seq_len(n - 1)]) & (V * Y[seq_len(n - 1)] < Z[(N * (i - 1) + 1):(N * (i - 1) + (n - 1))]))) {
          Z[,,, i] <- pmax(V * Y, Z[,,, i])
        }
        Exp <- Exp + rexp(1)
        V <- 1 / Exp
      }
    }
  }
  
  ## Plot the results
  plot_spatio_temporal(Z, lx, ly, lt)
  
  ## Return
  Z
}

plot_spatio_temporal <- function(Z, lx, ly, lt) {
  library(ggplot2)
  library(reshape2)
  
  for (i in seq_len(dim(Z)[4])) {
    data <- melt(Z[,,,i])
    colnames(data) <- c("x", "y", "t", "value")
    ggplot(data, aes(x = x, y = y, fill = value)) + 
      geom_tile() +
      facet_wrap(~ t) +
      scale_fill_viridis_c() +
      ggtitle(paste("Realization", i)) +
      theme_minimal()
  }
}

## Example usage:
x <- 1:5
y <- 1:5
t <- 1:10
beta1 <- 1
beta2 <- 1
alpha1 <- 0.5
alpha2 <- 0.5
n.BR <- 1
adv <- c(0.1, 0.2)
result <- sim_BR(beta1, beta2, alpha1, alpha2, x, y, t, n.BR, adv)

plot_spatio_temporal(result, length(x), length(y), length(t))
