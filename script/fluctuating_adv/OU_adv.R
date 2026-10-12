

# Simulate a 2D Ornstein-Uhlenbeck wind trajectory
simulate_OU_advection <- function(t, mu, kappa, sigma, V0 = NULL,
                                    stationary_init = TRUE) {

  t <- as.numeric(t)
  nt <- length(t)

  # Allow same parameter for x and y
  if (length(kappa) == 1) {
    kappa <- rep(kappa, 2)
  }

  if (length(sigma) == 1) {
    sigma <- rep(sigma, 2)
  }

  V <- matrix(NA_real_, nrow = nt, ncol = 2)
  colnames(V) <- c("Vx", "Vy")

  # Initial value v0
  if (!is.null(V0)) {
    V[1, ] <- V0
  } else if (stationary_init) {
    # stationary variance:
    # Var(V) = sigma^2 / (2*kappa)
    stationary_sd <- sigma / sqrt(2 * kappa)
    V[1, ] <- mu + stationary_sd * rnorm(2)

  } else {
    # Start exactly at episode mean
    V[1, ] <- mu
  }


  # Exact OU transitions
  if (nt >= 2) {
    for (j in 2:nt) {
      # Time increment
      dt <- t[j] - t[j - 1]
      # Phi = exp(-kappa * dt) is the autocorrelation at lag dt
      phi <- exp(-kappa * dt)
      
      # Standard deviation of the innovation term
      # Var(V_j | V_{j-1}) = sigma^2 * (1 - exp(-2 * kappa * dt)) / (2 * kappa)
      innovation_sd <- sigma * sqrt(
                        (1 - exp(-2 * kappa * dt)) /
                            (2 * kappa))
      V[j, ] <- mu + phi * (V[j - 1, ] - mu) +
                    innovation_sd * rnorm(2)
    }
  }
  return(V)
}



# Integrate the advection trajectory
# D(t_j) = integral from t[1] to t_j of V(u) du
integrate_advection <- function(V, t) {
  V <- as.matrix(V)
  t <- as.numeric(t)
  nt <- length(t)

  D <- matrix(0, nrow = nt, ncol = 2)
  colnames(D) <- c("Dx", "Dy")

  if (nt >= 2) {
    for (j in 2:nt) {
      dt <- t[j] - t[j - 1]
      # integral V(u)du
      # ~= dt/2 * (V_j-1 + V_j)
      D[j, ] <- D[j - 1, ] + 0.5 * dt * (V[j - 1, ] + V[j, ])
    }
  }
  return(D)
}



# make unique spatial shifted coordinates
unique_shifted_coordinates <- function(x, y) {
  key <- paste(sprintf("%.15g", x), sprintf("%.15g", y), sep = "_")
  unique_index <- !duplicated(key)
  list(
    x = x[unique_index],
    y = y[unique_index],
    map = match(key, key[unique_index])
  )
}



# r-Pareto simulator with OU advection
sim_rpareto_OU <- function(beta1, beta2, alpha1, alpha2, x, y, t,
            adv_mean = c(0, 0), kappa = 1, sigma = 0.1, V0 = NULL, 
            stationary_init = TRUE, s0 = c(1, 1), t0 = t[1],
            nres = 1, distance = "euclidean", seed = NULL,
            return_latent = TRUE) {

  # Init
  if (!is.null(seed)) {
    set.seed(seed)
  }
  distance <- match.arg(distance)
  RandomFields::RFoptions(spConform = FALSE,allow_duplicated_locations = TRUE, 
                            install = "no")
  x <- as.numeric(x)
  y <- as.numeric(y)
  t <- as.numeric(t)
  lx <- length(x)
  ly <- length(y)
  lt <- length(t)

  # Gaussian models for variograms
  #
  # gamma_s(h) = beta1 * ||h||^alpha1
  # gamma_t(tau) = beta2 * |tau|^alpha2
  #
  # RandomFields RMfbm uses var = 2 * beta
  # because RFvariogram returns the semivariogram.
  modelSpace <- RandomFields::RMfbm(
    alpha = alpha1,
    var = 2 * beta1
  )
  modelTime <- RandomFields::RMfbm(
    alpha = alpha2,
    var = 2 * beta2
  )


  # Original regular grid
  grid <- expand.grid(
    x = seq_len(lx),
    y = seq_len(ly),
    it = seq_len(lt)
  )
  grid$t <- t[grid$it]
  n_points <- nrow(grid)


  # Conditioning index
  t0_index <- match(t0, t)
  ind0 <- which(grid$x == s0[1] & grid$y == s0[2] &
                grid$it == t0_index)

  if (length(ind0) != 1) {
    stop("Conditioning point could not be identified uniquely.")
  }

  # final r-Pareto
  Z <- array(NA_real_, dim = c(lx, ly, lt, nres))
  # advection
  V_store <- array(NA_real_, dim = c(lt, 2, nres))
  # displacement
  D_store <- array(NA_real_, dim = c(lt, 2, nres))
  # variogram
  gamma_store <- array(NA_real_, dim = c(lx, ly, lt, nres))
  # gaussian field
  W_store <- array(NA_real_, dim = c(lx, ly, lt, nres))


  # Loop over r-Pareto simulations
  for (isim in seq_len(nres)) {

    # Simulate one OU advection trajectory
    V_t <- simulate_OU_advection(
      t = t,
      mu = adv_mean,
      kappa = kappa,
      sigma = sigma,
      V0 = V0,
      stationary_init = stationary_init
    )


    # Cumulative advective displacement
    # D(t) = int V(u) du
    D_t <- integrate_advection(V = V_t, t = t)

    V_store[, , isim] <- V_t
    D_store[, , isim] <- D_t


    # Advected coordinates
    #
    # Original constant-advection model:
    #     x* = x - t * V
    #     y* = y - t * V
    # becomes
    #     x*(t) = x - D_x(t)
    #     y*(t) = y - D_y(t)

    grid$shifted_x <- x[grid$x] - D_t[grid$it, 1]
    grid$shifted_y <- y[grid$y] - D_t[grid$it, 2]


    # Coordinates of conditioning point in transported frame
    shifted_x0 <- grid$shifted_x[ind0]
    shifted_y0 <- grid$shifted_y[ind0]


    # Variogram relative to (s0,t0)
    # gamma(s,t ; s0,t0) = beta1 ||[s - D(t)] - [s0 - D(t0)]||^alpha1
    #                      + beta2 |t-t0|^alpha2
    if (distance == "euclidean") {
      gamma_space_vec <-
        RandomFields::RFvariogram(
          modelSpace,
          x = grid$shifted_x - shifted_x0,
          y = grid$shifted_y - shifted_y0
        )

    } else {
      # L_alpha-type additive spatial model
      gamma_space_x <-
        RandomFields::RFvariogram(
          modelSpace,
          x = grid$shifted_x - shifted_x0
        )

      gamma_space_y <-
        RandomFields::RFvariogram(
          modelSpace,
          x = grid$shifted_y - shifted_y0
        )

      gamma_space_vec <- gamma_space_x + gamma_space_y
    }


    gamma_temp_vec <-
      RandomFields::RFvariogram(
        modelTime,
        x = grid$t - t0
      )


    gamma_vec <- as.numeric(gamma_space_vec) +
                 as.numeric(gamma_temp_vec)

    gamma_array <- array(gamma_vec, dim = c(lx, ly, lt))


    # Simulate spatial Gaussian process
    # on the ADVECTED coordinates
    # Duplicate coordinates are simulated only once.
    unique_coords <- unique_shifted_coordinates(
                                        grid$shifted_x,
                                        grid$shifted_y
                                    )

    if (distance == "euclidean") {
      W_s_unique <-
        RandomFields::RFsimulate(
          modelSpace,
          x = unique_coords$x,
          y = unique_coords$y,
          grid = FALSE
        )
      W_s_unique <- as.numeric(W_s_unique)
      W_s_vec <- W_s_unique[unique_coords$map]

    } else {
        # error message for L_alpha-type additive spatial model
        stop("L_alpha-type additive spatial model not implemented yet.")
    }

    # Temporal Gaussian process
    W_t <-
      RandomFields::RFsimulate(
        modelTime,
        x = t,
        grid = TRUE
      )
    W_t <- as.numeric(W_t)


    # Full spatio-temporal Gaussian field
    # W(s,t) = Ws(s - D(t)) + Wt(t)
    W_vec <- W_s_vec + W_t[grid$it]
    W_array <- array(W_vec, dim = c(lx, ly, lt))



    # Brown-Resnick / log-Gaussian spectral process
    W0 <- W_vec[ind0]
    Y_vec <- exp(W_vec - W0 - gamma_vec)

    # Pareto variable
    R <- 1 / runif(1)

    # r-Pareto field
    Z_vec <- R * Y_vec
    Z[, , , isim] <- array(Z_vec, dim = c(lx, ly, lt))
    gamma_store[, , , isim] <- gamma_array
    W_store[, , , isim] <- W_array
  }


  # Return
  out <- list(
    Z = Z,
    V = V_store, # OU wind
    D = D_store, # integrated displacement
    # conditioning information
    s0 = s0,
    t0 = t0,
    # parameters
    adv_mean = adv_mean,
    kappa = kappa,
    sigma = sigma,
    t = t,
    x = x,
    y = y
  )
  if (return_latent) {
    out$gamma <- gamma_store
    out$W <- W_store
  }
  return(out)
}


########################################################################


library(RandomFields)
library(RandomFieldsUtils)

x <- 1:20
y <- 1:20

# 5 minutes = 1/12 heure
# t <- seq(from = 0, to = 2, by = 5 / 60)
t <- seq(from = 0, to = 10, by = 1)

res <- sim_rpareto_OU(
  beta1 = 0.4, beta2 = 0.2, alpha1 = 1.5, alpha2 = 1.0,
  x = x, y = y, t = t,
  # mean wind/advection of the episode
  adv_mean = c(1, 0.4),
  # correlation time = 1/kappa
  kappa =  0.5,
  sigma = c(1, 0.4),
  s0 = c(1, 1),
  t0 = 0,
  nres = 1,
  seed = 123
)

plot(t, res$V[, 1, 1], type = "l", xlab = "Time", ylab = "Vx")
abline(h = res$adv_mean[1], lty = 2)

plot(t, res$V[, 2, 1], type = "l", xlab = "Time", ylab = "Vy")
abline(h = res$adv_mean[2], lty = 2)

plot(t, res$Z[1, 1, , 1],
  type = "l", xlab = "Time",ylab = "Z(1,1,t)",
  main = "r-Pareto process at (1,1)")


# constant
res0 <- sim_rpareto_OU(
  beta1 = 0.4,
  beta2 = 0.2,
  alpha1 = 1.5,
  alpha2 = 1,
  x = x,
  y = y,
  t = t,
  adv_mean = c(1.4, 0.4),
  kappa = 1,
  sigma = 0,
  s0 = c(1,1),
  t0 = 0,
  nres = 1,
  seed = 123
)

res1 <- sim_rpareto_OU(
  beta1 = 0.4,
  beta2 = 0.2,
  alpha1 = 1.5,
  alpha2 = 1,
  x = x,
  y = y,
  t = t,
  adv_mean = c(1.4, 0.4),
  kappa = 1,
  sigma = 0.4,
  s0 = c(1,1),
  t0 = 0,
  nres = 1,
  seed = 123
)

res2 <- sim_rpareto_OU(
  beta1 = 0.4,
  beta2 = 0.2,
  alpha1 = 1.5,
  alpha2 = 1,
  x = x,
  y = y,
  t = t,
  adv_mean = c(5, 4),
  kappa = 1,
  sigma = 0.4,
  s0 = c(1,1),
  t0 = 0,
  nres = 1,
  seed = 123
)

# image(res1$Z[,,1,1])
# image(res0$Z[,,1,1])

# compare res0 and res1 for one site in time
plot(t, res0$Z[1, 1, , 1], type = "l", xlab = "Time", ylab = "Z",
     main = "r-Pareto process at (1,1)", ylim = range(c(res0$Z[1, 1, , 1], res1$Z[1, 1, , 1])))
lines(t, res1$Z[1, 2, , 1], col = "red")
legend("topright", legend = c("constant advection", "OU advection"),
       col = c("black", "red"), lty = 1)

ymax <- max(res0$Z, res1$Z, na.rm = TRUE)
# plot every time series of res0 on the same plot
plot(t, res0$Z[1, 1, , 1], type = "l", xlab = "Time", ylab = "Z",
     main = "r-Pareto process", ylim = c(0, ymax))
for (i in 2:20) {
  for (j in 1:20) {
    lines(t, res0$Z[i, j, , 1], col = "black")
  }
}
for (i in 1:10) {
    for (j in 1:10) {
        lines(t, res1$Z[i, j, , 1], col = "red")
    }
}
legend("topright", legend = c("constant advection", "OU advection"),
       col = c("black", "red"), lty = 1)



ymax <- max(res2$Z, res1$Z, na.rm = TRUE)
# plot every time series  on the same plot
plot(t, res2$Z[1, 1, , 1], type = "l", xlab = "Time", ylab = "Z",
     main = "r-Pareto process", ylim = c(0, ymax))
for (i in 2:10) {
  for (j in 1:10) {
    lines(t, res2$Z[i, j, , 1], col = "black")
  }
}
for (i in 1:10) {
    for (j in 1:10) {
        lines(t, res1$Z[i, j, , 1], col = "red")
    }
}
