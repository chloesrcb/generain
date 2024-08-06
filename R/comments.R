

# empirical_excesses <- function(data_rain, quantile, h_vect,
#                                 nmin = 5) {
#   Tmax <- nrow(data_rain) # number of time steps
#   q <- quantile # quantile

#   unique_tau <- unique(h_vect$tau)

#   h_vect$N_vect <- NA
#   h_vect$n_vect <- NA
#   for (t in unique_tau) {
#     # df_dist_t <- df_dist[df_dist$tau == t, ]
#     df_h_t <- h_vect[h_vect$tau == t, ]
#     for (i in seq_len(nrow(df_h_t))) {
#       # get index pairs
#       ind_s2 <- as.numeric(as.character(df_h_t$s2[i]))
#       ind_s1 <- df_h_t$s1[i]
#       # get the couple of sites
#       rain_cp <- drop_na(data_rain[, c(ind_s1, ind_s2)])
#       colnames(rain_cp) <- c("s1", "s2")
#       rains1 <- rain_cp$s1
#       rains2 <- rain_cp$s2

#       rain_nolag <- rains1[1:(Tmax - t)] # without lag (in t_k)
#       rain_lag <- rains2[(1 + t):Tmax] # with lag t (in t_k + t)
#       data_cp <- cbind(rain_nolag, rain_lag) # get final couple
#       n <- nrow(data_cp)

#       if (n >= nmin) {
#         rain_unif <- cbind(rank(data_cp[, 1]) / (n + 1),
#                           rank(data_cp[, 2]) / (n + 1))

#         # check excess above a threshold q
#         cp_cond <- rain_unif[rain_unif[, 2] > q, ]

#         if (length(class(cp_cond)) == 1 && class(cp_cond) == "numeric") {
#           # if only one excess
#           cp_cond <- t(as.matrix(cp_cond))
#         }

#         # nb of conditional excesses
#         excess_count <- sum(cp_cond[, 1] > q)
#         # n_vect <- c(n_vect, excess_count)
#         # N_vect <- c(N_vect, nrow(cp_cond))
#         h_vect[h_vect$s1 == ind_s1 & h_vect$s2 == ind_s2 & h_vect$tau == t, ]$n_vect <- excess_count
#         h_vect[h_vect$s1 == ind_s1 & h_vect$s2 == ind_s2 & h_vect$tau == t, ]$N_vect <- nrow(cp_cond)
#       }
#     }
#   }
#   return(h_vect)
# }


# theorical_chi_ind <- function(params, h, tau, adv) {
#   # get variogram parameter
#   beta1 <- params[1]
#   beta2 <- params[2]
#   # beta3 <- params[3]
#   alpha1 <- params[3]
#   alpha2 <- params[4]
#   # alpha3 <- params[6]
#   h1_adv <- abs(h1 + adv[1] * tau)
#   h2_adv <- abs(h2 + adv[2] * tau)
#   hnorm <- norm_Lp(h1_adv, h2_adv, p = alpha1)
#   # Get vario and chi for each lagtemp
#   varioval <- 2 * (beta1 * hnorm^alpha1 + beta2 * tau^alpha2)
#   phi <- pnorm(sqrt(0.5 * varioval))
#   chival <- 2 * (1 - phi)

#   return(chival)
# }

# evaluate_optim <- function(list_simu, quantile, true_param, tau, hmax,
#                            locations, nmin = 5,
#                            parscale = c(1, 1, 1, 1), latlon = FALSE) {

#   lower.bound <- c(1e-6, 1e-6, 1e-6, 1e-6)
#   upper.bound <- c(Inf, Inf, 1.999, 1.999)
#   if (length(true_param) == 6) {
#     lower.bound <- c(lower.bound, -1e-6, -1e-6)
#     upper.bound <- c(upper.bound, Inf, Inf)
#     parscale <- c(parscale, 1, 1)
#   }
#   # get the number of simulations
#   n_res <- length(list_simu)
#   # create a dataframe to store the results
#   df_result <- data.frame(beta1 = rep(NA, n_res), beta2 = rep(NA, n_res),
#                           alpha1 = rep(NA, n_res), alpha2 = rep(NA, n_res))

#   # if there is advection
#   if (length(true_param) == 6) {
#     df_result$adv1 <- rep(NA, n_res)
#     df_result$adv2 <- rep(NA, n_res)
#   }

#   h_vect <- get_lag_vectors(locations, true_param, tau = tau, hmax = hmax)
#   count_cv <- 0
#   # for all simulations
#   for (n in 1:n_res) {
#     simu_df <- as.data.frame(list_simu[[n]]) # get the simulation dataframe
#     # get the empirical excesses
#     if (length(true_param) == 6) {
#       excesses <- NULL
#     } else {
#       excesses <- empirical_excesses(simu_df, quantile, tau, h_vect,
#                                    nmin)
#     }
#     # optimize the negative log-likelihood function
#     tryCatch({
#         # result <- optimr(par = true_param, method = "Rcgmin",
#         #           gr = "grfwd", fn = function(par) {
#         #           neg_ll(par, simu = simu_df, quantile = quantile,
#         #                 h_vect = h_vect, tau = tau,
#         #                 locations = locations, latlon = latlon,
#         #                 nmin = nmin, excesses = NULL)
#         #           }, lower = lower.bound, upper = upper.bound,
#         #           control = list(parscale = parscale,
#         #                          maxit = 10000))


#           result <- optim(par = true_param, fn = neg_ll,
#                       excesses = excesses, quantile = quantile,
#                       h_vect = h_vect, tau = tau,
#                       locations = locations, simu = simu_df,
#                       method = "CG",
#                       control = list(parscale = parscale,
#                                     maxit = 10000))

#         print(result$convergence)
#         if (result$convergence == 0) { # if it converges
#           count_cv <- count_cv + 1
#           params <- result$par
#           df_result$beta1[n] <- params[1]
#           df_result$beta2[n] <- params[2]
#           df_result$alpha1[n] <- params[3]
#           df_result$alpha2[n] <- params[4]
#           if (length(true_param) == 6) {
#               df_result$adv1[n] <- params[5]
#               df_result$adv2[n] <- params[6]
#           }
#         } else {
#           print(n)
#         }
#     }, error = function(e) {
#         # Handle the error (e.g., print an error message)
#         print(paste("Error occurred for simulation", n))
#     })
#   }
#   print(paste0("Number of convergence: ", count_cv))
#   return(df_result)
# }




# #' Compute the distance matrix between locations
# #'
# #' This function computes the distance matrix between a set of locations.
# #'
# #' @param locations A matrix or data frame containing the coordinates of the
# #' locations.
# #' @param dmax The maximum distance threshold if specified.
# #' @param latlon Logical indicating whether the coordinates are in latitude and
# #' longitude format.
# #' @return A distance matrix where each element represents the distance between
# #' two locations.
# #'
# #' @import geodist
# #' @importFrom stats dist
# #'
# #' @export
# get_dist_mat <- function(locations, dmax = NA, latlon = TRUE) {
#   # get longitude and latitude in a dataframe to get distance between points
#   loc <- data.frame(lat = locations$Latitude, lon = locations$Longitude)
#   # get distance matrix
#   if(latlon) {
#     dist_mat <- geodist(loc, measure = "haversine") # meters by default
#   } else {
#     dist_mat <- as.matrix(dist(loc))
#   }
#   if (!is.na(dmax)) {
#     dist_mat[dist_mat > dmax] <- 0
#   }
#   return(dist_mat)
# }



# reshape_distances <- function(dist_mat) {
#   # convert distance matrix into a dataframe
#   tmax <- length(names(dist_mat))
#   if (tmax == 0) {
#     df_dist <- as.data.frame(dist_mat)
#     n <- nrow(df_dist)
#     colnames(df_dist) <- c(1:n)
#     rownames(df_dist) <- c(1:n)

#     # Make a triangle
#     df_dist[lower.tri(df_dist)] <- NA

#     # Convert to a data frame, and add tenure labels
#     df_dist <- as.data.frame(df_dist)
#     df_dist$Y <- 1:n

#     # Reshape to suit ggplot, remove NAs, and sort the labels
#     df_dist <- na.omit(reshape2::melt(df_dist, "Y", variable_name = "X"))
#     colnames(df_dist) <- c("Y", "X", "value")
#     df_dist$X <- factor(df_dist$X, levels = rev(levels(df_dist$X)))
#   } else { # with advection
#     df_dist <- data.frame()
#     for (i in 1:tmax) {
#       df_dist_t <- as.data.frame(dist_mat[i])
#       n <- nrow(df_dist_t)
#       colnames(df_dist_t) <- c(1:n)
#       rownames(df_dist_t) <- c(1:n)

#       # Make a triangle
#       df_dist_t[lower.tri(df_dist_t)] <- NA

#       # Convert to a data frame, and add tenure labels
#       df_dist_t <- as.data.frame(df_dist_t)
#       df_dist_t$Y <- 1:n

#       # Reshape to suit ggplot, remove NAs, and sort the labels
#       df_dist_t <- na.omit(reshape2::melt(df_dist_t, "Y", variable_name = "X"))
#       colnames(df_dist_t) <- c("Y", "X", "value")
#       df_dist_t$X <- factor(df_dist_t$X, levels = rev(levels(df_dist_t$X)))
#       df_dist_t$tau <- i
#       df_dist <- rbind(df_dist, df_dist_t)
#       # df_dist[[paste0("t", t)]] <- df_dist_t
#     }
#   }
#   return(df_dist)
# }


# get_lag_vectors <- function(df_coords, params, hmax = NA, tau_vect = 1:10) {
#   alpha_spa <- params[3]
#   if (length(params) != 6) {
#     adv <- c(0, 0)
#   } else {
#     adv <- params[5:6]
#   }

#   n <- nrow(df_coords)
#   lags <- data.frame()

#   for (i in 1:(n - 1)) { # Loop over all pairs of points
#     for (j in (i + 1):n) {
#       for (tau in tau_vect) { # Loop over tau values
#         # Calculate lag vector without advection
#         lag_latitude <- df_coords$Latitude[j] - df_coords$Latitude[i]
#         lag_longitude <- df_coords$Longitude[j] - df_coords$Longitude[i]

#         # Calculate lag vector with advection
#         adv_latitude <- lag_latitude - adv[2] * tau
#         adv_longitude <- lag_longitude - adv[1] * tau
#         # hnorm <- sqrt(adv_latitude^2 + adv_longitude^2) # Distance
#         hnorm <- norm_Lp(adv_latitude, adv_longitude, alpha_spa)

#         # Add results to a new row in the dataframe
#         if (hnorm <= hmax) {
#           new_row <- data.frame(
#             s1 = i,
#             s2 = j,
#             h1 = lag_latitude,
#             h2 = lag_longitude,
#             h1_adv = adv_latitude,
#             h2_adv = adv_longitude,
#             tau = tau,
#             hnorm = hnorm
#           )
#           lags <- bind_rows(lags, new_row)
#         }
#       }
#     }
#   }
#   return(lags)
# }



# sim_rpareto <- function(beta1, beta2, alpha1, alpha2, x, y, t, n.res,
#                         adv = c(0, 0)) {
#   # beta1, beta2, alpha1, alpha2 are variogram parameters
#   # x is the first dimension (spatial x in our case)
#   # y is the second dimension (spatial y in our case)
#   # z is the third dimension (time in our case)
#   # (adv1, adv2) advection coordinates vector
#   ## Setup
#   RandomFields::RFoptions(spConform = FALSE, install = "no")
#   lx <- length(sx <- seq_along(x))  # spatial
#   ly <- length(sy <- seq_along(y))  # spatial
#   lt <- length(st <- seq_along(t))  # temporal
  
#   ## Model-Variogram BuhlCklu
#   modelBuhlCklu <- RandomFields::RMfbm(alpha = alpha1, var = beta1, proj = 1) +
#                    RandomFields::RMfbm(alpha = alpha1, var = beta1, proj = 2) +
#                    RandomFields::RMfbm(alpha = alpha2, var = beta2, proj = 3)

#   ## Construct grid
#   Nxy <- lx * ly # spatial grid size
#   N <- Nxy * lt # spatio-temporal grid size
#   grid <- expand.grid(x = x, y = y, t = t)

#   ## Construct shifted variogram
#     if (all(adv == c(0, 0))) {

#     Varm1 <- vapply(seq_len(N), function(n) {
#       dx <- sx - grid[n, 1] # spatial lags
#       dy <- sy - grid[n, 2] # spatial lags
#       dt <- st - grid[n, 3] # temporal lags
#       result <- RandomFields::RFvariogram(modelBuhlCklu,
#                   x = dx,
#                   y = dy,
#                   z = dt
#                 )
#       return(result)
#     }, array(NA_real_, dim = c(lx, ly, lt))) ## => (lx, ly, lt, N)-array
#   } else {
#     ## Construct shifted variogram
#     Varm1 <- vapply(seq_len(N), function(n) {
#       dx <- sx - grid[n, 1] # spatial lags
#       dy <- sy - grid[n, 2] # spatial lags
#       dt <- st - grid[n, 3] # temporal lags
#       # norm_h <- sqrt(combi$dx_adv^2 + combi$dy_adv^2)
#       # compute variogram for each combination
#       combi <- expand.grid(dx = dx, dy = dy, dt = dt) # combinations of lags
#       combi$dx_adv <- combi$dx - adv[1] * combi$dt # adding advection on x
#       combi$dy_adv <- combi$dy - adv[2] * combi$dt # adding advection on y
#       result <- RandomFields::RFvariogram(modelBuhlCklu,
#                   x = combi$dx_adv,
#                   y = combi$dy_adv,
#                   z = combi$dt
#                 )
#       return(result)
#     }, array(NA_real_, dim = c(lx, ly, lt))) ## => (lx, ly, lt, N)-array

#   }

#   # Main
#   Z <- array(, dim = c(lx, ly, lt, n.res)) # 3d array
#   for (i in seq_len(n.res)) {
#     W <- RandomFields::RFsimulate(modelBuhlCklu, x, y, t) # GP
#     Y <- exp(W - W[1] - Varm1[,,, 1])
#     R <- evd::rgpd(n = 1, loc = 1, scale = 1, shape = 1)
#     Z[,,, i] <- R * Y
#   }
#   # Return
#   Z
# }


# sim_BR <- function(beta1, beta2, alpha1, alpha2, x, y, t, n.BR, adv = c(0, 0)) {
#   # beta1, beta2, alpha1, alpha2 are variogram parameters
#   # x is the first dimension (spatial x in our case)
#   # y is the second dimension (spatial y in our case)
#   # z is the third dimension (time in our case)
#   # (adv1, adv2) advection coordinates vector
#   ## Setup
#   RandomFields::RFoptions(spConform = FALSE, install = "no")
#   lx <- length(sx <- seq_along(x))  # spatial
#   ly <- length(sy <- seq_along(y))  # spatial
#   lt <- length(st <- seq_along(t))  # temporal

#   ## Model-Variogram BuhlCklu (fractional Brownian motion)
#   modelBuhlCklu <- RandomFields::RMfbm(alpha = alpha1, var = beta1, proj = 1) +
#                    RandomFields::RMfbm(alpha = alpha1, var = beta1, proj = 2) +
#                    RandomFields::RMfbm(alpha = alpha2, var = beta2, proj = 3)

#   ## Construct grid
#   Nxy <- lx * ly # spatial grid size
#   N <- Nxy * lt # spatio-temporal grid size
#   grid <- matrix(0, nrow=N, ncol=3) # (N,3)-matrix

#   for (i in sx)
#     for (j in seq_len(ly*lt))
#       grid[i+(j-1)*ly, 1] <- i

#   for (i in sy)
#     for (j in sx)
#       for(k in st)
#         grid[j+lx*(i-1)+(k-1)*Nxy, 2] <- i

#   for (i in st)
#     for (j in seq_len(Nxy))
#       grid[j+Nxy*(i-1), 3] <- i

#   if (all(adv == c(0, 0))) {

#     Varm1 <- vapply(seq_len(N), function(n) {
#       dx <- sx - grid[n, 1] # spatial lags
#       dy <- sy - grid[n, 2] # spatial lags
#       dt <- st - grid[n, 3] # temporal lags
#       result <- RandomFields::RFvariogram(modelBuhlCklu,
#                   x = dx,
#                   y = dy,
#                   z = dt
#                 )
#       return(result)
#     }, array(NA_real_, dim = c(lx, ly, lt)))
#   } else {
#     ## Construct shifted variogram
#     Varm1 <- vapply(seq_len(N), function(n) {
#       dx <- sx - grid[n, 1] # spatial lags
#       dy <- sy - grid[n, 2] # spatial lags
#       dt <- st - grid[n, 3] # temporal lags
#       # norm_h <- sqrt(combi$dx_adv^2 + combi$dy_adv^2)
#       # compute variogram for each combination
#       combi <- expand.grid(dx = dx, dy = dy, dt = dt) # combinations of lags
#       combi$dx_adv <- combi$dx - adv[1] * combi$dt # adding advection on x
#       combi$dy_adv <- combi$dy - adv[2] * combi$dt # adding advection on y
#       result <- RandomFields::RFvariogram(modelBuhlCklu,
#                   x = combi$dx_adv,
#                   y = combi$dy_adv,
#                   z = combi$dt
#                 )
#       return(result)
#     }, array(NA_real_, dim = c(lx, ly, lt)))

#   }

#   ## Main
#   Z <- array(, dim = c(lx, ly, lt, n.BR)) # 3d array
#   E <- matrix(rexp(n.BR * N), nrow = n.BR, ncol = N)

#   for (i in seq_len(n.BR)) {
#     # for n = 1
#     V <- 1 / E[i, 1] # poisson process
#     W <- RandomFields::RFsimulate(modelBuhlCklu, x, y, t) # GP
#     Y <- exp(W - W[1] - Varm1[,,, 1])
#     Z[,,,i] <- V * Y

#     # n in {2,..,N}
#     for (n in 2:N) {
#       Exp <- E[i, n]
#       V <- 1 / Exp
#       while(V > Z[N * (i - 1) + n]) {
#         W <- RandomFields::RFsimulate(modelBuhlCklu, x, y, t)
#         Y <- exp(W - W[n] - Varm1[,,, n])
#         if (all(V * Y[seq_len(n - 1)] < Z[(N * (i - 1) + 1):(N * (i - 1) + (n - 1))])) {
#           Z[,,, i] <- pmax(V * Y, Z[,,, i])
#         }
#         Exp <- Exp + rexp(1)
#         V <- 1 / Exp
#       }
#     }
#   }
#   ## Return
#   Z
# }


