# EXCESSES ---------------------------------------------------------------------

#' empirical_excesses function
#'
#' This function calculates the empirical excesses based on indicators above a
#' quantile threshold.
#'
#' @param data_rain The rainfall data.
#' @param quantile The quantile value.
#' @param tau The temporal lag vector.
#' @param h_vect The spatial lag vectors and distances dataframe.
#' @param nmin The minimum number of observations to consider. Default is 5.
#'
#' @return A list with n_vect the number of excesses and N_vect the number of
#' possible excesses ie the number of observations for each pair of sites.
#'
#' @import tidyr
#'
#' @export
empirical_excesses <- function(data_rain, quantile, tau, h_vect,
                                nmin = 5) {
  Tmax <- nrow(data_rain) # number of time steps
  q <- quantile # quantile
  # h_norms <- unique(h_vectors$hnorm)
  # df_dist$h <- ifelse(df_dist$value %in% h_norms, df_dist$value, NA)
  # N_vect <- c() # number of observations
  # n_vect <- c() # number of excesses (sum)
  # if it is an only spatial matrix
  # df_dist$tau <- ifelse(is.null(df_dist$tau), 0, df_dist$tau)
  tau_vect <- h_vect$tau
  h_vect$N_vect <- NA
  h_vect$n_vect <- NA
  for (t in unique(tau_vect)) {
    # df_dist_t <- df_dist[df_dist$tau == t, ]
    df_h_t <- h_vect[h_vect$tau == t, ]
    for (i in seq_len(nrow(df_h_t))) {
      # get index pairs
      ind_s2 <- as.numeric(as.character(df_h_t$s2[i]))
      ind_s1 <- df_h_t$s1[i]
      # get the couple of sites
      rain_cp <- drop_na(data_rain[, c(ind_s1, ind_s2)])
      colnames(rain_cp) <- c("s1", "s2")
      rains1 <- rain_cp$s1
      rains2 <- rain_cp$s2

      rain_nolag <- rains1[1:(Tmax - t)] # without lag (in t_k)
      rain_lag <- rains2[(1 + t):Tmax] # with lag t (in t_k + t)
      data_cp <- cbind(rain_nolag, rain_lag) # get final couple
      n <- nrow(data_cp)

      if (n >= nmin) {
        rain_unif <- cbind(rank(data_cp[, 1]) / (n + 1),
                          rank(data_cp[, 2]) / (n + 1))

        # check excess above a threshold q
        cp_cond <- rain_unif[rain_unif[, 2] > q, ]

        if (length(class(cp_cond)) == 1 && class(cp_cond) == "numeric") {
          # if only one excess
          cp_cond <- t(as.matrix(cp_cond))
        }

        # nb of conditional excesses
        excess_count <- sum(cp_cond[, 1] > q)
        # n_vect <- c(n_vect, excess_count)
        # N_vect <- c(N_vect, nrow(cp_cond))
        h_vect[h_vect$s1 == ind_s1 & h_vect$s2 == ind_s2 & h_vect$tau == t, ]$n_vect <- excess_count
        h_vect[h_vect$s1 == ind_s1 & h_vect$s2 == ind_s2 & h_vect$tau == t, ]$N_vect <- nrow(cp_cond)
      }
    }
  }
  return(h_vect)
}

# THEORICAL CHI ----------------------------------------------------------------

#' Calculate the theoretical chi-indicator value.
#'
#' This function calculates the individual theoretical spatio-temporal chi value
#' based on a the variogram parameters, a spatial lag and a temporal lag.
#'
#' @param params The variogram parameter values (beta1, beta2, alpha1, alpha2).
#' @param h The spatial lag value.
#' @param tau The temporal lag value.
#'
#' @return The theoretical chi value.
#'
#' @import stats
#'
#' @export
theorical_chi_ind <- function(params, h, tau) {
  # get variogram parameter
  beta1 <- params[1]
  beta2 <- params[2]
  alpha1 <- params[3]
  alpha2 <- params[4]

  # Get vario and chi for each lagtemp
  varioval <- 2 * (beta1 * h^alpha1 + beta2 * tau^alpha2)
  phi <- pnorm(sqrt(0.5 * varioval))
  chival <- 2 * (1 - phi)

  return(chival)
}

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


#' Compute the theoretical chi matrix.
#'
#' This function calculates the theoretical chi matrix based on the given
#' parameters.
#'
#' @param params A vector of parameters.
#' @param h_vect A vector of spatial lag values.
#' @param tau A vector of temporal lag values.
#'
#' @return The theoretical chi matrix.
#'
#' @export
theorical_chi_mat <- function(params, h_vect, tau) {
  # Init matrix
  chi <- matrix(0, nrow = length(tau), ncol = length(h_vect))
  adv <- params[5:6]
  # for each spatial lag
  for (j in seq_along(h_vect)) {
      # Get vario and chi for each lagtemp
      h <- h_vect[j]
      for (t in seq_along(tau)) {
          chi[t, j] <- theorical_chi_ind(params, h, tau[t], adv)
      }
  }
  return(chi)
}


#' Compute the theoretical chi dataframe.
#'
#' This function calculates the theoretical chi dataframe. based on the given
#' parameters.
#'
#' @param params A vector of parameters.
#' @param h_vectors A dataframe with spatial lag values.
#' @param tau A vector of temporal lag values.
#'
#' @return The theoretical chi matrix.
#'
#' @export
theorical_chi <- function(params, h_vectors, tau) {
  # Init data frame
  lx <- nrow(h_vectors)
  lt <- length(tau)

  chi_df <- data.frame(
    hnorm = rep(h_vectors$hnorm, each = lt),
    tau = rep(tau, times = lx),
    chi = rep(0, times = lx * lt)
  )

  # for each spatial lag
  for (i in 1:nrow(h_vectors)) {
    h <- h_vectors$hnorm[i]
    tau <- h_vectors$tau[i]
    for (t in seq_along(tau)) {
      chi_value <- theorical_chi_ind(params, h, tau[t])
      chi_df$chi[(i - 1) * lt + t] <- chi_value
    }
  }

  return(chi_df)
}

#' get_chi_vect function
#'
#' This function calculates the theorical spatio-temporal extremogram vector
#' based on the given theorical spatio-temporal extremogram matrix.
#'
#' @param chi_mat The theorical spatio-temporal extremogram matrix.
#' @param h_vect The spatial lags vector.
#' @param tau The temporal lag vector.
#' @param df_dist The distances long dataframe.
#'
#' @return The theorical spatio-temporal extremogram vector.
#'
#' @export
get_chi_vect <- function(chi_mat, h_vect, tau, df_dist) {
  df_dist$h <- ifelse(df_dist$value %in% h_vect, df_dist$value, NA)
  # create excesses matrix
  npairs <- nrow(df_dist) # number of pairs
  nconfig <- npairs * length(tau) # number of configurations
  h_vect_all <- rep(df_dist$h, times = length(tau))
  tau_all <- rep(tau, each = npairs) # repeat tau for each pair
  chi_vect <- numeric(length = nconfig) # init chi vector
  for (p in 1:nconfig) {
    h <- h_vect_all[p] # get h for each configuration
    if (is.na(h)) { # if NA then chi = NA
      chi_vect[p] <- NA
    } else { # else get chi value
      t <- tau_all[p] # get tau for each configuration
      h_ind <- which(h_vect == h) # get index of h
      chi <- chi_mat[t, h_ind] # get chi value from matrix
      # chi inside ]0, 1] to avoid future error
      chi_vect[p] <- ifelse(chi <= 0, 0.000001, chi) # get chi value in vector
    }
  }
  return(chi_vect)
}

# NEGATIVE LOG-LIKELIHOOD ------------------------------------------------------

#' neg_ll function
#'
#' Calculate the negative log-likelihood for a given set of variogram
#' parameters by considering excess indicators following a binomial distribution
#' with a probability parameter equals to the theorical spatio-temporal chi
#' value. The number observations is different for each spatio-temporal
#' configuration.
#'
#' @param params Vector of variogram parameters (beta1, beta2, alpha1, alpha2).
#' @param excesses List of excesses (N_vect, n_vect).
#' @param h_vect The spatial lag vector.
#' @param tau The temporal lag vector.
#' @param df_dist The distances long dataframe.
#' @param locations The locations dataframe.
#' @param latlon A boolean value to indicate if the locations are in latitude
#'               and longitude. Default is FALSE.
#' @param simu_exp A boolean value to indicate if the data is simulated with
#'                 an exponential distribution.
#'
#' @return The negative log-likelihood value.
#'
#' @export
neg_ll <- function(params, simu, h_vect, tau, locations, # nolint
                  latlon = FALSE, quantile = 0.9, nmin = 5,
                  simu_exp = FALSE, excesses = NULL) {
  hmax <- max(h_vect$hnorm)
  # df_dist_new <- df_dist
  print(params)
  if (is.null(excesses)) {
    excesses <- empirical_excesses(simu, quantile, tau, h_vect,
                                  nmin)
  }

  lower.bound <- c(1e-6, 1e-6, 1e-6, 1e-6)
  upper.bound <- c(Inf, Inf, 1.999, 1.999)
  if (length(params) == 6) {
    lower.bound <- c(lower.bound, -1e-6, -1e-6)
    upper.bound <- c(upper.bound, Inf, Inf)
  }
  # Check if the parameters are in the bounds
  if (any(params < lower.bound) || any(params > upper.bound)) {
    print("out of bounds")
    return(1e8)
  }

  # h_vect_new <- h_vect

  if (length(params) == 6) {
    adv <- params[5:6]
    # change for each advection
    # dist_mat <- get_dist_mat(locations, adv = adv, tau = tau, latlon = latlon)
    # df_dist_new <- reshape_distances(dist_mat) # reshape the distance matrix
    # Create dataframe from h_adv
    # h_vect <- get_h_vect(df_dist_new, hmax)
    h_vect <- get_lag_vectors(locations, params, tau = tau, hmax = hmax)
    if (nrow(h_vect) <= 10) {
      return(1e8)
    }
    excesses <- empirical_excesses(simu, quantile, tau, h_vect, nmin = nmin)
  }

  N_vect <- excesses$N_vect # number of observations
  n_vect <- excesses$n_vect # number of excesses
  chi <- theorical_chi(params, h_vect, tau) # get chi matrix
  # transform in chi vector
  # chi_vect <- get_chi_vect(chi, h_vect$hnorm, tau)
  chi_vect <- as.vector(chi$chi)
  chi_vect <- ifelse(chi_vect <= 0, 0.000001, chi_vect)

  # logC <- lchoose(N_vect, n_vect) # log binomial coefficient
  non_excesses <- N_vect - n_vect # number of non-excesses
  # log-likelihood vector
  # ll_vect <- logC + n_vect * log(chi_vect) + non_excesses * log(1 - chi_vect)
  ll_vect <- n_vect * log(chi_vect) + non_excesses * log(1 - chi_vect)

  # negative log-likelihood
  nll <- -sum(ll_vect, na.rm = TRUE)
  return(nll)
}


# GRADIENT ---------------------------------------------------------------------

#' Compute the derivative of the extremogram with respect to beta1.
#'
#' This function calculates the derivative of the extremogram with respect to 
#' beta1.
#'
#' @param params a vector of variogram parameters.
#' @param hnorm a vector of spatial lags.
#' @param tau a vector of temporal lags.
#'
#' @return The derivative of the extremogram with respect to beta1.
#'
#' @examples
#' dchi_dbeta1(params = c(1, 2, 3), hnorm = 0.5, tau = 0.1)
#'
dchi_dbeta1 <- function(params, hnorm, tau) {
  # if adv is present it appaers in hnorm
  beta1 <- params[1]
  beta2 <- params[2]
  alpha1 <- params[3]
  alpha2 <- params[4]

  semivario <- beta1 * hnorm^alpha1 + beta2 * tau^alpha2
  phi <- pnorm(sqrt(semivario))
  dphi_dvarioval <- dnorm(sqrt(semivario)) * (1 / sqrt(semivario))

  dvarioval_dbeta1 <- 2 * hnorm^alpha1
  dchi_dbeta1 <- -2 * dphi_dvarioval * dvarioval_dbeta1
  return(dchi_dbeta1)
}

#' Compute the derivative of the extremogram function with respect to beta2.
#'
#' This function calculates the derivative of the extremogram function with
#' respect to beta2.
#'
#' @param params a vector of variogram parameters.
#' @param hnorm a vector of spatial lags.
#' @param tau a vector of temporal lags.
#'
#' @return The derivative of the extremogram function with respect to beta2.
#'
#' @examples
#' params <- c(1, 2)
#' hnorm <- c(0.1, 0.2, 0.3, 0.4)
#' tau <- 0.5
#' dchi_dbeta2(params, hnorm, tau)
#'
#' @export
dchi_dbeta2 <- function(params, hnorm, tau) {
  # if adv is present it appaers in hnorm
  beta1 <- params[1]
  beta2 <- params[2]
  alpha1 <- params[3]
  alpha2 <- params[4]

  semivario <- beta1 * hnorm^alpha1 + beta2 * tau^alpha2
  phi <- pnorm(sqrt(semivario))
  dphi_dvarioval <- dnorm(sqrt(semivario)) * (1 / sqrt(semivario))

  dvarioval_dbeta2 <- 2 * tau^alpha2
  dchi_dbeta2 <- -2 * dphi_dvarioval * dvarioval_dbeta2
  return(dchi_dbeta2)
}

#' Compute the derivative of the extremogram with respect to alpha1.
#'
#' This function calculates the derivative of theextremogram with respect to
#' alpha1.
#'
#' @param params a vector of variogram parameters.
#' @param hnorm a vector of spatial lags.
#' @param tau a vector of temporal lags.
#'
#' @return The derivative of the extremogram with respect to alpha1.
#'
#' @examples
#' params <- c(0.5, 0.8)
#' hnorm <- 0.2
#' tau <- 0.1
#' dchi_dalpha1(params, hnorm, tau)
#'
#' @export
dchi_dalpha1 <- function(params, hnorm, tau) {
  # if adv is present it appaers in hnorm
  beta1 <- params[1]
  beta2 <- params[2]
  alpha1 <- params[3]
  alpha2 <- params[4]

  semivario <- beta1 * hnorm^alpha1 + beta2 * tau^alpha2
  phi <- pnorm(sqrt(semivario))
  dphi_dvarioval <- dnorm(sqrt(semivario)) * (1 / sqrt(semivario))

  dvarioval_dalpha1 <- 2 * beta1 * hnorm^alpha1 * log(hnorm)
  dchi_dalpha1 <- -2 * dphi_dvarioval * dvarioval_dalpha1
  return(dchi_dalpha1)
}


#' Compute the derivative of the extremogram with respect to alpha2.
#'
#' This function calculates the derivative of the extremogram with respect
#' to alpha2.
#'
#' @param params a vector of variogram parameters.
#' @param hnorm a vector of spatial lags.
#' @param tau a vector of temporal lags.
#'
#' @return The derivative of the extremogram with respect to alpha2
#'
#' @examples
#' params <- c(0.5, 0.2)
#' hnorm <- c(0.1, 0.3, 0.2, 0.4)
#' tau <- 0.5
#' dchi_dalpha2(params, hnorm, tau)
#'
#' @export
dchi_dalpha2 <- function(params, hnorm, tau) {
  # if adv is present it appaers in hnorm
  beta1 <- params[1]
  beta2 <- params[2]
  alpha1 <- params[3]
  alpha2 <- params[4]

  semivario <- beta1 * hnorm^alpha1 + beta2 * tau^alpha2
  phi <- pnorm(sqrt(semivario))
  dphi_dvarioval <- dnorm(sqrt(semivario)) * (1 / sqrt(semivario))

  dvarioval_dalpha2 <- 2 * beta2 * tau^alpha2 * log(tau)
  dchi_dalpha2 <- -2 * dphi_dvarioval * dvarioval_dalpha2
  return(dchi_dalpha2)
}


#' Compute the gradient of the negative log-likelihood function
#'
#' This function computes the gradient of the negative log-likelihood function
#' given a set of parameters, simulated data, bandwidth vector, tau value,
#' degrees of freedom for the distribution, and locations.
#'
#' @param params a vector of variogram parameters
#' @param simu a dataframe of simulated data
#' @param h_vect a vector of spatial lag values
#' @param tau a vector of temporal lag values
#' @param df_dist a numeric value representing the degrees of freedom for the
#'                distribution
#' @param locations a matrix of locations
#'
#' @return a vector representing the gradient of the negative
#'         log-likelihood function
#'
#' @export
grad_neg_ll <- function(params, simu, h_vect, tau, df_dist, locations, # nolint
                       latlon = FALSE, quantile = 0.9, nmin = 5,
                       simu_exp = FALSE, excesses = NULL) {
  grad <- numeric(length(params))
  hmax <- max(h_vect)

  lower.bound <- c(1e-6, 1e-6, 1e-6, 1e-6)
  upper.bound <- c(Inf, Inf, 1.999, 1.999)
  if (length(params) == 6) {
    lower.bound <- c(lower.bound, -Inf, -Inf)
    upper.bound <- c(upper.bound, Inf, Inf)
  }

  if (any(params < lower.bound) || any(params > upper.bound)) {
    return(rep(0, length(params)))
  }

  if (is.null(excesses)) {
    excesses <- empirical_excesses(simu, quantile, tau, h_vect, df_dist,
                              nmin)
  }

  df_dist_new <- df_dist

  if (length(params) == 6) { # if advection is present
    adv <- params[5:6]
    dist_mat <- get_dist_mat(locations, adv = adv, tau = tau, latlon = latlon)
    df_dist_new <- reshape_distances(dist_mat)
    h_vect <- get_h_vect(df_dist_new, hmax)
    excesses <- empirical_excesses(simu, quantile, tau, h_vect,
                                   df_dist_new, nmin = nmin)
  }

  N_vect <- excesses$N_vect
  n_vect <- excesses$n_vect
  chi <- theorical_chi_mat(params, h_vect, tau)
  chi_vect <- get_chi_vect(chi, h_vect, tau, df_dist_new)

  non_excesses <- N_vect - n_vect
  # ll_vect <- n_vect * log(chi_vect) + non_excesses * log(1 - chi_vect)

  # Compute the gradient
  for (i in 1:length(params)) {
    dchi_dparam <- matrix(0, nrow = length(tau), ncol = length(h_vect))

    for (j in seq_along(h_vect)) {
      for (t in seq_along(tau)) {
        if (i == 1) {
          dchi_dparam[t, j] <- dchi_dbeta1(params, h_vect[j], tau[t])
        } else if (i == 2) {
          dchi_dparam[t, j] <- dchi_dbeta2(params, h_vect[j], tau[t])
        } else if (i == 3) {
          dchi_dparam[t, j] <- dchi_dalpha1(params, h_vect[j], tau[t])
        } else if (i == 4) {
          dchi_dparam[t, j] <- dchi_dalpha2(params, h_vect[j], tau[t])
        }
      }
    }

    dchi_dparam_vect <- get_chi_vect(dchi_dparam, h_vect, tau, df_dist_new)
    grad[i] <- -sum((n_vect / chi_vect - non_excesses / (1 - chi_vect)) * dchi_dparam_vect, na.rm = TRUE) # nolint
  }

  return(grad)
}


# SIMU -------------------------------------------------------------------------

#' simulate_excess_ind function
#'
#' Simulate individual excess following a binomial distribution with a given
#' probability parameter equals to the chi value.
#'
#' @param Tmax The maximum time observation.
#' @param chi_h_t The spatio-temporal extremogram value at a given spatial lag h
#' and a given temporal lag t.
#'
#' @return The simulated individual excess.
#'
#' @import stats
#'
#' @export
simulate_excess_ind <- function(Tmax, chi_h_t) {
  # simulate excesses following a binomial distribution
  vect_E <- rbinom(Tmax, 1, chi_h_t)
  return(vect_E)
}

#' Simulate excesses
#'
#' This function simulates excesses following a binomial distribution with given
#' probability parameters equal to the spatio-temporal extremogram values.
#'
#' @param Tmax The maximum time observation.
#' @param tau The temporal lag vector.
#' @param h_vect The spatial lag vector.
#' @param chi The spatio-temporal extremogram matrix.
#' @param df_dist The distances long dataframe.
#'
#' @return The simulated excesses matrix and the chi vector.
#'
#' @export
simulate_excesses <- function(Tmax, tau, h_vect, chi, df_dist) {
  df_dist$h <- ifelse(df_dist$value %in% h_vect, df_dist$value, NA)
  # create excesses matrix
  npairs <- nrow(df_dist) # number of pairs
  nconfig <- npairs * length(tau) # number of configurations
  mat_E <- matrix(nrow = Tmax, ncol = nconfig) # init excesses matrix
  h_vect_all <- rep(df_dist$h, times = length(tau)) # repeat h for each tau
  tau_all <- rep(tau, each = npairs) # repeat tau for each pair
  chi_vect <- numeric(length = nconfig) # init chi vector
  # for each configuration
  for (p in 1:nconfig) {
    h <- h_vect_all[p] # get h for each configuration
    if (is.na(h)) { # if NA then chi = NA and excesses = NA
      mat_E[, p] <- rep(NA, Tmax)
      chi_vect[p] <- NA
    } else {
      tau <- tau_all[p] # get tau for each configuration
      h_ind <- which(h_vect == h) # get index of h
      if (chi[tau, h_ind] <= 0) {
        # if chi <= 0 then chi = 0.000001 and excesses = 0
        mat_E[, p] <- rep(0, Tmax)
        chi_vect[p] <- 0.000001 # avoid future error
      } else {
        vect_E <- simulate_excess_ind(Tmax, chi[tau, h_ind]) # simulate excesses
        mat_E[, p] <- vect_E # add excesses to matrix
        chi_vect[p] <- chi[tau, h_ind] # add chi value to vector
      }
    }
  }
  return(list(mat_E = mat_E, chi_vect = chi_vect))
}

# VALIDATION -------------------------------------------------------------------

#' Evaluate optimization using simulated binomial experiments
#'
#' This function evaluates the optimization process using simulated experiments.
#' Simulation of excesses is performed considering a binomial distribution with
#' a probability parameter equals to the spatio-temporal extremogram value.
#'
#' @param n_res The number of realizations
#' @param Tmax The maximum time observation
#' @param tau_vect The temporal lag vector
#' @param h_vect The spatial lag vector
#' @param chi The spatio-temporal extremogram matrix
#' @param df_dist The distances long dataframe
#' @param nconfig The number of configurations
#'
#' @import spam
#' @import stats
#'
#' @return The result of the optimization process as a dataframe.
#'
#' @export
evaluate_optim_simuExp <- function(n_res, Tmax, tau_vect, h_vect, chi, df_dist,
                                   nconfig) {
  df_result <- data.frame(beta1 = rep(NA, n_res), beta2 = rep(NA, n_res),
                      alpha1 = rep(NA, n_res), alpha2 = rep(NA, n_res))

  for (n in 1:n_res){
    simu_E <- simulate_excesses(Tmax, tau_vect, h_vect, chi, df_dist)
    mat_E <- simu_E$mat_E
    N_vect <- rep(Tmax, nconfig)
    n_vect <- colSums(mat_E)
    excesses <- list(N_vect = N_vect, n_vect = n_vect) # simu of excesses

    # Conjugate gradient method
    result <- optim(par = c(0.4, 0.2, 1.5, 1), fn = neg_ll,
                    excesses = excesses,
                    h_vect = h_vect, tau = tau_vect, df_dist = df_dist,
                    method = "CG")

    params <- result$par
    df_result$beta1[n] <- params[1]
    df_result$beta2[n] <- params[2]
    df_result$alpha1[n] <- params[3]
    df_result$alpha2[n] <- params[4]
  }

    return(df_result)
}

#' evaluate_optim function
#'
#' This function evaluates the optimization process on simulated data (for
#' example, following a Brown-Resnick process) using the negative
#' log-likelihood function.
#'
#' @param list_simu A list of simulated data
#' @param quantile The quantile value
#' @param true_param The true variogram parameter (beta1, beta2, alpha1, alpha2)
#' @param tau The temporal lag vector
#' @param df_dist The distances long dataframe
#' @param locations The locations dataframe
#' @param hmax The maximum spatial lag value. Default is sqrt(17).
#' @param method The optimization method to use "CG", "Nelder-Mead",
#'               "BFGS", "SANN", "Brent" (default is "CG")
#' @param nmin The minimum number of observations to consider
#' @param parscale The scaling parameter for the optimization process. Default
#'                 is c(1, 1, 1, 1).
#' @param latlon A boolean value to indicate if the locations are in latitude
#'              and longitude. Default is FALSE.
#'
#' @return The result of the optimization process as a dataframe.
#'
#' @import spam
#' @import stats
#' @import optimx
#'
#' @export
evaluate_optim <- function(list_simu, quantile, true_param, tau, hmax,
                           locations, nmin = 5,
                           parscale = c(1, 1, 1, 1), latlon = FALSE) {

  lower.bound <- c(1e-6, 1e-6, 1e-6, 1e-6)
  upper.bound <- c(Inf, Inf, 1.999, 1.999)
  if (length(true_param) == 6) {
    lower.bound <- c(lower.bound, -1e-6, -1e-6)
    upper.bound <- c(upper.bound, Inf, Inf)
    parscale <- c(parscale, 1, 1)
  }
  # get the number of simulations
  n_res <- length(list_simu)
  # create a dataframe to store the results
  df_result <- data.frame(beta1 = rep(NA, n_res), beta2 = rep(NA, n_res),
                          alpha1 = rep(NA, n_res), alpha2 = rep(NA, n_res))

  # if there is advection
  if (length(true_param) == 6) {
    df_result$adv1 <- rep(NA, n_res)
    df_result$adv2 <- rep(NA, n_res)
  }

  h_vect <- get_lag_vectors(locations, true_param, tau = tau, hmax = hmax)
  count_cv <- 0
  # for all simulations
  for (n in 1:n_res) {
    simu_df <- as.data.frame(list_simu[[n]]) # get the simulation dataframe
    # get the empirical excesses
    # excesses <- empirical_excesses(simu_df, quantile, tau, h_vect, df_dist,
    #                                nmin)
    # optimize the negative log-likelihood function
    tryCatch({
        result <- optimr(par = true_param, method = "Rcgmin",
                  gr = "grfwd", fn = function(par) {
                  neg_ll(par, simu = simu_df, quantile = quantile,
                        h_vect = h_vect, tau = tau,
                        locations = locations, latlon = latlon,
                        nmin = nmin, excesses = NULL)
                  }, lower = lower.bound, upper = upper.bound,
                  control = list(parscale = parscale,
                                 maxit = 10000))

        if (result$convergence == 0) { # if it converges
          count_cv <- count_cv + 1
          params <- result$par
          df_result$beta1[n] <- params[1]
          df_result$beta2[n] <- params[2]
          df_result$alpha1[n] <- params[3]
          df_result$alpha2[n] <- params[4]
          if (length(true_param) == 6) {
              df_result$adv1[n] <- params[5]
              df_result$adv2[n] <- params[6]
          }
        } else {
          print(n)
        }
    }, error = function(e) {
        # Handle the error (e.g., print an error message)
        print(paste("Error occurred for simulation", n))
    })
  }
  print(paste0("Number of convergence: ", count_cv))
  return(df_result)
}


#' get_criterion function
#'
#' Calculates the mean, RMSE and MAE criterion values for the given dataframe
#' result and the true parameter values.
#'
#' @param df_result The data frame containing the result values
#' @param true_param The true variogram parameter (beta1, beta2, alpha1, alpha2)
#' @return The calculated criterion value in a dataframe.
#'
#' @export
get_criterion <- function(df_result, true_param) {
    # remove NA values
    df_result <- na.omit(df_result)
    # get the mean, RMSE and MAE values for each parameter
    mean_beta1 <- mean(df_result$beta1)
    mean_alpha1 <- mean(df_result$alpha1)
    mean_beta2 <- mean(df_result$beta2)
    mean_alpha2 <- mean(df_result$alpha2)

    rmse_beta1 <- sqrt(mean((true_param[1] - df_result$beta1)^2))
    rmse_alpha1 <- sqrt(mean((true_param[3] - df_result$alpha1)^2))
    rmse_beta2 <- sqrt(mean((true_param[2] - df_result$beta2)^2))
    rmse_alpha2 <- sqrt(mean((true_param[4] - df_result$alpha2)^2))

    mae_beta1 <- mean(abs(true_param[1] - df_result$beta1))
    mae_alpha1 <- mean(abs(true_param[3] - df_result$alpha1))
    mae_beta2 <- mean(abs(true_param[2] - df_result$beta2))
    mae_alpha2 <- mean(abs(true_param[4] - df_result$alpha2))

    # create a dataframe to store the results
    df_valid <- data.frame(mean = c(mean_beta1, mean_beta2, mean_alpha1,
                                mean_alpha2),
                        rmse = c(rmse_beta1, rmse_beta2, rmse_alpha1,
                                rmse_alpha2),
                        mae = c(mae_beta1, mae_beta2, mae_alpha1, mae_alpha2),
                        row.names = c("beta1", "beta2", "alpha1", "alpha2"))
    return(df_valid)
}