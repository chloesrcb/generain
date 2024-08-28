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
empirical_excesses <- function(data_rain, quantile, h_vect) {
  q <- quantile # quantile

  unique_tau <- unique(h_vect$tau) # unique temporal lags

  for (t in unique_tau) { # loop over temporal lags
    df_h_t <- h_vect[h_vect$tau == t, ] # get the dataframe for each lag

    for (i in seq_len(nrow(df_h_t))) { # loop over each pair of sites
      # get the indices of the sites
      ind_s2 <- as.numeric(as.character(df_h_t$s2[i]))
      ind_s1 <- df_h_t$s1[i]

      # get the data for the pair of sites
      rain_cp <- data_rain[, c(ind_s1, ind_s2), drop = FALSE]
      rain_cp <- na.omit(rain_cp)
      colnames(rain_cp) <- c("s1", "s2")

      Tmax <- nrow(rain_cp) # number of time steps
      rain_nolag <- rain_cp$s1[1:(Tmax - t)] # get the data without lag
      rain_lag <- rain_cp$s2[(1 + t):Tmax] # get the data with lag

      n <- length(rain_nolag) # number of observations
      # transform the data in uniform data
      rain_unif <- cbind(rank(rain_nolag) / (n + 1), rank(rain_lag) / (n + 1))
      # get the conditional excesses on s2
      cp_cond <- rain_unif[rain_unif[, 2] > q, , drop = FALSE]
      excess_count <- sum(cp_cond[, 1] > q) # number of excesses for s1 given
                                            # those of s2
      num_cond_excesses <- nrow(cp_cond) # number excesses for s2

      # store the number of excesses
      h_vect$N_vect[h_vect$s1 == ind_s1 & h_vect$s2 == ind_s2 & h_vect$tau == t] <- num_cond_excesses
      h_vect$n_vect[h_vect$s1 == ind_s1 & h_vect$s2 == ind_s2 & h_vect$tau == t] <- excess_count
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
  # if (length(params) == 6) {
  #   hnorm <- 
  # }
  # Get vario and chi for each lagtemp
  varioval <- 2 * (beta1 * h^alpha1 + beta2 * tau^alpha2)
  phi <- pnorm(sqrt(0.5 * varioval))
  chival <- 2 * (1 - phi)

  return(chival)
}

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

  # for each spatial lag
  h_vect_unique <- unique(h_vect$hnorm)
  for (j in seq_along(h_vect_unique)) {
      # Get vario and chi for each lagtemp
      h <- h_vect_unique[j]
      for (t in seq_along(tau)) {
          chi[t, j] <- theorical_chi_ind(params, h, tau[t])
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
#' @param df_lags A dataframe with spatial lag values.
#' @param tau A vector of temporal lag values.
#'
#' @return The theoretical chi matrix.
#'
#' @export
theorical_chi <- function(params, df_lags) {
  chi_df <- df_lags # copy the dataframe
  chi_df$chi <- theorical_chi_ind(params, df_lags$hnorm, df_lags$tau)
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
                  simu_exp = FALSE, excesses = NULL, adv1 = 0, adv2 = 0) {
  # params <- c(beta1, beta2, alpha1, alpha2)
  hmax <- max(h_vect$hnorm)
  # df_dist_new <- df_dist
  # print(params)
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
    message("out of bounds")
    return(1e8)
  }

  # h_vect_new <- h_vect

  if (length(params) == 6) {
    # adv <- params[5:6]
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
  chi <- theorical_chi(params, h_vect) # get chi matrix
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
simulate_excesses <- function(Tmax, chi) {

  excesses_df <- chi
  excesses_df$probaN <- runif(nrow(excesses_df), 0, 1)
  excesses_df$N_vect <- rep(NA, nrow(excesses_df))
  excesses_df$n_vect <- rep(NA, nrow(excesses_df))
  for (i in 1:nrow(excesses_df)) {

    excesses_df$N_vect[i] <- rbinom(1, Tmax, excesses_df$probaN[i])
    excesses_df$n_vect[i] <- rbinom(1, excesses_df$N_vect[i],
                                    excesses_df$chi[i])
  }
  return(excesses_df)
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
evaluate_optim_simuExp <- function(n_res, Tmax, tau_vect, h_vect, chi,
                                   locations, nconfig) {
  df_result <- data.frame(beta1 = rep(NA, n_res), beta2 = rep(NA, n_res),
                      alpha1 = rep(NA, n_res), alpha2 = rep(NA, n_res))

  for (n in 1:n_res){
    simu_E <- simulate_excesses(Tmax, chi)
    N_vect <-  simu_E$N_vect
    n_vect <-  simu_E$n_vect
    excesses <- list(N_vect = N_vect, n_vect = n_vect) # simu of excesses
    # Conjugate gradient method
    result <- optim(par = c(0.4, 0.2, 1.5, 1), fn = neg_ll,
                    excesses = excesses, quantile = 0.9,
                    h_vect = h_vect, tau = tau_vect,
                    locations = locations,
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
#' @param tau_vect The temporal lag vector
#' @param locations The locations dataframe
#' @param hmax The maximum spatial lag value. Default is sqrt(17).
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
#'
#' @export
evaluate_optim <- function(list_simu, quantile, true_param, tau_vect, hmax,
                           locations, nmin = 5, parscale = c(1, 1, 1, 1),
                           latlon = FALSE) {

  # lower.bound <- c(1e-6, 1e-6, 1e-6, 1e-6)
  # upper.bound <- c(Inf, Inf, 1.999, 1.999)
  if (length(true_param) == 6 && length(parscale) == 4) {
    # lower.bound <- c(lower.bound, -1e-6, -1e-6)
    # upper.bound <- c(upper.bound, Inf, Inf)
    parscale <- c(parscale, 1, 1)
  }

  n_res <- length(list_simu)
  df_result <- data.frame(beta1 = rep(NA, n_res), beta2 = rep(NA, n_res),
                          alpha1 = rep(NA, n_res), alpha2 = rep(NA, n_res))

  if (length(true_param) == 6) {
    df_result$adv1 <- rep(NA, n_res)
    df_result$adv2 <- rep(NA, n_res)
  }

  h_vect <- get_lag_vectors(locations, true_param, tau = tau_vect, hmax = hmax)
  count_cv <- 0

  cl <- makeCluster(detectCores() - 1)
  clusterEvalQ(cl, library(generain))
  # clusterExport(cl, c("list_simu", "neg_ll", "true_param",
  #                     "quantile", "tau_vect", "hmax", "locations", "nmin",
  #                     "parscale", "latlon", "h_vect", "empirical_excesses",
  #                     "optim"), envir = NULL)

  results <- parLapply(cl, 1:n_res, function(n) {
    simu_df <- as.data.frame(list_simu[[n]])

    if (length(true_param) == 6) {
      excesses <- NULL
    } else {
      excesses <- empirical_excesses(simu_df, quantile, tau_vect, h_vect, nmin)
    }

    result <- tryCatch({
      optim(par = true_param, fn = neg_ll, excesses = excesses,
            quantile = quantile,
            h_vect = h_vect, tau = tau_vect, locations = locations,
            simu = simu_df,
            method = "CG", control = list(parscale = parscale, maxit = 10000))
    }, error = function(e) {
      NULL
    })

    if (!is.null(result) && result$convergence == 0) {
      params <- result$par
      return(c(params, TRUE))
    } else {
      return(rep(NA, length(true_param) + 1))
    }
  })

  stopCluster(cl)

  for (n in 1:n_res) {
    res <- results[[n]]
    if (!is.na(res[length(res)])) {
      count_cv <- count_cv + 1
      df_result[n, 1:length(true_param)] <- res[1:length(true_param)]
    }
  }

  print(paste0("Number of convergences: ", count_cv))
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