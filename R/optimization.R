# EXCESSES ---------------------------------------------------------------------

#' empirical_excesses function
#'
#' This function calculates the empirical excesses based on indicators above a
#' quantile threshold.
#'
#' @param data_rain The rainfall data.
#' @param quantile The quantile value.
#' @param tau The temporal lag vector.
#' @param h_vect The spatial lag vector.
#' @param df_dist The distances inside a long dataframe.
#' @param nmin The minimum number of observations to consider.
#'
#' @return A list with n_vect the number of excesses and N_vect the number of
#' possible excesses ie the number of observations for each pair of sites.
#' 
#' @import tidyr
#' 
#' @export
empirical_excesses <- function(data_rain, quantile, tau, h_vect, df_dist,
                                nmin = 5) {
  Tmax <- nrow(data_rain) # number of time steps
  q <- quantile # quantile
  df_dist$h <- ifelse(df_dist$value %in% h_vect, df_dist$value, NA)
  npairs <- nrow(df_dist)
  N_vect <- c() # number of observations
  n_vect <- c() # number of excesses (sum)
  for (t in tau) {
    for (i in 1:npairs) {
      # get index pairs
      ind_s2 <- as.numeric(as.character(df_dist$X[i]))
      ind_s1 <- df_dist$Y[i]
      rain_cp <- drop_na(data_rain[, c(ind_s1, ind_s2)])
      colnames(rain_cp) <- c("s1", "s2")
      rains1 <- rain_cp$s1
      rains2 <- rain_cp$s2

      rain_nolag <- rains1[1:(Tmax - t)] # without lag (in t_k)
      rain_lag <- rains2[(1 + t):Tmax] # with lag t (in t_k + t)
      data_cp <- cbind(rain_nolag, rain_lag) # get final couple
      n <- nrow(data_cp)
      if (n < nmin) { # if not enough data
        n_vect <- c(n_vect, NA)
        N_vect <- c(N_vect, NA)
      } else {
        rain_unif <- cbind(rank(data_cp[, 1]) / (n + 1),
                          rank(data_cp[, 2]) / (n + 1))

        # check excess above a threshold q
        cp_cond <- rain_unif[rain_unif[, 2] > q, ]

        # nb of simultaneous excesses
        excess_count <- sum(cp_cond[, 1] > q)

        n_vect <- c(n_vect, excess_count)
        N_vect <- c(N_vect, nrow(cp_cond))
      }
    }
  }
  return(list(n_vect = n_vect, N_vect = N_vect))
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
  for (j in seq_along(h_vect)) {
      # Get vario and chi for each lagtemp
      h <- h_vect[j]
      for (t in seq_along(tau)) {
          chi[t, j] <- theorical_chi_ind(params, h, tau[t])
      }
  }
  return(chi)
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
#' @param simu_exp A boolean value to indicate if the data is simulated with
#'                 an exponential distribution.
#'
#' @return The negative log-likelihood value.
#'
#' @export
neg_ll <- function(params, excesses, h_vect, tau, df_dist, simu_exp = FALSE) {
  # params conditions
  beta1 <- params[1]
  beta2 <- params[2]
  alpha1 <- params[3]
  alpha2 <- params[4]

#   if (length(params) != 4) {
#     adv <- params[5:6]
#   }

  if (abs(alpha1) < 0.001 || abs(alpha2) < 0.001 ||  abs(alpha1 - 2) < 0.001 ||
    abs(alpha2 - 2) < 0.001 || abs(beta1) < 0.00001 || abs(beta2) <= 0.00001 ||
    alpha1 <= 0 || alpha2 <= 0 || beta1 <= 0 || beta2 <= 0 ||
    alpha1 >= 2 || alpha2 >= 2) {
      return(Inf)
  }

  N_vect <- excesses$N_vect # number of observations
  n_vect <- excesses$n_vect # number of excesses
  chi <- theorical_chi_mat(params, h_vect, tau) # get chi matrix
  chi_vect <- get_chi_vect(chi, h_vect, tau, df_dist) # transform in chi vector

  logC <- lchoose(N_vect, n_vect) # log binomial coefficient
  non_excesses <- N_vect - n_vect # number of non-excesses
  # log-likelihood vector
  ll_vect <- logC + n_vect * log(chi_vect) + non_excesses * log(1 - chi_vect)
  # negative log-likelihood
  nll <- -sum(ll_vect, na.rm = TRUE)
#   print(nll)
  return(nll)
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
#' @param method The optimization method to use "CG", "Nelder-Mead",
#'               "BFGS", "SANN", "Brent" (default is "CG")
#' @param nmin The minimum number of observations to consider
#'
#' @return The result of the optimization process as a dataframe.
#'
#' @import spam
#' @import stats
#'
#' @export
evaluate_optim <- function(list_simu, quantile, true_param, tau, df_dist,
                           method = "CG", nmin = 5) {
  # get the number of simulations
  n_res <- length(list_simu)
  # create a dataframe to store the results
  df_result <- data.frame(beta1 = rep(NA, n_res), beta2 = rep(NA, n_res),
                          alpha1 = rep(NA, n_res), alpha2 = rep(NA, n_res))
  h_vect <- get_h_vect(df_dist, sqrt(17))
  # for all simulations
  for (n in 1:n_res) {
    simu_df <- as.data.frame(list_simu[[n]]) # get the simulation dataframe
    # get the empirical excesses
    excesses <- empirical_excesses(simu_df, quantile, tau, h_vect, df_dist,
                                   nmin)
    # optimize the negative log-likelihood function
    tryCatch({
        result <- optim(par = c(0.4, 0.2, 1.5, 1), fn = neg_ll,
            excesses = excesses, h_vect = h_vect, tau = tau,
            df_dist = df_dist, method = method)
        params <- result$par
        df_result$beta1[n] <- params[1]
        df_result$beta2[n] <- params[2]
        df_result$alpha1[n] <- params[3]
        df_result$alpha2[n] <- params[4]
    }, error = function(e) {
        # Handle the error (e.g., print an error message)
        print(paste("Error occurred for simulation", n))
    })
  }
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
#' @import terra
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