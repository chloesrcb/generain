# EXCESSES ---------------------------------------------------------------------

#' get_marginal_excess function
#'
#' This function calculates the number of marginal excesses above a given
#' quantile threshold.
#'
#' @param data_rain The rainfall data.
#' @param quantile The quantile value.
#'
#' @return The number of marginal excesses.
#'
#' @import stats
#'
#' @export
get_marginal_excess <- function(data_rain, quantile, ind_s0 = NA, t0 = NA) {
  Tmax <- nrow(data_rain)
  if (is.na(ind_s0) || is.na(t0)){
    ind_s0 <- 1
    t0 <- 1
  }
  # shifted data
  data <- data_rain[t0:Tmax,]
  Tobs <- nrow(data)
  rain_unif <- rank(data[, ind_s0]) / (Tobs + 1)
  marginal_excesses <- sum(rain_unif > quantile)
  return(marginal_excesses)
}


#' empirical_excesses function
#'
#' This function calculates the empirical excesses based on indicators above a
#' quantile threshold.
#'
#' @param data_rain The rainfall data.
#' @param quantile The quantile value.
#' @param df_lags The dataframe with spatial and temporal lag values.
#'
#' @return The empirical excesses dataframe with the number of excesses kij and
#' the number of possible excesses Tobs, with the lag values.
#'
#' @import tidyr
#'
#' @export
empirical_excesses <- function(data_rain, quantile, df_lags) {
  excesses <- df_lags # copy the dataframe
  unique_tau <- unique(df_lags$tau) # unique temporal lags

  for (t in unique_tau) { # loop over temporal lags
    df_h_t <- df_lags[df_lags$tau == t, ] # get the dataframe for each tau lag

    for (i in seq_len(nrow(df_h_t))) { # loop over each pair of sites
      # get the indices of the sites
      ind_s2 <- as.numeric(as.character(df_h_t$s2[i]))
      ind_s1 <- df_h_t$s1[i]

      # get the data for the pair of sites
      rain_cp <- data_rain[, c(ind_s1, ind_s2), drop = FALSE]
      rain_cp <- as.data.frame(na.omit(rain_cp))
      colnames(rain_cp) <- c("s1", "s2")

      Tmax <- nrow(rain_cp) # number of total observations
      rain_nolag <- rain_cp$s1[1:(Tmax - t)] # get the data without lag
      rain_lag <- rain_cp$s2[(1 + t):Tmax] # get the data with lag

      Tobs <- length(rain_nolag) # number of observations for the lagged pair
                                 # i.e. T - tau

      # transform the data in uniform data
      rain_unif <- cbind(rank(rain_nolag) / (Tobs + 1),
                         rank(rain_lag) / (Tobs + 1))

      # get the conditional excesses on s2
      cp_cond <- rain_unif[rain_unif[, 2] > quantile, ]
      # number of joint excesses
      joint_excesses <- sum(cp_cond[, 1] > quantile)

      # store the number of excesses and T - tau
      excesses$Tobs[excesses$s1 == ind_s1
                      & excesses$s2 == ind_s2
                      & excesses$tau == t] <- Tobs

      excesses$kij[excesses$s1 == ind_s1
                    & excesses$s2 == ind_s2
                    & excesses$tau == t] <- joint_excesses
    }
  }
  return(excesses)
}


# THEORICAL CHI ----------------------------------------------------------------

#' Calculate the theoretical chi value.
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
#' This function calculates the theoretical chi dataframe, based on the given
#' variogram parameters.
#'
#' @param params A vector of variogram parameters.
#' @param df_lags A dataframe with spatial and temporal lag values.
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
#' value.
#'
#' @param params Vector of variogram parameters (beta1, beta2, alpha1, alpha2).
#' @param data The data dataframe.
#' @param df_lags The dataframe with spatial and temporal lag values.
#' @param excesses The excesses dataframe with the number of excesses kij and
#'                 the number of possible excesses Tobs.
#' @param locations The locations dataframe.
#' @param quantile The quantile value.
#' @param latlon A boolean value to indicate if the locations are in latitude
#'               and longitude. Default is FALSE.
#'
#' @return The negative log-likelihood value.
#'
#' @export
neg_ll <- function(params, data, df_lags, locations, quantile, excesses,
                   latlon = FALSE, hmax = NA, s0 = NA, t0 = NA) {
  if (is.na(hmax)) {
    hmax <- max(df_lags$hnorm)
  }

  tau <- unique(df_lags$tau)

  print(params)

  adv <- if (length(params) == 6) params[5:6] else c(0, 0)
  ind_s0 <- if (all(is.na(s0))) 1 else which(locations$Latitude == s0[1] &&
                                             locations$Longitude == s0[2])


  # Bounds for the parameters
  lower.bound <- c(1e-6, 1e-6, 1e-6, 1e-6)
  upper.bound <- c(Inf, Inf, 1.999, 1.999)
  if (length(params) == 6) {
    lower.bound <- c(lower.bound, 1e-6, 1e-6)
    upper.bound <- c(upper.bound, Inf, Inf)
  }

  # Check if the parameters are in the bounds
  if (any(params < lower.bound) || any(params > upper.bound)) {
    return(1e50)
  }

  if (!all(adv == c(0, 0))) { # if we have the advection parameters
    # then the lag vectors are different
    if (is.na(s0) && is.na(t0)) {
      df_lags <- get_lag_vectors(locations, params, hmax = hmax, tau_vect = tau)
    } else {
      df_lags <- get_conditional_lag_vectors(locations, params, hmax = hmax,
                                       tau_vect = tau, s0 = s0, t0 = t0)
      excesses <- empirical_excesses(data, quantile, df_lags)
    }
  }

  # number of marginal excesses
  T_marg <- get_marginal_excess(data, quantile, ind_s0, t0)
  Tmax <- if (is.na(t0)) nrow(data) else nrow(data) + 1 - t0
  p <- T_marg / Tmax # probability of marginal excesses
  chi <- theorical_chi(params, df_lags) # get chi matrix
  ll_df <- excesses
  ll_df$chi <- chi$chi
  ll_df$chi <- ifelse(ll_df$chi <= 0, 1e-10, ll_df$chi)
  ll_df$pchi <- 1 - p * ll_df$chi

  ll_df$non_excesses <- ll_df$Tobs - ll_df$kij # number of non-excesses
  ll_df$ll <- ll_df$kij * log(ll_df$chi) +
              ll_df$non_excesses * log(ll_df$pchi)

  nll <- -sum(ll_df$ll, na.rm = TRUE)
  return(nll)
}

#' neg_ll_composite function
#'
#' Calculate the negative log-likelihood for a list of simulations.
#'
#' @param params Vector of variogram parameters (beta1, beta2, alpha1, alpha2).
#' @param list_simu A list of simulated data.
#' @param df_lags The dataframe with spatial and temporal lag values.
#' @param locations The locations dataframe.
#' @param quantile The quantile value.
#' @param list_excesses A list of excesses dataframes.
#' @param latlon A boolean value to indicate if the locations are in latitude
#'              and longitude. Default is FALSE.
#' @param s0 The starting location.
#' @param t0 The starting time.
#' @param hmax The maximum spatial lag value.
#'
#' @return The negative log-likelihood value.
#'
#' @import stats
#'
#' @export
neg_ll_composite <- function(params, list_simu, df_lags, locations, quantile,
                  list_excesses, latlon = FALSE, s0 = NA, t0 = NA, hmax = NA) {

  nll_composite <- 0
  # number of simulations  in list_simu
  nsim <- length(list_simu)
  for (i in 1:nsim) {
    simu <- list_simu[[i]]
    excesses <- list_excesses[[i]]
    nll_i <- neg_ll(params, simu, df_lags, locations, quantile,
                    latlon = latlon, excesses = excesses, hmax = hmax, s0 = s0,
                    t0 = t0)
    nll_composite <- nll_composite + nll_i
  }
  return(nll_composite)
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
#' @param p_marg The probability of marginal excesses.
#'
#' @return The simulated individual excess.
#'
#' @import stats
#'
#' @export
simulate_excess_ind <- function(Tmax, chi_h_t, p_marg) {
  # simulate excesses following a binomial distribution
  vect_E <- rbinom(Tmax, 1, p_marg * chi_h_t)
  return(vect_E)
}

# VALIDATION -------------------------------------------------------------------

#' evaluate_optim function
#'
#' This function evaluates the optimization process on simulated data (for
#' example, following a Brown-Resnick process) using the negative
#' log-likelihood function.
#'
#' @param list_simu A list of simulated data
#' @param quantile The quantile value
#' @param true_param The true variogram parameter (beta1, beta2, alpha1, alpha2)
#' @param df_lags The dataframe with spatial and temporal lag values
#' @param locations The locations dataframe
#' @param parscale The scaling parameter for the optimization process. Default
#'                 is c(1, 1, 1, 1).
#' @param latlon A boolean value to indicate if the locations are in latitude
#'              and longitude. Default is FALSE.
#'
#' @return The result of the optimization process as a dataframe.
#'
#' @import spam
#' @import stats
#' @import doParallel
#' @import foreach
#'
#' @export
evaluate_optim <- function(list_simu, quantile, true_param, df_lags,
                           locations, parscale = c(1, 1, 1, 1),
                           latlon = FALSE) {

  if (length(true_param) == 6 && length(parscale) == 4) {
    parscale <- c(parscale, 1, 1)
  }

  n_res <- length(list_simu)
  df_result <- data.frame(beta1 = rep(NA, n_res), beta2 = rep(NA, n_res),
                          alpha1 = rep(NA, n_res), alpha2 = rep(NA, n_res))

  if (length(true_param) == 6) {
    df_result$adv1 <- rep(NA, n_res)
    df_result$adv2 <- rep(NA, n_res)
  }

  # df_lags <- get_lag_vectors(locations, true_param, tau = tau_vect, hmax = hmax)
  count_cv <- 0

  cl <- parallel::makeCluster(parallel::detectCores() - 1)
  # clusterEvalQ(cl, library(generain))
  parallel::clusterExport(cl, c("list_simu", "neg_ll", "true_param",
                      "quantile", "df_lags", "locations", "nmin",
                      "parscale", "latlon", "empirical_excesses",
                      "optim"), envir = NULL)

  results <- parallel::parLapply(cl, 1:n_res, function(n) {
    simu_df <- as.data.frame(list_simu[[n]])

    if (length(true_param) == 6) {
      excesses <- NULL
    } else {
      excesses <- empirical_excesses(simu_df, quantile, df_lags)
    }

    result <- tryCatch({
      optim(par = true_param, fn = neg_ll, excesses = excesses,
            quantile = quantile,
            df_lags = df_lags, locations = locations,
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

  parallel::stopCluster(cl)

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


#' save_results_optim function
#'
#' Save the results of the optimization process in a CSV file.
#'
#' @param result The result of the optimization process.
#' @param true_param The true variogram parameter (beta1, beta2, alpha1, alpha2)
#' @param filename The name of the file to save the results.
#'
#' @import utils
#' 
#' @export
save_results_optim <- function(result, true_param, filename) {
  if (result$convergence == 0) {
    rmse <- sqrt((result$par - true_param)^2)
    df_rmse <- data.frame(estim = result$par, rmse = rmse)
    rownames(df_rmse) <- c("beta1", "beta2", "alpha1", "alpha2", "Vx", "Vy")
    # save the results
    utils::write.csv(t(df_rmse), file = paste0("../data/simulations_BR/results/",
                          filename, ".csv"))
  } else {
    print("No convergence")
  }
}

#' get_results_optim function
#'
#' Get the results of the optimization process from a CSV file.
#'
#' @param filename The name of the file to get the results.
#'
#' @import utils
#' 
#' @export
get_results_optim <- function(filename) {
  df_rmse <- utils::read.csv(paste0("../data/simulations_BR/results/", filename,
                      ".csv"))
  return(df_rmse)
}