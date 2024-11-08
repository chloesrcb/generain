# EXCESSES ---------------------------------------------------------------------

#' get_marginal_excess function
#'
#' This function calculates the number of marginal excesses above a given
#' quantile threshold.
#'
#' @param data_rain The rainfall data.
#' @param quantile The quantile value.
#' @param threshold A boolean value to indicate if the quantile variable is a
#'                threshold value and not a uniform quantile. Default is FALSE.
#' @param index The index of the site. Default is 1.
#' @param time The starting time index. Default is 1.
#'
#' @return The number of marginal excesses.
#'
#' @import stats
#'
#' @export
get_marginal_excess <- function(data_rain, quantile, threshold = FALSE,
                                index = 1, time = 1) {
  # shifted data
  Tmax <- nrow(data_rain)
  data <- data_rain[time:Tmax, ]

  Tobs <- nrow(data)
  if(!threshold) {
    rain_unif <- rank(data[, index]) / (Tobs + 1)
  } else {
    rain_unif <- data[, index]
  }
  marginal_excesses <- sum(rain_unif > quantile)
  return(marginal_excesses)
}



#' empirical_excesses_rpar function
#'
#' This function calculates the empirical excesses based on indicators above a
#' quantile threshold for the r-Pareto process.
#'
#' @param data_rain The rainfall data.
#' @param quantile The quantile value.
#' @param df_lags The dataframe with spatial and temporal lag values.
#' @param threshold A boolean value to indicate if the quantile variable is a
#'                threshold value and not a uniform quantile. Default is FALSE.
#' @param t0 The starting time. Default is 1.
#'
#' @return The empirical excesses dataframe with the number of excesses kij and
#' the number of possible excesses Tobs, with the lag values.
#'
#' @export
empirical_excesses_rpar <- function(data_rain, quantile, df_lags,
                                    threshold = FALSE, t0 = 1) {
  excesses <- df_lags # copy the dataframe
  unique_tau <- unique(df_lags$tau) # unique temporal lags
  for (t in unique_tau) { # loop over temporal lags
    df_h_t <- df_lags[df_lags$tau == t, ] # get the dataframe for each tau lag

    for (i in seq_len(nrow(df_h_t))) { # loop over each pair of sites
      # get the indices of the sites
      ind_s2 <- as.numeric(as.character(df_h_t$s2[i]))
      ind_s1 <- 1 # s0

      # get the data for the pair of sites
      rain_cp <- data_rain[, c(ind_s1, ind_s2), drop = FALSE]
      rain_cp <- as.data.frame(na.omit(rain_cp))
      colnames(rain_cp) <- c("s1", "s2")

      # shifted data
      X_s_t <- rain_cp$s2[(t0 + abs(t))] # X_{s,t0 + tau}
      nmargin <- sum(X_s_t > quantile) # 0 or 1
      # store the number of excesses and T - tau
      excesses$Tobs[excesses$s1 == ind_s1
                      & excesses$s2 == ind_s2
                      & excesses$tau == t] <- 1

      excesses$kij[excesses$s1 == ind_s1
                    & excesses$s2 == ind_s2
                    & excesses$tau == t] <- nmargin
    }
  }
  return(excesses)
}


#' empirical_excesses function
#'
#' This function calculates the empirical excesses based on indicators above a
#' quantile threshold.
#'
#' @param data_rain The rainfall data.
#' @param quantile The quantile value.
#' @param df_lags The dataframe with spatial and temporal lag values.
#' @param threshold A boolean value to indicate if the quantile variable is a
#'                 threshold value and not a uniform quantile. Default is FALSE.
#' @param type The type of the process, "rpareto" or "brownresnick".
#'             Default is "rpareto".
#' @param t0 The conditioning time, for r-pareto process. Default is 1.
#'
#' @return The empirical excesses dataframe with the number of excesses kij and
#' the number of possible excesses Tobs, with the lag values.
#'
#' @import tidyr
#'
#' @export
empirical_excesses <- function(data_rain, quantile, df_lags, threshold = FALSE,
                               type = "rpareto", t0 = 1) {
  if (type == "rpareto") {
    excesses <- empirical_excesses_rpar(data_rain, quantile, df_lags, threshold,
                t0)
  } else if (type == "brownresnick")  {
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
        rain_nolag <- rain_cp$s1[1:(Tmax - abs(t))] # get the data without lag
        rain_lag <- rain_cp$s2[(1 + abs(t)):Tmax] # get the data with lag
        Tobs <- length(rain_nolag) # number of observations for the lagged pair
                                  # i.e. T - tau

        # transform the data in uniform data
        if (!threshold) {
          rain_unif <- cbind(rank(rain_nolag) / (Tobs + 1),
                            rank(rain_lag) / (Tobs + 1))
        } else {
          rain_unif <- cbind(rain_nolag, rain_lag)
        }

        # number of joint excesses
        joint_excesses <- sum(rain_unif[, 2] > quantile &
                              rain_unif[, 1] > quantile)

        # store the number of excesses and T - tau
        excesses$Tobs[excesses$s1 == ind_s1
                        & excesses$s2 == ind_s2
                        & excesses$tau == t] <- Tobs

        excesses$kij[excesses$s1 == ind_s1
                      & excesses$s2 == ind_s2
                      & excesses$tau == t] <- joint_excesses
      }
    }
  } else {
    print("The variable 'type' is not valid. It has to be 'rpareto' 
          or 'brownresnick'.")
  }
  return(excesses)
}

# THEORICAL CHI ----------------------------------------------------------------

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
#' @return The theoretical chi
#'
#' @export
theorical_chi <- function(params, df_lags) {
  beta1 <- params[1]
  beta2 <- params[2]
  alpha1 <- params[3]
  alpha2 <- params[4]
  if (length(params) == 6) {
    adv <- params[5:6]
  } else {
    adv <- c(0, 0)
  }

  chi_df <- df_lags[c("s1", "s2", "tau")]
  # Get vario and chi for each lagtemp
  chi_df$hx <- df_lags$hx - adv[1] * df_lags$tau
  chi_df$hy <- df_lags$hy - adv[2] * df_lags$tau
  chi_df$hnorm <- sqrt(chi_df$hx^2 + chi_df$hy^2)
  # chi_df$hnorm <- norm_Lp(chi_df$hy, chi_df$hx, p = alpha1)

  chi_df$vario <- (2 * beta1) * chi_df$hnorm^alpha1 +
                  (2 * beta2) * abs(chi_df$tau)^alpha2

  chi_df$chi <- 2 * (1 - pnorm(sqrt(0.5 * chi_df$vario)))
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
#' @param quantile The quantile value.
#' @param excesses The excesses dataframe with the number of excesses kij and
#'                 the number of possible excesses Tobs.
#' @param hmax The maximum spatial lag value. Default is NA.
#' @param s0 The conditioning location. Default is NA.
#' @param t0 The conditioning time. Default is NA.
#' @param threshold A boolean value to indicate if the quantile variable is a
#'                threshold value and not a uniform quantile. Default is FALSE.
#'
#' @return The negative log-likelihood value.
#'
#' @export
neg_ll <- function(params, data, df_lags, quantile, excesses, hmax = NA,
                  s0 = NA, t0 = NA, threshold = FALSE, pmarg = NA) {
  Tmax <- nrow(data) # number of total observations
  # print(params)
  if (all(!is.na(s0))) { # if we have a conditioning location
    p <- 1 # sure excess for r-Pareto process in (s0,t0)
  } else {
    if (all(!is.na(pmarg))) {
      p <- pmarg
    } else {
      # number of marginal excesses
      nmarg <- get_marginal_excess(data, quantile, threshold)
      p <- nmarg / Tmax # probability of marginal excesses;
    }
  }

  # Bounds for the parameters
  lower.bound <- c(1e-6, 1e-6, 1e-6, 1e-6)
  upper.bound <- c(Inf, Inf, 1.999, 1.999)
  if (length(params) == 6) {
    lower.bound <- c(lower.bound, -Inf, -Inf)
    upper.bound <- c(upper.bound, Inf, Inf)
  }

  # Check if the parameters are in the bounds
  if (any(params < lower.bound) || any(params > upper.bound)) {
    return(1e50)
  }

  chi <- theorical_chi(params, df_lags) # get chi matrix
  ll_df <- df_lags # copy the dataframe
  ll_df$kij <- excesses$kij # number of excesses
  ll_df$Tobs <- excesses$Tobs
  ll_df$hnorm <- chi$hnorm
  ll_df$chi <- chi$chi
  ll_df$chi <- ifelse(ll_df$chi <= 0, 1e-10, ll_df$chi)
  ll_df$pchi <- 1 - p * ll_df$chi
  ll_df$pchi <- ifelse(ll_df$pchi <= 0, 1e-10, ll_df$pchi)

  # number of non-excesses
  ll_df$non_excesses <- ll_df$Tobs - ll_df$kij
  ll_df$ll <- ll_df$kij * log(ll_df$chi) +
              ll_df$non_excesses * log(ll_df$pchi)
  if (!is.na(hmax)) {
    ll_df <- ll_df[ll_df$hnorm <= hmax, ]
  }

  nll <- -sum(ll_df$ll, na.rm = TRUE)
  return(nll)
}


#' neg_ll_par function
#'
#' Calculate the negative log-likelihood for a given set of variogram
#' parameters by considering excess indicators following a binomial distribution
#' with a probability parameter equals to the theorical spatio-temporal chi
#' value. Used for fixing parameters in the optimization process.
#'
#' @param beta1 The beta1 parameter.
#' @param beta2 The beta2 parameter.
#' @param alpha1 The alpha1 parameter.
#' @param alpha2 The alpha2 parameter.
#' @param adv1 The advection parameter 1.
#' @param adv2 The advection parameter 2.
#' @param data The data dataframe.
#' @param df_lags The dataframe with spatial and temporal lag values.
#' @param excesses The excesses dataframe with the number of excesses kij and
#'                 the number of possible excesses Tobs.
#' @param quantile The quantile value.
#' @param hmax The maximum spatial lag value. Default is NA.
#' @param s0 The conditioning location. Default is NA.
#' @param t0 The conditioning time. Default is NA.
#' @param threshold A boolean value to indicate if the quantile variable is a
#'                threshold value and not a uniform quantile. Default is FALSE.
#'
#' @return The negative log-likelihood value.
#'
#' @export
neg_ll_par <- function(beta1, beta2, alpha1, alpha2, adv1, adv2, data, df_lags,
                   quantile, excesses, hmax = NA, s0 = NA, t0 = NA,
                   threshold = FALSE) {
  params <- c(beta1, beta2, alpha1, alpha2, adv1, adv2)
  nll <- neg_ll(params, data, df_lags, quantile, excesses,
                hmax = hmax, s0 = s0, t0 = t0,
                threshold = threshold)
  return(nll)
}


#' neg_ll_composite function
#'
#' Calculate the negative log-likelihood for a list of simulations.
#'
#' @param params Vector of variogram parameters (beta1, beta2, alpha1, alpha2).
#' @param list_simu A list of simulated data.
#' @param df_lags The dataframe with spatial and temporal lag values.
#' @param quantile The quantile value.
#' @param list_excesses A list of excesses dataframes.
#' @param s0 The starting location.
#' @param t0 The starting time.
#' @param hmax The maximum spatial lag value.
#'
#' @return The negative log-likelihood value.
#'
#' @import stats
#'
#' @export
neg_ll_composite <- function(params, list_simu, df_lags, quantile,
                    list_excesses, hmax = NA, s0 = NA,
                    t0 = NA, threshold = FALSE) {
  print(params)
  # Bounds for the parameters
  lower.bound <- c(1e-6, 1e-6, 1e-6, 1e-6)
  upper.bound <- c(Inf, Inf, 1.999, 1.999)
  if (length(params) == 6) {
    lower.bound <- c(lower.bound, -Inf, -Inf)
    upper.bound <- c(upper.bound, Inf, Inf)
  }

  # Check if the parameters are in the bounds
  if (any(params < lower.bound) || any(params > upper.bound)) {
    return(1e50)
  }

  m <- length(list_simu) # number of replicates
  nll_composite <- 0 # composite negative log-likelihood
  for (i in 1:m) {
    # extract simulation data from i-th simulation
    simu <- list_simu[[i]]
    excesses <- list_excesses[[i]]
    nll_i <- neg_ll(params, simu, df_lags, quantile, hmax = hmax,
                    excesses = excesses, s0 = s0, t0 = t0,
                    threshold = threshold)
    nll_composite <- nll_composite + nll_i
  }
  return(nll_composite)
}


#' neg_ll_composite_par function
#'
#' Calculate the negative log-likelihood for a list of simulations.
#' Used for fixing parameters in the optimization process.
#'
#' @param beta1 The beta1 parameter.
#' @param beta2 The beta2 parameter.
#' @param alpha1 The alpha1 parameter.
#' @param alpha2 The alpha2 parameter.
#' @param adv1 The advection parameter 1.
#' @param adv2 The advection parameter 2.
#' @param list_simu A list of simulated data.
#' @param df_lags The dataframe with spatial and temporal lag values.
#' @param quantile The quantile value.
#' @param list_excesses A list of excesses dataframes.
#' @param s0 The starting location.
#' @param t0 The starting time.
#' @param hmax The maximum spatial lag value.
#'
#' @return The negative log-likelihood value.
#'
#' @import stats
#'
#' @export
neg_ll_composite_par <- function(beta1, beta2, alpha1, alpha2, adv1, adv2,
                    list_simu, df_lags, quantile, list_excesses, hmax = NA, 
                    s0 = NA, t0 = NA, threshold = FALSE) {
  params <- c(beta1, beta2, alpha1, alpha2, adv1, adv2)
  nll_composite <- neg_ll_composite(params, list_simu, df_lags,
                                    quantile, list_excesses,
                                    hmax = hmax, s0 = s0, t0 = t0,
                                    threshold = threshold)
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

    if (length(true_param) == 6) {
      mean_adv1 <- mean(df_result$adv1)
      mean_adv2 <- mean(df_result$adv2)

      rmse_adv1 <- sqrt(mean((true_param[5] - df_result$adv1)^2))
      rmse_adv2 <- sqrt(mean((true_param[6] - df_result$adv2)^2))

      mae_adv1 <- mean(abs(true_param[5] - df_result$adv1))
      mae_adv2 <- mean(abs(true_param[6] - df_result$adv2))

      df_valid <- data.frame(mean = c(mean_beta1, mean_beta2, mean_alpha1,
                                  mean_alpha2, mean_adv1, mean_adv2),
                          rmse = c(rmse_beta1, rmse_beta2, rmse_alpha1,
                                  rmse_alpha2, rmse_adv1, rmse_adv2),
                          mae = c(mae_beta1, mae_beta2, mae_alpha1, mae_alpha2,
                                  mae_adv1, mae_adv2),
                          row.names = c("beta1", "beta2", "alpha1", "alpha2",
                                        "adv1", "adv2"))
    } else {
      df_valid <- data.frame(mean = c(mean_beta1, mean_beta2, mean_alpha1,
                                  mean_alpha2),
                          rmse = c(rmse_beta1, rmse_beta2, rmse_alpha1,
                                  rmse_alpha2),
                          mae = c(mae_beta1, mae_beta2, mae_alpha1, mae_alpha2),
                          row.names = c("beta1", "beta2", "alpha1", "alpha2"))
    }
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