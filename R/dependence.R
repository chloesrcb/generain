# CHI --------------------------------------------------------------------------

#' This function calculates the chi-square statistic for a given data set and
#' quantile.
#'
#' @param data The data set for which the chi-square statistic is calculated.
#' @param quantile The quantile value used in the calculation.
#'
#' @return The chi-square statistic.
#'
#' @import terra
#' @import dplyr
#' @import tidyr
#'
#' @examples
#' data1 <- c(1, 2, 3, 4, 5)
#' data2 <- c(2, 3, 4, 5, 6)
#' data <- cbind(data1, data2)
#' quantile <- 0.95
#' get_chiq(data, quantile)
#'
#' @export
get_chiq <- function(data, quantile) {
  n <- nrow(data)
  # Transform data into uniform ranks
  data_unif <- cbind(rank(data[, 1]) / (n + 1),
                     rank(data[, 2]) / (n + 1))
  # Calculate the maximum rank for each row
  rowmax <- apply(data_unif, 1, max)
  # Calculate the chi value
  u <- quantile
  cu <- mean(rowmax < u)
  chiu <- 2 - log(cu) / log(u)
  # Calculate the lower bound for chi
  chiulb <- 2 - log(pmax(2 * u - 1, 0)) / log(u)
  # Take the maximum of chi and chiulb
  chiu <- pmax(chiu, chiulb)
  return(chiu)
}



# TEMPORAL CHI -----------------------------------------------------------------


#' Calculate the temporal chi statistic
#'
#' This function calculates the temporal extremogram for a given dataset.
#'
#' @param data_rain The rainfall data.
#' @param tmax The maximum threshold value.
#' @param quantile The quantile value.
#' @param zeros Logical indicating whether to include zero values in the
#'              calculation.
#' @param mean Logical indicating whether to return the mean value.
#'
#' @return The temporal chi statistic.
#'
#' @import terra
#' @import dplyr
#' @import tidyr
#'
#' @examples
#' data <- c(1, 2, 3, 4, 5)
#' temporal_chi(data, 3, 0.5)
#'
#' @export
temporal_chi <- function(data_rain, tmax, quantile, zeros = TRUE, mean = TRUE) {
  # number of stations
  nsites <- ncol(data_rain)
  Tmax <- nrow(data_rain)
  # get maximum number of observations
  chi_s_temp <- matrix(1, nrow = nsites, ncol = tmax)
  q <- quantile
  for (s in 1:nsites) {
    rain_cp <- drop_na(data_rain[s]) # for one fixed station
    rain_Xs_unif <- data.frame(rank(rain_cp) / (nrow(rain_cp) + 1))
    mean_excess_Xs <- mean(rain_Xs_unif > q)
    for (t in 1:tmax){
      rain_nolag <- rain_cp[1:(Tmax - t), ] # without lag (in t_k)
      rain_lag <- rain_cp[(1 + t):Tmax, ] # with lag t (in t_k + t)
      data_cp <- cbind(rain_nolag, rain_lag) # get couple
      n <- nrow(data_cp)
      rain_unif <- cbind(rank(data_cp[, 1]) / (n + 1),
                        rank(data_cp[, 2]) / (n + 1))

      # check excess above a threshold q
      excess_t <- mean(rain_unif[, 1] > q & rain_unif[, 2] > q)
      chival <- excess_t / mean_excess_Xs
      chi_s_temp[s, t] <- chival
    }
    print(paste0(s, "/", nsites)) # print progress
  }
  if (mean) {
    chi_temp <- colMeans(chi_s_temp, na.rm = TRUE)
  } else {
    chi_temp <- chi_s_temp
  }
  chi_temp[chi_temp <= 0] <- 0.0000001
  return(chi_temp)
}


#' Calculate temporal chi using WLSE method
#'
#' This function calculates the temporal chi using the Weighted Least Squares
#' Estimation (WLSE) method.
#'
#' @param dftemp The input data frame containing the temporal chi values.
#' @param weights The weights to be used in the calculation: "residuals", "exp",
#'               or "none"
#' @return The calculated temporal chi value.
#'
#' @import terra
#' @import dplyr
#' @import tidyr
#' @import lmtest
#'
#' @examples
#' data <- data.frame(time = c(1, 2, 3), value = c(10, 20, 30))
#' weights <- c(0.5, 0.3, 0.2)
#' result <- temporal_chi_WLSE(data, weights)
#' print(result)
#'
#' @export
temporal_chi_WLSE <- function(dftemp, weights){
    # define weights to use (classic with residuals and Buhl one with exp)
    if (weights == "residuals") {
      # LS reg (clasical model)
      model <- lm(tchi ~ lagtemp, data = dftemp)
      wt <- 1 / lm(abs(model$residuals) ~ model$fitted.values)$fitted.values^2
    } else if (weights == "exp") {
      wt <- exp(-(dftemp$lagtemp)^2)
    } else if (weights == "none") {
      wt <- rep(1, length(dftemp$lagtemp))
    }

    # perform weighted least squares regression
    wls_model_temp <- lm(tchi ~ lagtemp, data = dftemp, weights = wt)

    # view summary of model
    sum_wls_temp <- summary(wls_model_temp)
    return(sum_wls_temp)
}



#' Calculate the estimate of variotemp
#'
#'
#' This function calculates the estimate of variotemp based on the given
#' parameters.
#'
#' @param chitemp The chitemp parameter.
#' @param tmax The tmax parameter.
#' @param npoints The number of points parameter.
#' @param weights The weights parameter.
#'
#' @return The estimate of variotemp.
#'
#' @import terra
#' @import dplyr
#' @import tidyr
#' @import lmtest
#'
#' @examples
#' get_estimate_variotemp(chitemp = 0.5, tmax = 100, npoints = 10,
#'                        weights = c(0.2, 0.3, 0.5))
#' 
#' @export
get_estimate_variotemp <- function(chitemp, tmax, npoints, weights,
                                  summary = FALSE) {
  # every chi lagged mean
  chi_df_dt <- as.data.frame(chitemp)
  # mean_chi_lag <- apply(chi_df_dt, 2, mean)

  # regression on temporal chi
  chitemp <- c(t(chi_df_dt))
  lagtemp <- c(1:tmax) # minutes

  # need to add a eta continuous function for the WLSE with log
  dftemp <- data.frame(tchi = eta(chitemp), lagtemp = log(lagtemp))

  sum_wls_temp <- temporal_chi_WLSE(dftemp, weights)

  df_wls_temp <- data.frame(sum_wls_temp$coefficients)

  # final temporal variogram parameters
  alpha_temp <- df_wls_temp$Estimate[2]
  c_temp <- df_wls_temp$Estimate[1]
  theta_temp <- exp(c_temp)
  if (summary) {
    return(list(theta_temp, alpha_temp, sum_wls_temp))
  } else {
    return(c(theta_temp, alpha_temp))
  }
}


# SPATIAL CHI ------------------------------------------------------------------

#' Calculate the spatial mean lags
#'
#' This function calculates the spatial mean lags based on a given radius.
#' The function takes an optional argument 'mid' which determines whether
#' the midpoint of the lags should be returned instead of the full lags.
#'
#' @param radius The radius used to calculate the spatial mean lags.
#' @param mid Logical value indicating whether to return the midpoint of the 
#'            lags. Default is FALSE.
#'
#' @return A vector of spatial mean lags or the midpoint of the lags if 'mid' 
#'         is TRUE.
#' 
#' @import terra
#' @import dplyr
#' @import tidyr
#' 
#' @examples
#' spatial_mean_lags(5)
#' spatial_mean_lags(10, mid = TRUE)
#'
#' @export
spatial_mean_lags <- function(radius, mid = FALSE) {
  h_vect <- c()
  for (i in c(2:length(radius))) {
    h_mean <- (radius[i - 1] + radius[i]) / 2
    h_vect <- c(h_vect, h_mean)
  }
  return(h_vect)
}



#' Function to estimate the spatial extremogram
#'
#' This function calculates the spatial extremogram for a given lags vector and
#' and radius matrix and rainfall data.
#' 
#' @param lags The lags vector.
#' @param rad_mat The matrix of radii.
#' @param data_rain The rainfall data.
#' @param quantile The quantile value.
#' @param zeros Logical value indicating whether to include zero values in the
#'             calculation.
#' @param mid Logical value indicating whether to return the mean values for
#'           the plot.
#'
#' @return The spatial extremogram.
#'
#' @import terra
#' @import dplyr
#' @import tidyr
#'
#' @examples
#' lags <- c(1, 2, 3, 4, 5)
#' rad_mat <- matrix(c(0, 1, 2, 1, 0, 3, 2, 3, 0), nrow = 3)
#' data_rain <- c(10, 20, 30)
#' quantile <- 0.5
#' spatial_chi(lags, rad_mat, data_rain, quantile)
#' 
#' @export
spatial_chi <- function(lags, rad_mat, data_rain, quantile, zeros = TRUE, 
                        mid = TRUE) {
  # initialize values
  chi_slag <- c()
  chi_val <- c()
  # get mid values for each intervals for the plot
  if (mid) {
    h_vect <- spatial_mean_lags(lags) # get mean values for the plot
  } else {
    h_vect <- lags
  }
  lags_sup <- lags[-1] # we don't need the 0, only superior born
  for (h in h_vect){
    # station number inside h lag
    ind_h <- which(h_vect == h)
    h_sup <- lags_sup[ind_h]
    indices <- data.frame(which(rad_mat == h_sup, arr.ind = TRUE))
    print(paste0("h = ", h))
    nb_pairs <- dim(indices)[1]
    if (nb_pairs == 0) {
      chi_slag <- c(chi_slag, NA)
    } else {
      # get index pairs
      ind_s1 <- indices$row
      ind_s2 <- indices$col
      chi_val <- c()
      for (i in 1:nb_pairs){
        rain_cp <- drop_na(data_rain[, c(ind_s1[i], ind_s2[i])])
        colnames(rain_cp) <- c("s1", "s2")
        if (!zeros) {
          rain_cp <- rain_cp[rowSums(rain_cp != 0, na.rm = TRUE) > 0, ]
        }
        print(c(ind_s1[i], ind_s2[i]))
        if (length(quantile)>1) {
          q <- quantile[ind_s1[i], ind_s2[i]]
        }
        chi_val <- c(chi_val, get_chiq(rain_cp, q))
      }
      chi_slag <- c(chi_slag, mean(na.omit(chi_val)))
    }
    print(chi_slag)
  }
  chispa_df <- data.frame(chi = chi_slag, lagspa = h_vect)
  return(chispa_df)
}


#' Calculate spatial chi-squared distances
#'
#' This function calculates the spatial extremogram for a given
#' distance matrix and rainfall data.
#'
#' @param df_dist The distance dataframe.
#' @param data_rain The rainfall data.
#' @param quantile The quantile value.
#' @param hmax The maximum spatial lag value (optional).
#' @param comephore Logical value indicating whether we work on the COMEPHORE 
#'                  dataset. Default is FALSE.
#'
#' @return The spatial chi-squared distances.
#'
#' @examples
#' df_dist <- matrix(c(0, 1, 2, 1, 0, 3, 2, 3, 0), nrow = 3)
#' data_rain <- c(10, 20, 30)
#' quantile <- 0.5
#' spatial_chi_alldist(df_dist, data_rain, quantile)
#'
#' @import terra
#' @import dplyr
#' @import tidyr
#'
#' @export
spatial_chi_alldist <- function(df_dist, data_rain, quantile, hmax = NA,
                                comephore = FALSE) {
  chi_slag <- c()
  # initialize values
  if (comephore) {
    df_dist$value <- ceiling(df_dist$value / 100) * 100 / 1000 # in km
  }
  # get unique distances from rad_mat
  h_vect <- sort(df_dist$value)
  h_vect <- unique(h_vect[h_vect > 0])
  if (!is.na(hmax)) {
    h_vect <- h_vect[h_vect <= hmax]
  }

  for (h in h_vect) {
    print(paste0("h = ", h))
    chi_val <- c()
    # get index pairs
    df_dist_h <- df_dist[df_dist$value == h, ]
    ind_s2 <- as.numeric(as.character(df_dist_h$X))
    ind_s1 <- df_dist_h$Y
    for (i in seq_along(ind_s1)){
      rain_cp <- drop_na(data_rain[, c(ind_s1[i], ind_s2[i])])
      colnames(rain_cp) <- c("s1", "s2")
      print(c(ind_s1[i], ind_s2[i]))
      if (length(quantile) > 1) {
        q <- quantile[ind_s1[i], ind_s2[i]]
      }
      chi_val <- c(chi_val, get_chiq(rain_cp, q))
    }
    chi_slag <- c(chi_slag, mean(na.omit(chi_val)))
  }

  chispa_df <- data.frame(chi = chi_slag, lagspa = h_vect)
  return(chispa_df)
}


#' Calculate the estimate of the spatial variogram parameters
#'
#' This function calculates the estimate of the spatial variogram parameters
#' based on the given spatial chi values. It uses the Weighted Least Squares
#' Estimation (WLSE) method to estimate the parameters.
#'
#' @param chispa The spatial chi values.
#' @param weights The weights to be used in the calculation: "residuals", "exp",
#'               or "none"
#' @param summary Logical value indicating whether to return the summary of the
#'               model. Default is FALSE.
#'
#' @return The estimated spatial variogram parameters.
#'
#' @import terra
#' @import dplyr
#' @import tidyr
#' @import lmtest
#' 
#' @examples
#' chispa <- data.frame(chi = c(0.5, 0.6, 0.7), lagspa = c(1, 2, 3))
#' weights <- "residuals"
#' get_estimate_variospa(chispa, weights)
#' @export
get_estimate_variospa <- function(chispa, weights, summary = FALSE) {
  # eta transformation
  etachispa_df <- data.frame(chi = eta(chispa$chi),
                           lagspa = log(chispa$lagspa))

  if (weights == "residuals") {
    # LS reg (clasical model)
    model <- lm(chi ~ lagspa, data = etachispa_df)
    # define weights to use
    ws <- 1 / lm(abs(model$residuals) ~ model$fitted.values)$fitted.values^2
  } else if (weights == "exp") {
    ws <- exp(-(etachispa_df$lagspa)^2)
  } else if (weights == "none") {
    ws <- rep(1, length(etachispa_df$lagspa))
  }

  # perform weighted least squares regression
  wls_model_spa <- lm(chi ~ lagspa, data = etachispa_df, weights = ws)

  # view summary of model
  sum_wls_lag <- summary(wls_model_spa)

  wls_coef <- data.frame(sum_wls_lag$coefficients)

  alpha_spa <- wls_coef$Estimate[2]
  c_spa <- wls_coef$Estimate[1]
  theta_spa <- exp(c_spa)
  if (summary) {
    return(list(theta_spa, alpha_spa, sum_wls_lag))
  } else {
    return(c(theta_spa, alpha_spa))
  }
}


# VARIOGRAM --------------------------------------------------------------------

#' Transforms the extremogram using the eta transformation.
#'
#' The eta function transforms the input vector by applying a logarithmic
#' transformation using the quantile function of the standard normal
#' distribution.
#'
#' @param chi The input extremogram to be transformed.
#' @return The transformed vector.
#'
#' @examples
#' eta(c(-1, 0, 1, 2))
#' # Output: [1]  2.000000e+00  1.000000e-06 -1.151293e-01 -2.000000e-01
#'
#' @import stats
#'
#' @export
eta <- function(chi) {
  chi[chi <= 0] <- 0.000001 # to avoid log(0)
  stdnorm <- qnorm(1 - 0.5 * chi)
  chi <- 2 * log(stdnorm)
  return(chi)
}

#' Calculate the theorical spatio-temporal variogram
#'
#' This function calculates the theorical spatio-temporal variogram based on the
#' given extremogram.
#'
#' @param chi The extremogram.
#'
#' @return The theorical spatio-temporal variogram.
#'
#' @import stats
#' 
#' @examples
#' chi <- 0.5
#' vario_spatemp(chi)
#'
#' @export
vario_spatemp <- function(chi) {
  invphi <- qnorm(1 - 0.5 * chi)
  return(2 * invphi^2)
}

#' Calculate the spatial or temporal variogram
#'
#' This function calculates the spatial or temporal variogram based on the given
#' WLSE estimated parameters c and alpha. (c = log(beta))
#'
#' @param x The corresponding lags.
#' @param c The c variogram parameter. (log(beta))
#' @param alpha The alpha variogram parameter.
#'
#' @return The calculated spatial or temporal variogram.
#'
#' @examples
#' x <- c(1, 2, 3, 4, 5)
#' c <- 1
#' alpha <- 0.5
#' vario <- vario(x, c, alpha)
#' print(vario)
#' 
#' @export
vario <- function(x, c, alpha) {
  return(2 * exp(c) * x^alpha)
}

#' Calculate the chi-squared variogram
#'
#' This function calculates the chi-squared variogram based on the given
#' variogram value.
#'
#' @param vario The variogram value.
#'
#' @return The theorical spatio-temporal extremogram.
#'
#' @examples
#' vario <- 0.5
#' chi_vario(vario)
#'
#' @export
chi_vario <- function(vario) {
  chi <- 2 * (1 - sqrt(vario / 2))
  return(chi)
}

#' Calculate the standard deviation of a variogram
#'
#' This function calculates the standard deviation of a variogram based on
#' the given parameters.
#'
#' @param x The input data.
#' @param vario The variogram function.
#' @param sd_c The constant term in the standard deviation formula.
#' @param sd_alpha The alpha term in the standard deviation formula.
#'
#' @return The standard deviation of the variogram.
#'
#' @examples
#' x <- c(1, 2, 3, 4, 5)
#' vario <- function(h) h^2
#' sd_c <- 1
#' sd_alpha <- 0.5
#' sd_vario(x, vario, sd_c, sd_alpha)
#'
#' @export
sd_vario <- function(x, vario, sd_c, sd_alpha) {
  sd <- vario * sqrt(sd_c^2 + sd_alpha^2) # with delta-method
  return(sd)
}

# SPATIO-TEMPORAL CHI ----------------------------------------------------------

#' chispatemp_dt function
#'
#' This function calculates the chi-squared spatial-temporal dependence test.
#'
#' @param lagt The temporal lag.
#' @param lags The spatial lags.
#' @param rad_mat The matrix of radii.
#' @param data_rain The rainfall data.
#' @param quant_mat The matrix of quantiles (optional).
#' @param q_fixed The fixed quantile value (optional).
#'
#' @return The result of the chi-squared spatial-temporal dependence test.
#'
#' @import terra
#' @import dplyr
#' @import tidyr
#' 
#' @examples
#' chispatemp_dt(1, 2, rad_mat, data_rain)
#' 
#' @export
chispatemp_dt <- function(lagt, lags, rad_mat, data_rain, quant_mat = NA,
                         q_fixed = NA) {
  if (!is.na(q_fixed)) {
    q <- q_fixed
  }
  chi_slag <- c()
  # get mid values for each intervals for the plot
  h_vect <- spatial_mean_lags(c(0, lags))
  lags_sup <- lags #[-1] # we don't need the 0, only superior born
  for (h in h_vect){
    # station number inside h lag
    ind_h <- which(h_vect == h)
    h_sup <- lags_sup[ind_h]
    indices <- data.frame(which(rad_mat == h_sup, arr.ind = TRUE))
    print(paste0("h = ", h))
    nb_pairs <- dim(indices)[1]
    if (nb_pairs == 0) {
      chi_slag <- c(chi_slag, NA)
    } else {
      # get index pairs
      ind_s1 <- indices$row
      ind_s2 <- indices$col
      chi_val_h <- c()
      for (i in 1:nb_pairs){
        rain_cp <- drop_na(data_rain[, c(ind_s1[i], ind_s2[i])])
        colnames(rain_cp) <- c("s1", "s2")
        if (is.na(q_fixed)) {
          q <- quant_mat[ind_s1[i], ind_s2[i]]
        }
        rains1 <- rain_cp$s1
        rains2 <- rain_cp$s2
        t <- 0
        while (t < lagt) {
          rains1 <- rains1[-length(rains1)] # remove the last value
          rains2 <- rains2[-1] # remove the first row
          t <- t + 1
        }
        data <- cbind(rains1, rains2) # get couple
        chi_val_h <- c(chi_val_h, get_chiq(data, q))
      }
      chi_slag <- c(chi_slag, mean(na.omit(chi_val_h)))
    }
  }
  return(chi_slag)
}


# VALIDATION -------------------------------------------------------------------

#' Calculate the mean absolute error (MAE) between true values and
#' predicted values.
#'
#' @param true_values A numeric vector of true values.
#' @param predicted_values A numeric vector of predicted values.
#' @return The mean absolute error between true values and predicted values.
#'
#' @examples
#' true <- c(1, 2, 3)
#' predicted <- c(1.5, 2.5, 3.5)
#' mae(true, predicted)
#'
#' @export
mae <- function(true_values, predicted_values) {
  return(mean(abs(true_values - predicted_values)))
}


#' Calculate the root mean squared error (RMSE) between true values and
#' predicted values.
#'
#' @param true_values The vector of true values.
#' @param predicted_values The vector of predicted values.
#' @return The root mean squared error (RMSE) between true values and
#'         predicted values.
#'
#' @examples
#' true_values <- c(1, 2, 3)
#' predicted_values <- c(1.5, 2.5, 3.5)
#' rmse <- calculate_rmse(true_values, predicted_values)
#' print(rmse)
#' # Output: 0.5
#'
#' @export
calculate_rmse <- function(true_values, predicted_values) {
  # Calculate the squared differences between true values and predicted values
  squared_diff <- (true_values - predicted_values)^2

  # Calculate the mean of squared differences
  mean_squared_diff <- mean(squared_diff)

  # Calculate the square root of mean squared differences
  rmse <- sqrt(mean_squared_diff)

  return(rmse)
}

#' Evaluate variogram estimates
#'
#' This function evaluates variogram estimates over all simulations and returns
#' the mean, the RMSE and the MAE of the estimates.
#' 
#' @param list_simu The list of dataframes containing the simulations.
#' @param quantile The quantile for the chiplot.
#' @param true_param The true parameters of the simulations (beta, alpha).
#' @param spatial Logical indicating whether to use spatial chi.
#' @param df_dist The dataframe containing the distances between the sites for
#'               spatial chi.
#' @param tmax The maximum temporal lag for temporal chi.
#' @param hmax The maximum spatial lag for spatial chi.
#'
#' @return A list containing the dataframe of valid results and the results.
#'
#' @import terra
#' @import dplyr
#' @import tidyr
#'
#' 
#' @export
evaluate_vario_estimates <- function(list_simu, quantile, true_param,
                                    spatial = TRUE, df_dist = NA,
                                    tmax = NA, hmax = NA) {
  # get the number of simulations
  n_res <- length(list_simu)
  # create a dataframe to store the results
  df_result <- data.frame(beta = rep(NA, n_res), alpha = rep(NA, n_res))
  # for all simulations
  for (n in 1:n_res) {
    simu_df <- as.data.frame(list_simu[[n]])
    if (spatial) {
      chi <- spatial_chi_alldist(df_dist, simu_df, quantile = quantile,
                                 hmax = hmax)
      params <- get_estimate_variospa(chi, weights = "exp", summary = FALSE)
    } else {
      chi <- temporal_chi(simu_df, tmax = tmax, quantile = quantile)
      params <- get_estimate_variotemp(chi, tmax, npoints = ncol(simu_df),
                                      weights = "exp", summary = FALSE)
    }
    df_result$beta[n] <- params[1]
    df_result$alpha[n] <- params[2]
  }

  mean_beta <- mean(df_result$beta)
  mean_alpha <- mean(df_result$alpha)

  rmse_beta <- sqrt(mean((true_param[1] - df_result$beta)^2))
  rmse_alpha <- sqrt(mean((true_param[2] - df_result$alpha)^2))

  mae_beta <- mean(abs(true_param[1] - df_result$beta))
  mae_alpha <- mean(abs(true_param[2] - df_result$alpha))

  df_valid <- data.frame(mean = c(mean_beta, mean_alpha),
                        rmse = c(rmse_beta, rmse_alpha),
                        mae = c(mae_beta, mae_alpha),
                        row.names = c("beta", "alpha"))

  return(list(df_valid = df_valid, results = df_result))
}
