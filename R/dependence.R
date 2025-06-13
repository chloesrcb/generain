# CHI --------------------------------------------------------------------------

#' This function calculates the chi extremogram for a given dataset
#' and a quantile value.
#'
#' @param data The data set.
#' @param quantile The quantile value.
#'
#' @return The chi value.
#'
#' @import dplyr
#' @import tidyr
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


get_chi_empirical <- function(data, q) {
  if (nrow(data) == 0) return(NA)

  nb_joint_excess <- sum(data[,1] > q & data[,2] > q)
  nb_excess_total <- sum(data[,1] > q)

  if (nb_excess_total == 0) return(NA)

  return(nb_joint_excess / nb_excess_total)
}



#' chispatemp_empirical function
#' 
#' Calculate the spatio-temporal chi extremogram for a given dataset
#' and a quantile value.
#' 
#' @param data_rain The data set.
#' @param df_lags The dataframe containing the spatial and temporal lags.
#' @param quantile The quantile value.
#' 
#' @return The spatio-temporal chi value.
#' 
#' @import stats
#' @import dplyr
#' @import tidyr
#' 
#' @export 
chispatemp_empirical <- function(data_rain, df_lags, quantile) {
  chi_st <- df_lags
  chi_st$chiemp <- NA
  chi_st$chiemp2 <- NA
  Tmax <- nrow(data_rain)
  for (t in unique(df_lags$tau)) {
    df_h_t <- df_lags[df_lags$tau == t, ] # get the dataframe for each lag
    for (i in seq_len(nrow(df_h_t))) { # loop over each pair of sites
      # get index pairs
      ind_s1 <- df_h_t$s1[i]
      ind_s2 <- as.numeric(as.character(df_h_t$s2[i]))
      rain_cp <- na.omit(data_rain[, c(ind_s1, ind_s2)])
      colnames(rain_cp) <- c("s1", "s2")
      rain_lag <- rain_cp$s1[(t + 1):Tmax]
      rain_nolag <- rain_cp$s2[1:(Tmax - t)]
      data <- cbind(rain_lag, rain_nolag) # get couple
      Tobs <- nrow(data) # T - tau
      data_unif <- cbind(rank(data[, 2]) / (Tobs + 1),
                          rank(data[, 1]) / (Tobs + 1))

      # get the number of marginal excesses
      n_marg <- sum(data_unif[, 1] > quantile)

      # get the number of joint excesses
      cp_cond <- data_unif[data_unif[, 1] > quantile,]
      joint_excesses <- sum(cp_cond[, 2] > quantile)
      chi_hat_h_tau <- joint_excesses / n_marg
      chi_st$chiemp[chi_st$s1 == ind_s1 & chi_st$s2 == ind_s2 &
                    chi_st$tau == t] <- chi_hat_h_tau
      chi_val <- get_chiq(data_unif, quantile)
      chi_st$chiemp2[chi_st$s1 == ind_s1 & chi_st$s2 == ind_s2 &
                      chi_st$tau == t] <- chi_val
    }
  }
  # if chi <= 0 then set it to 0.000001
  chi_st$chiemp <- ifelse(chi_st$chiemp <= 0, 0.000001, chi_st$chiemp)
  chi_st$chiemp2 <- ifelse(chi_st$chiemp2 <= 0, 0.000001, chi_st$chiemp2)
  return(chi_st)
}



spatial_chi_pair <- function(rain_pair, q, zeros = TRUE) {
  if (!zeros) {
    nonzero_idx <- rowSums(rain_pair != 0, na.rm = TRUE) > 0
    rain_pair <- rain_pair[nonzero_idx, , drop = FALSE]
  }
  
  complete_idx <- complete.cases(rain_pair)
  rain_cp <- rain_pair[complete_idx, , drop = FALSE]
  
  if (nrow(rain_cp) == 0) return(NA)
  
  # Calculate empirical chi
  get_chi_empirical(rain_cp, q)
}


compute_all_pairs_chi <- function(data_rain, dist_mat, q = 0.99, zeros = TRUE) {
  n_sites <- ncol(data_rain)
  results <- data.frame(s1 = integer(0), s2 = integer(0), distance = numeric(0), chi = numeric(0))

  for (s1 in 1:(n_sites - 1)) {
    for (s2 in (s1 + 1):n_sites) {
      dist_val <- dist_mat[s1, s2]
      rain_pair <- data_rain[, c(s1, s2), drop = FALSE]

      # Filtrage des dates communes
      valid_idx <- complete.cases(rain_pair)
      if (!zeros) {
        valid_idx <- valid_idx & rowSums(rain_pair != 0, na.rm = TRUE) > 0
      }

      rain_valid <- rain_pair[valid_idx, , drop = FALSE]

      if (nrow(rain_valid) == 0) {
        chi_val <- NA
      } else {
        num_joint_excess <- sum(rain_valid[,1] > q & rain_valid[,2] > q)
        num_marginal_excess <- sum(rain_valid[,1] > q)
        chi_val <- if (num_marginal_excess == 0) NA else num_joint_excess / num_marginal_excess
      }

      results <- rbind(results, data.frame(s1 = s1, s2 = s2, distance = dist_val, chi = chi_val))
    }
  }

  return(results)
}

aggregate_chi_by_distance <- function(chi_table, nb_bins = 10, method = c("fixed", "quantile")) {
  method <- match.arg(method)
  
  dist_vals <- chi_table$distance
  if (method == "quantile") {
    breaks <- unique(quantile(dist_vals, probs = seq(0, 1, length.out = nb_bins + 1)))
  } else {
    breaks <- seq(min(dist_vals), max(dist_vals), length.out = nb_bins + 1)
  }

  # Création de la classe pour chaque paire
  chi_table$bin <- cut(chi_table$distance, breaks = breaks, include.lowest = TRUE, labels = FALSE)

  # Moyenne de chi par bin
  agg_result <- aggregate(chi ~ bin, data = chi_table, FUN = mean, na.rm = TRUE)

  # Ajouter le centre des classes de distance
  bin_mids <- (breaks[-1] + breaks[-length(breaks)]) / 2
  agg_result$distance_mid <- bin_mids[agg_result$bin]

  return(agg_result)
}

spatial_chi <- function(data_rain, dist_mat, q, zeros = TRUE,
                        nb_bins = 10, bin_method = "quantile") {
  bin_method <- match.arg(bin_method)

  chi_by_pair <- compute_all_pairs_chi(data_rain, dist_mat, q = q,
                                       zeros = zeros)

  agg_chi <- aggregate_chi_by_distance(chi_by_pair, nb_bins = nb_bins, 
                                        method = bin_method)
  chispa_df <- data.frame(lagspa = agg_chi$distance_mid,
                     chi = agg_chi$chi)
  return(chispa_df)
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
#' @import dplyr
#' @import tidyr
#'
#' @export
temporal_chi <- function(data_rain, tmax, quantile, zeros = TRUE, mean = TRUE) {
  nsites <- ncol(data_rain)
  tau_vect <- 0:tmax
  chi_s_temp <- matrix(1, nrow = nsites, ncol = length(tau_vect))
  for (s in 1:nsites) {
    rain_Xs <- na.omit(data_rain[, s])
    
    if (!zeros) {
      rain_Xs <- rain_Xs[rain_Xs > 0]
    } else if (is.list(rain_Xs)) {
      rain_Xs <- rain_Xs[[1]]
    } else {
      rain_Xs <- as.vector(rain_Xs)
    }
    
    if (length(rain_Xs) <= tmax) next  # skip if not enough data

    # transform to uniform using rank
    ranks <- rank(rain_Xs) / (length(rain_Xs) + 1)
    
    # get quantile (possibly site-specific)
    q <- if (is.matrix(quantile)) quantile[s, s] else quantile
    
    # count denominator (number of exceedances at lag 0)
    denom <- sum(ranks > q)

    for (tau in tau_vect) {
      if (tau >= length(ranks)) {
        chi_s_temp[s, tau + 1] <- NA
        next
      }
      
      R1 <- ranks[1:(length(ranks) - tau)]
      R2 <- ranks[(1 + tau):length(ranks)]
      
      numerator <- sum(R1 > q & R2 > q)
      
      chival <- if (denom > 0) numerator / denom else NA
      chi_s_temp[s, tau + 1] <- chival
    }
  }
  
  if (mean) {
    chi_temp <- colMeans(chi_s_temp, na.rm = TRUE)
  } else {
    chi_temp <- chi_s_temp
  }
  
  # restrict values to [1e-8, 1 - 1e-8]
  chi_temp[chi_temp <= 0] <- 1e-8
  chi_temp[chi_temp >= 1] <- 1 - 1e-8
  
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
#' @import dplyr
#' @import tidyr
#' @import lmtest
#'
#'
#' @export
temporal_chi_WLSE <- function(dftemp, weights){
    # check if weights is valid
    if (!(weights %in% c("residuals", "exp", "none"))) {
      stop("Invalid weights parameter. Must be 'residuals', 'exp', or 'none'.")
    }

    # define weights to use
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

    # get summary of model
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
#' @param weights The weights parameter, "residuals", "exp", or "none".
#' @param summary Logical value indicating whether to return the summary of the
#'              model. Default is FALSE.
#'
#' @return The estimate of variotemp.
#'
#' @import dplyr
#' @import tidyr
#' @import lmtest
#'
#' @export
get_estimate_variotemp <- function(df_chi, weights, summary = FALSE) {
  # regression on temporal chi without 0 lag
  if (0 %in% df_chi$lag) {
    # remove line with 0 lag
    df_chi <- df_chi[df_chi$lag != 0, ]
  }

  chitemp <- df_chi$chi
  lagtemp <- df_chi$lag

  # need to add a eta continuous function for the WLSE with log
  dftemp <- data.frame(tchi = eta(chitemp), lagtemp = log(lagtemp))

  sum_wls_temp <- temporal_chi_WLSE(dftemp, weights)

  df_wls_temp <- data.frame(sum_wls_temp$coefficients)

  # final temporal variogram parameters
  alpha_temp <- df_wls_temp$Estimate[2]
  c_temp <- df_wls_temp$Estimate[1]
  beta_temp <- exp(c_temp) # beta

  if (summary) {
    print(sum_wls_temp)
  }

  # get significance of parameters
  p_values <- df_wls_temp$`Pr...t..`
  # get significance symbols for parameters
  significance <- ifelse(p_values < 0.001, "***",
                         ifelse(p_values < 0.01, "**",
                                ifelse(p_values < 0.05, "*", "")))
  return(c(c_temp, beta_temp, alpha_temp, significance))
}

get_estimate_variotemp <- function(df_chi, weights, summary = FALSE, lag_unit = 1) {
  # regression on temporal chi without 0 lag
  if (0 %in% df_chi$lag) {
    # remove line with 0 lag
    df_chi <- df_chi[df_chi$lag != 0, ]
  }

  chitemp <- df_chi$chi
  lagtemp <- df_chi$lag

  # regression with log lag
  dftemp <- data.frame(tchi = eta(chitemp), lagtemp = log(lagtemp))

  sum_wls_temp <- temporal_chi_WLSE(dftemp, weights)
  df_wls_temp <- data.frame(sum_wls_temp$coefficients)

  alpha_temp <- df_wls_temp$Estimate[2]
  c_temp <- df_wls_temp$Estimate[1]

  # Correction de l'intercept pour que beta soit à la minute
  c_temp_minute <- c_temp + alpha_temp * log(lag_unit)
  beta_temp_minute <- exp(c_temp_minute)

  if (summary) {
    print(sum_wls_temp)
  }

  # significance
  p_values <- df_wls_temp$`Pr...t..`
  significance <- ifelse(p_values < 0.001, "***",
                         ifelse(p_values < 0.01, "**",
                                ifelse(p_values < 0.05, "*", "")))

  # retourne le c non corrigé (pour info), beta corrigé, alpha, et significativité
  return(c(c_temp, beta_temp_minute, alpha_temp, significance))
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
#' @import dplyr
#' @import tidyr
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
#' @import dplyr
#' @import tidyr
#'
#' @export
spatial_chi <- function(rad_mat, data_rain, quantile, zeros = TRUE,
                        mid = TRUE) {
  # initialize values
  chi_slag <- c()
  chi_val <- c()
  q <- quantile
  lags <- get_h_vect(rad_mat, NA, intervals = TRUE)
  # get mid values for each intervals for the plot
  if (mid) {
    h_vect <- spatial_mean_lags(lags, mid = mid) # get mean values for the plot
  } else {
    h_vect <- lags
  }
  lags_sup <- lags # we don't need the 0, only superior born
  for (h in h_vect){
    # station number inside h lag
    ind_h <- which(h_vect == h)
    h_sup <- lags_sup[ind_h]
    indices <- data.frame(which(rad_mat == h_sup, arr.ind = TRUE))
    # print(paste0("h = ", h))
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
        # print(c(ind_s1[i], ind_s2[i]))
        if (length(quantile) > 1) {
          q <- quantile[ind_s1[i], ind_s2[i]]
        }
        chi_val <- c(chi_val, get_chiq(rain_cp, q))
      }
      chi_slag <- c(chi_slag, mean(na.omit(chi_val)))
    }
    # print(chi_slag)
  }
  chispa_df <- data.frame(chi = chi_slag, lagspa = h_vect)
  return(chispa_df)
}


spatial_chi <- function(rad_mat, data_rain, breaks, quantile, zeros = TRUE, mid = TRUE) {
  chi_slag <- numeric(0)
  q <- quantile
  
  lags_inf <- breaks[-length(breaks)]
  lags_sup <- breaks[-1]
  h_vect <- if (mid) (lags_inf + lags_sup) / 2 else lags_sup

  for (i in seq_along(h_vect)) {
    cat(sprintf("h = %.2f\n", h_vect[i]))
    h_sup <- lags_sup[i]
    indices <- which(rad_mat == h_sup, arr.ind = TRUE)

    nb_pairs <- nrow(indices)
    cat(sprintf("nb pairs: %d\n", nb_pairs))

    if (nb_pairs == 0) {
      chi_slag <- c(chi_slag, NA)
    } else {
      chi_vals <- numeric(nb_pairs)

      for (j in seq_len(nb_pairs)) {
        s1 <- indices[j, 1]
        s2 <- indices[j, 2]
        rain_pair <- data_rain[, c(s1, s2), drop = FALSE]
        complete_idx <- complete.cases(rain_pair)

        if (!zeros) {
          nonzero_idx <- rowSums(rain_pair != 0, na.rm = TRUE) > 0
          complete_idx <- complete_idx & nonzero_idx
        }

        rain_cp <- rain_pair[complete_idx, , drop = FALSE]

        q_val <- if (length(q) > 1) q[s1, s2] else q
        chi_vals[j] <- get_chiq(rain_cp, q_val)
      }

      chi_mean <- mean(chi_vals, na.rm = TRUE)
      cat(sprintf("chi = %.4f\n", chi_mean))
      chi_slag <- c(chi_slag, chi_mean)
    }
  }

  data.frame(chi = chi_slag, lagspa = h_vect)
}


#' Calculate spatial chi for all lags
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
#' @import dplyr
#' @import tidyr
#'
#' @export
spatial_chi_alldist <- function(df_dist, data_rain, quantile, hmax = NA,
                                comephore = FALSE) {
  chi_slag <- c()
  q <- quantile
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
      if (!zeros) {
        rain_cp <- rain_cp[rowSums(rain_cp != 0, na.rm = TRUE) > 0, ]
      }
      print(c(ind_s1[i], ind_s2[i]))
      if (length(quantile) > 1) {
        q <- quantile[ind_s1[i], ind_s2[i]]
      }
      chi_val <- c(chi_val, get_chiq(rain_cp, q))
    }
    chi_val[chi_val <= 0] <- 1e-8
    chi_val[chi_val == 1] <- 1 - 1e-8
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
#' @param df_dist The distance dataframe.
#' @param data_rain The rainfall data.
#' @param quantile The quantile value.
#' @param hmax The maximum spatial lag value (optional).
#' @param zeros Logical value indicating whether to include zero values in the
#'             calculation.
#'
#' @return The estimated spatial variogram parameters.
#'
#' @export
spatial_chi_alldist <- function(df_dist, data_rain, quantile, hmax = NA,
                                zeros = FALSE) {
  chi_slag <- c()
  q <- quantile
  # get unique distances from rad_mat
  h_vect <- sort(df_dist$value)
  h_vect <- unique(h_vect[h_vect > 0]) # remove 0
  if (!is.na(hmax)) {
    h_vect <- h_vect[h_vect <= hmax]
  }

  for (h in h_vect) {
    chi_val <- c()
    print(h)
    # get index pairs
    df_dist_h <- df_dist[df_dist$value == h, ]
    ind_s2 <- as.numeric(as.character(df_dist_h$X)) # get index of s2
    ind_s1 <- df_dist_h$Y # get index of s1
    for (i in seq_along(ind_s1)){ # loop over each pair of sites
      rain_cp <- drop_na(data_rain[, c(ind_s1[i], ind_s2[i])]) # get couple
      colnames(rain_cp) <- c("s1", "s2") # rename columns
      if (!zeros) {
        rain_cp <- rain_cp[rowSums(rain_cp) > 0, ] # remove zeros
        # rain_cp <- rain_cp[rowSums(rain_cp != 0, na.rm = TRUE) > 0, ]
      }
      if (length(quantile) > 1) {
        q <- quantile[ind_s1[i], ind_s2[i]]
      }
      chi_val <- c(chi_val, get_chiq(rain_cp, q))
    }
    chi_val[chi_val <= 0] <- 1e-6
    chi_val[chi_val == 1] <- 1 - 1e-6
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
#' @import dplyr
#' @import tidyr
#' @import lmtest
#'
#' @export
get_estimate_variospa <- function(chispa, weights, summary = FALSE) {
  # remove 0 lag
  if (0 %in% chispa$lag) {
    # remove line with 0 lag
    chispa <- chispa[chispa$lag != 0, ]
      # df_chi$lag <- df_chi$lag + 0.0001
  }

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
  beta_spa <- exp(c_spa)
  if (summary) {
    print(sum_wls_lag)
  }

  # get significance of parameters
  p_values <- wls_coef$`Pr...t..`
  # get significance symbols for parameters
  significance <- ifelse(p_values < 0.001, "***",
                         ifelse(p_values < 0.01, "**",
                                ifelse(p_values < 0.05, "*", "")))
  return(c(c_spa, beta_spa, alpha_spa, significance))
}



#' Estimate spatial variogram parameters using WLSE
#'
#' @param chispa Data frame with columns: lagspa (in meters), chi
#' @param weights Type of weighting: "residuals", "exp", or "none"
#' @param summary Logical, print model summary or not
#' @param lag_unit Unit conversion factor for lags (e.g. 1000 for km, 1 for m)
#'
#' @return A data.frame with columns: beta, alpha, signif_beta, signif_alpha
#' @export
get_estimate_variospa <- function(chispa, weights = "exp", summary = FALSE, lag_unit = 1) {
  # Remove lag 0 if present
  if (0 %in% chispa$lagspa) {
    chispa <- chispa[chispa$lagspa != 0, ]
  }

  # Convert lag to desired unit
  chispa$lagspa <- chispa$lagspa / lag_unit

  # Eta transformation and log-lag
  df_eta <- data.frame(
    chi = eta(chispa$chi),
    lagspa = log(chispa$lagspa)
  )

  # Define weights
  ws <- switch(weights,
               "residuals" = {
                 model <- lm(chi ~ lagspa, data = df_eta)
                 1 / lm(abs(model$residuals) ~ model$fitted.values)$fitted.values^2
               },
               "exp" = exp(-df_eta$lagspa^2),
               "none" = rep(1, nrow(df_eta)),
               stop("Unknown weight type")
  )

  # Fit WLSE
  model <- lm(chi ~ lagspa, data = df_eta, weights = ws)
  model_sum <- summary(model)

  # Extract coefficients
  coefs <- data.frame(model_sum$coefficients)
  c_m <- coefs$Estimate[1]
  alpha <- coefs$Estimate[2]
  beta <- exp(c_m)
  c_km <- c_m - alpha * log(lag_unit)  # Adjust c for km unit
  beta_km <- exp(c_km)

  # Significance stars
  p_values <- coefs$`Pr...t..`
  signif <- ifelse(p_values < 0.001, "***",
                   ifelse(p_values < 0.01, "**",
                          ifelse(p_values < 0.05, "*", "")))

  if (summary) print(model_sum)

  return(data.frame(
    c = c_km,
    beta = beta_km,
    alpha = alpha,
    signif_beta = signif[1],
    signif_alpha = signif[2]
  ))
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
#' @import stats
#'
#' @export
eta <- function(chi) {
  chi[chi <= 0] <- 1e-6 # to avoid log(0)
  stdnorm <- qnorm(1 - 0.5 * chi)
  chi <- 2 * log(stdnorm)
  return(chi)
}

eta <- function(chi) {
  chi <- pmin(pmax(chi, 1e-6), 1 - 1e-6)  # clamp chi into (0,1)
  stdnorm <- qnorm(1 - 0.5 * chi)
  transformed <- 2 * log(stdnorm)
  return(transformed)
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
#' @export
sd_vario <- function(x, vario, sd_c, sd_alpha) {
  sd <- vario * sqrt(sd_c^2 + sd_alpha^2) # with delta-method
  return(sd)
}

# SPATIO-TEMPORAL CHI ----------------------------------------------------------

#' chispatemp_dt function
#'
#' This function calculates the chi spatial-temporal dependence
#'
#' @param lagt The temporal lag.
#' @param lags The spatial lags.
#' @param rad_mat The matrix of radii.
#' @param data_rain The rainfall data.
#' @param quant_mat The matrix of quantiles (optional).
#' @param q_fixed The fixed quantile value (optional).
#'
#' @return The result of the chi spatial-temporal dependence.
#'
#' @import dplyr
#' @import tidyr
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
#' @param spatial Logical indicating whether to use spatial chi.
#' @param df_dist The dataframe containing the distances between the sites for
#'               spatial chi.
#' @param tmax The maximum temporal lag for temporal chi.
#' @param hmax The maximum spatial lag for spatial chi.
#'
#' @return A list containing the dataframe of valid results and the results.
#'
#' @import lmtest
#'
#' @export
evaluate_vario_estimates <- function(list_simu, quantile,
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

  return(df_result)
}



#' Estimate spatial variogram parameters using WLSE
#'
#' This function estimates the spatial variogram parameters using a transformed
#' spatial extremogram and weighted least squares estimation (WLSE).
#'
#' @param chispa A data frame with columns:
#'        - `chi`: empirical extremogram values
#'        - `lagspa`: spatial lags (in km, m, etc.)
#' @param weights Type of weights to use: "residuals", "exp", or "none"
#' @param summary Logical, whether to print regression summary
#'
#' @return A numeric vector: c(intercept, beta, alpha)
#'         - intercept: log(beta)
#'         - beta: scale parameter
#'         - alpha: rate of spatial decay (in original lag units)
#'
#' @export
get_estimate_variospa <- function(chispa, weights, summary = FALSE) {
  # eta transformation
  # etachispa_df <- data.frame(chi = eta(chispa$chi),
  #                          lagspa = log(chispa$lagspa))

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
  beta_spa <- exp(c_spa)
  if (summary) {
    print(sum_wls_lag)
  }
  return(c(c_spa, beta_spa, alpha_spa))
}



# get_estimate_variospa <- function(chispa, weights = "none", summary = FALSE) {
#   # Remove lag = 0 if present (undefined or unstable)
#   chispa <- chispa[chispa$lagspa != 0, ]

#   # Normalize lags (to avoid instability with meters vs km)
#   max_lag <- max(chispa$lagspa)
#   chispa$lag_scaled <- chispa$lagspa / max_lag

#   # Apply eta transformation to chi values
#   etachispa_df <- data.frame(
#     chi_eta = eta(chispa$chi),
#     lag_scaled = chispa$lag_scaled
#   )

#   # Define weights
#   if (weights == "residuals") {
#     # Initial OLS to estimate variance of residuals
#     ols_model <- lm(chi_eta ~ lag_scaled, data = etachispa_df)
#     resid_model <- lm(abs(ols_model$residuals) ~ ols_model$fitted.values)
#     ws <- 1 / (fitted(resid_model)^2)
#   } else if (weights == "exp") {
#     # Exponential decay on normalized lag
#     ws <- exp(-(etachispa_df$lag_scaled)^2)
#   } else if (weights == "none") {
#     ws <- rep(1, nrow(etachispa_df))
#   } else {
#     stop("Invalid 'weights' value. Choose 'residuals', 'exp', or 'none'.")
#   }

#   # Weighted least squares regression
#   wls_model <- lm(chi_eta ~ lag_scaled, data = etachispa_df, weights = ws)

#   if (summary) {
#     print(summary(wls_model))
#   }

#   # Extract and rescale parameters
#   intercept <- coef(wls_model)[1]
#   slope <- coef(wls_model)[2]

#   beta <- exp(intercept)
#   alpha <- slope / max_lag  # Rescale to original lag units

#   return(c(intercept = intercept, beta = beta, alpha = alpha))
# }

