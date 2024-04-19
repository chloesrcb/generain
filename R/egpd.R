#' Function to initialize values for EGPD fitting
#'
#' This function is used to initialize values for the EGPD fitting with
#' the estimates of the classical GPD parameters.
#' It performs the necessary setup and initialization steps given a classical
#' GPD distribution with a threshold equals to u. For EGPD, classical threshold
#' will be u=0.
#'
#' @param y The data to be fitted.
#' @param u The threshold value.
#'
#' @return GPD estimates for sigma and xi.
#'
#' @import mev
#'
#' @export
init_values <- function(y, u) {
  gpfit <- gp.fit(y, u) # Fit the GPD
  # Get the estimates of the GPD parameters
  sigma_0 <- gpfit$estimate[1]
  xi_0 <- gpfit$estimate[2]
  return(c(sigma_0, xi_0))
}

#' get_egpd_estimates function
#'
#' This function calculates the estimates of the Extreme Generalized Pareto
#' Distribution (EGPD) parameters.
#'
#' @param rain_df A data frame containing the rainfall data.
#' @param left_censoring The threshold for left censoring. Default is 0.
#'
#' @return A list containing the estimates of the EGPD parameters.
#'
#' @export
get_egpd_estimates <- function(rain_df, left_censoring = 0) {
  kappa <- c()
  sigma <- c()
  xi <- c()
  for (col in 1:ncol(rain_df)) {
    y <- as.data.frame(na.omit(rain_df[, col]))
    y <- y[y > 0]
    kappa_0 <- 2
    inits <- init_values(y, 0)
    sigma_0 <- inits[1]
    xi_0 <- inits[2]
    if (length(left_censoring) != 1) {
      censore <- left_censoring[col]
    } else {
      censore <- left_censoring
    }
    egpd.fit <- fit.extgp(y, model = 1, method = "mle",
                          init = c(kappa_0, sigma_0, xi_0),
                          censoring = c(censore, Inf), plots = FALSE,
                          confint = FALSE, ncpus = 7, R = 1000)
    param <- egpd.fit$fit$mle
    kappa <- c(kappa, param[1])
    sigma <- c(sigma, param[2])
    xi <- c(xi, param[3])
  }

  return(list(kappa = kappa, sigma = sigma, xi = xi))
}


#' get_df_long_params_egpd function
#'
#' This function takes a dataframe of parameters and returns a long format 
#' dataframe for EGPD estimates.
#'
#' @param df_params A dataframe containing parameters.
#'
#' @import dplyr
#' @import tidyr
#'
#' @return A long format dataframe.
#'
#' @export
get_df_long_params_egpd <- function(df_params) {
  df_estimates <- data.frame(Xi = df_params$xi,
                             Sigma = df_params$sigma,
                             Kappa = df_params$kappa)

  colnames(df_estimates) <- c("Xi", "Sigma", "Kappa")
  df_long <- gather(df_estimates, key = "Variable", value = "Value")
  return(df_long)
}


#' dgpdExt1 function
#'
#' This function calculates the density distribution of the EGPD with the first 
#' model $G(v) = v^\kappa$ based on the given parameters.
#'
#' @param x The input value.
#' @param kappa The kappa parameter (power transform).
#' @param sigma The sigma parameter (scale).
#' @param gamma The gamma parameter (shape).
#'
#' @return The calculated dgpdExt1 value.
#' 
#' @import POT
#'
#' @export
dgpdExt1 <- function(x, kappa, sigma, gamma){
  h <- dgpd(x / sigma, loc=0, scale=sigma, shape=gamma)
  H <- pgpd(x / sigma, loc=0, scale=sigma, shape=gamma)
  dens <- (kappa / sigma) * h * H^(kappa - 1)
  return(dens)
}


#' choose_censore function
#'
#' This function selects the censored data from a rain_df dataframe for each
#' site according to the NRMSE and RMSE.
#'
#' @param rain_df The input dataframe containing the rain data.
#' @param censore The censored vector.
#' @param nb_simu The number of simulations (default is 100).
#'
#' @return The selected censored data in a datafram according to the NRMSE and
#'         RMSE. # TODO: add CRPS
#'
#' @importFrom mev qextgp
#' @importFrom mev rextgp
#' @importFrom scoringRules crps_sample
#' @importFrom scoringRules crps
#'
#'
#' @export
choose_censore <- function(rain_df, censore, nb_simu = 100) {
  df_score <- data.frame(locations = seq_along(rain_df))
  df_score$RMSE <- Inf
  df_score$NRMSE <- Inf
  NRMSEs <- c()
  # df_score$CRPS <- Inf
  df_score$censoreRMSE <- NA
  df_score$censoreNRMSE <- NA
  # df_score$censoreCRPS <- NA
  for (c in seq_along(censore)) {
    params <- as.data.frame(get_egpd_estimates(rain_df,
                            left_censoring = censore[c]))
    # RMSEs <- c()
    for (i in seq_along(rain_df)){
      y <- na.omit(rain_df[, i])
      y <- y[y > censore[c]]
      # Quantile
      qextgp <- qextgp(p = c(1:length(y))/(length(y) + 1), type = 1,
                    kappa = params$kappa[i], sigma = params$sigma[i],
                    xi = params$xi[i])
      # sort values
      y.sort <- sort(y)
      n <- length(y)
      # RMSE
      RMSE <- sqrt(mean((y.sort - qextgp)^2))
      # RMSEs <- c(RMSEs, RMSE)
      if (RMSE <= df_score$RMSE[i]) {
        df_score$RMSE[i] <- RMSE
        df_score$censoreRMSE[i] <- censore[c]
      }

      NRMSE <- RMSE / mean(y)
      # NRMSEs <- c(NRMSEs, NRMSE)
      if (NRMSE <= df_score$NRMSE[i]) {
        df_score$NRMSE[i] <- NRMSE
        df_score$censoreNRMSE[i] <- censore[c]
      }

    #   CRPSs <- rep(NA, nb_simu)
    #   for (j in 1:nb_simu) {
    #     # simulation with param estimators
    #     X.mle <- rextgp(n, kappa = params$kappa[i], sigma = params$sigma[i],
    #                 xi = params$xi[i])

    #     # CRPS
    #     dat <- t(matrix(rep(X.mle, n), ncol=n))
    #     crps_val <- crps_sample(y, dat)
    #     CRPSs[j] <- mean(crps_val)
    #   }
    #   # get mean of CRPSs
    #   mean_CRPS <- mean(CRPSs)
    #   print(mean_CRPS)
    #   print(censore[c])
    #   if (mean_CRPS <= df_score$CRPS[i]) {
    #      df_score$CRPS[i] <- mean_CRPS
    #     df_score$censoreCRPS[i] <- censore[c]
    #   }
    }
  }
  return(df_score)
}
