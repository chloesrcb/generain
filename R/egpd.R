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
  for (col in seq_len(ncol(rain_df))) {
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
                          censoring = c(censore, Inf), plots = F,
                          confint = F, ncpus = 7, R = 1000)
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
#' model G(v) = v^kappa based on the given parameters.
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
dgpdExt1 <- function(x, kappa, sigma, gamma) {
  h <- dgpd(x / sigma, loc = 0, scale = sigma, shape = gamma)
  H <- pgpd(x / sigma, loc = 0, scale = sigma, shape = gamma)
  dens <- (kappa / sigma) * h * H^(kappa - 1)
  return(dens)
}


#' choose_censore function
#'
#' This function selects the censored data from a rain_df dataframe for each
#' site according to the NRMSE and RMSE.
#'
#' @param rain_df The input dataframe containing the rain data.
#' @param censores The censored vector.
#' @param n_samples The number of samples to generate. Default is 100.
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
choose_censore <- function(rain_df, censores, n_samples = 100) {
  df_score <- data.frame(locations = seq_along(rain_df))
  df_score$RMSE <- Inf
  df_score$censoreRMSE <- NA

  for (c in seq_along(censores)) {
    censore <- censores[c]
    params <- as.data.frame(get_egpd_estimates(rain_df, left_censoring = censore))

    prev_rmse <- Inf
    
    for (i in seq_along(rain_df)) {
      y <- na.omit(rain_df[, i])
      y <- y[y > 0]
      y.sort <- sort(y)
      print(i)
      # Quantiles
      p <- (1:length(y)) / (length(y) + 1)
      qextgp <- qextgp(p = p, type = 1, kappa = params$kappa[i], sigma = params$sigma[i], xi = params$xi[i])

      RMSE <- sqrt(mean((y.sort - qextgp)^2))
      print(RMSE)
      # if (RMSE > prev_rmse & RMSE != Inf) {
      #   next
      # }

      if (RMSE < df_score$RMSE[i]) {
        df_score$RMSE[i] <- RMSE
        df_score$censoreRMSE[i] <- censore
      }

      prev_rmse <- RMSE
    }
  }
  return(df_score)
}


# Get left censoring
compute_rmse <- function(y, left_censore) {
  inits <- init_values(y, 0)
  sigma_0 <- inits[1]
  xi_0 <- inits[2]
  kappa_0 <- 1

  fit <- fit.extgp(y, model = 1, method = "mle",
                   init = c(kappa_0, sigma_0, xi_0),
                   censoring = c(left_censore, Inf),
                   plots = FALSE, confint = FALSE, ncpus = 1)

  param <- fit$fit$mle

  probs <- ppoints(length(y))
  q_theo <- qextgp(probs, type = 1, kappa = param[1], sigma = param[2], xi = param[3])
  q_emp <- sort(y)

  rmse <- sqrt(mean((q_emp - q_theo)^2))
  return(rmse)
}



boot_fun <- function(data, indices) {
  y_boot <- data[indices]
  inits_boot <- init_values(y_boot, 0)
  fit_boot <- tryCatch({
    fit <- fit.extgp(y_boot, model = 1, method = "mle",
                     init = c(kappa_0, inits_boot[1], inits_boot[2]),
                     censoring = c(0.22, Inf), confint = FALSE, plots = FALSE)
    return(as.numeric(fit$fit$mle))
  }, error = function(e) return(rep(NA, 3)))

  return(fit_boot)
}


compute_rmse <- function(y, left_censore) {
  inits <- init_values(y, 0)
  sigma_0 <- inits[1]
  xi_0 <- inits[2]
  kappa_0 <- 1

  fit <- tryCatch({
    fit.extgp(y, model = 1, method = "mle",
              init = c(kappa_0, sigma_0, xi_0),
              censoring = c(left_censore, Inf),
              plots = FALSE, confint = FALSE, ncpus = 1)
  }, error = function(e) return(NULL))

  if (is.null(fit)) return(Inf)

  param <- fit$fit$mle
  probs <- ppoints(length(y))
  q_theo <- qextgp(probs, type = 1, kappa = param[1], sigma = param[2], xi = param[3])
  q_emp <- sort(y)

  sqrt(mean((q_emp - q_theo)^2))
}



process_site <- function(y, site_name, save_path) {
  y <- y[y > 0]
  if (length(y) < 20) return(NULL)  # trop peu de données

  censures <- seq(0.22, 0.25, by = 0.01)
  rmse_vals <- sapply(censures, function(c) compute_rmse(y, c))
  best_cens <- censures[which.min(rmse_vals)]
  # best_cens <- 0.22
  # Fit EGPD avec meilleure censure
  inits <- init_values(y, 0)
  fit <- fit.extgp(y, model = 1, method = "mle",
                   init = c(1, inits[1], inits[2]),
                   censoring = c(best_cens, Inf),
                   plots = FALSE, confint = FALSE, ncpus = 7, R = 1000)

  param_mle <- fit$fit$mle
  probs <- ppoints(length(y))
  q_mle <- qextgp(probs, type = 1, kappa = param_mle[1], sigma = param_mle[2], xi = param_mle[3])
  y_sorted <- sort(y)

  # Bootstrap for CI
  boot_fun <- function(data, indices) {
    y_boot <- data[indices]
    inits_boot <- init_values(y_boot, 0)
    fit_boot <- tryCatch({
      fit <- fit.extgp(y_boot, model = 1, method = "mle",
                      init = c(1, inits_boot[1], inits_boot[2]),
                      censoring = c(best_cens, Inf), confint = FALSE, 
                      plots = FALSE)
      return(as.numeric(fit$fit$mle))
    }, error = function(e) return(rep(NA, 3)))

    return(fit_boot)
  }

  boot_res <- boot(data = y, statistic = boot_fun, R = R)

  # Compute quantiles for bootstrap samples
  set.seed(123)
  probs <- c(1:length(y)) / (length(y) + 1)

  # Only keep valid rows (no NA/Inf/NaN)
  valid_rows <- apply(boot_res$t, 1, function(x) all(is.finite(x)))
  boot_t_valid <- boot_res$t[valid_rows, ]

  q_mle_boot <- apply(boot_t_valid, 1, function(params) {
    qextgp(p = probs, type = 1, kappa = params[1], sigma = params[2], 
            xi = params[3])
  })
  q_mle_boot <- t(q_mle_boot)

  q_lower <- apply(q_mle_boot, 2, quantile, probs = 0.025, na.rm = TRUE)
  q_upper <- apply(q_mle_boot, 2, quantile, probs = 0.975, na.rm = TRUE)

  # Graphique
  df_plot <- data.frame(empirical = y_sorted, fitted = q_mle, lower = q_lower, 
                        upper = q_upper)
  p <- ggplot(df_plot, aes(x = empirical, y = fitted)) +
    geom_ribbon(aes(ymin = lower, ymax = upper), fill = "#69b3a2", alpha = 0.2) +
    geom_point(color = "#69b3a2", size = 1.2) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "#a33737") +
    xlab("Empirical quantiles") + ylab("Fitted quantiles") +
    theme_minimal() +
    theme(
        panel.grid.major = element_line(color = "#d6d1d1"),
        panel.grid.minor = element_line(color = "#d6d1d1")
    )

  filename <- paste0(save_path, "qqplot_egpd_", site_name, "_lcensoring_", best_cens, ".png")
  ggsave(filename = filename, p, width = 6, height = 6, dpi = 400, device = "png", bg = "transparent")

  return(data.frame(Site = site_name, BestCens = best_cens, RMSE = round(min(rmse_vals), 5)))
}


fit_egpd_site <- function(y, site_name) {
  y <- y[y > 0]
  if (length(y) < 20) return(NULL) 

  censures <- seq(0.22, 0.26, by = 0.01)
  rmse_vals <- sapply(censures, function(c) compute_rmse(y, c))
  best_cens <- censures[which.min(rmse_vals)]

  inits <- init_values(y, 0)
  fit <- fit.extgp(y, model = 1, method = "mle",
                   init = c(1, inits[1], inits[2]),
                   censoring = c(best_cens, Inf),
                   plots = FALSE, confint = FALSE, ncpus = 7, R = 1000)

  param_mle <- fit$fit$mle
  kappa <- param_mle[1]
  sigma <- param_mle[2]
  xi <- param_mle[3]

  return(data.frame(Site = site_name, BestCens = best_cens, RMSE = round(min(rmse_vals), 5),
                    kappa = kappa, sigma = sigma, xi = xi))
}

