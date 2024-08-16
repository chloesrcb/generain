

#### libraries ####
# remove objects from the environment
rm(list = ls())

# clear console
cat("\014")

# change working directory
setwd("./script")

# load libraries
source("load_libraries.R")
library(generain)
library(optimx)
library(bbmle)

# Simulation
num_iterations <- 5
ngrid <- 5
temp <- 1:300
list_BR <- list()
for (i in 1:num_iterations) {
  file_path <- paste0("../data/simulations_BR/sim_25s_300t/br_",
                      ngrid^2, "s_", length(temp), "t_", i, ".csv")
  df <- read.csv(file_path)
  list_BR[[i]] <- df
}

true_param <- c(0.4, 0.2, 1.5, 1)


simu_df <- list_BR[[1]] # first simulation
nsites <- ncol(simu_df) # number of sites
plot(simu_df[, 8])

# get grid coordinates
sites_coords <- generate_grid_coords(sqrt(nsites))
dist_mat <- get_dist_mat(sites_coords,
                         latlon = FALSE) # distance matrix
df_dist <- reshape_distances(dist_mat) # reshape the distance matrix

nsites <- ncol(simu_df) # number of sites
# get grid coordinates
sites_coords <- generate_grid_coords(sqrt(nsites))
# sites_coords$id <- 1:nsites
h_vect <- get_lag_vectors(sites_coords, true_param,
                          hmax = sqrt(17), tau_vect = 0:10)


BR_df <- list_BR[[1]]


q <- 0.96
excesses <- empirical_excesses(BR_df, q, h_vect)
excesses_new <- excesses[excesses$n_vect > 0, ]
chi <- theorical_chi(true_param, excesses_new)
chi$chiemp <- chi$n_vect / chi$N_vect
head(chi)

result <- optim(par = c(true_param), fn = neg_ll,
                        simu = BR_df,
                        quantile = q,
                        excesses = excesses,
                        h_vect = h_vect, tau = 0:10,
                        locations = sites_coords,
                        method = "CG",
                        control = list(parscale = c(1, 1, 1, 1),
                                        maxit = 10000))


neg_ll2 <- function(beta1, beta2, alpha1, alpha2, simu, h_vect, tau, locations,
                  latlon = FALSE, quantile = 0.9, nmin = 5,
                  simu_exp = FALSE, excesses = NULL, adv1 = 0, adv2 = 0) {
  params <- c(beta1, beta2, alpha1, alpha2)
  # hmax <- max(h_vect$hnorm)

  print(params)
  if (is.null(excesses)) {
    excesses <- empirical_excesses(simu, quantile, tau, h_vect,
                                  nmin)
  }

  lower.bound <- c(1e-6, 1e-6, 1e-6, 1e-6)
  upper.bound <- c(Inf, Inf, 1.999, 1.999)

  # Check if the parameters are in the bounds
  if (any(params < lower.bound) || any(params > upper.bound)) {
    print("out of bounds")
    return(1e8)
  }

  N_vect <- excesses$N_vect # number of observations
  n_vect <- excesses$n_vect # number of excesses
  chi <- theorical_chi(params, h_vect) # get chi matrix
  # transform in chi vector
  chi_vect <- as.vector(chi$chi)
  chi_vect <- ifelse(chi_vect <= 0, 0.000001, chi_vect)

  non_excesses <- N_vect - n_vect # number of non-excesses
  # log-likelihood vector
  ll_vect <- n_vect * log(chi_vect) + non_excesses * log(1 - chi_vect)

  # negative log-likelihood
  nll <- -sum(ll_vect, na.rm = TRUE)
  return(nll)
}


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

      excess_count <- 0
      num_cond_excesses <- 0

      for (tk in 1:min(Tmax - 2 * t, 10)) {
        tl <- tk + t
        if (rain_cp$s1[tk] > q && rain_cp$s2[tl] > q) {
          excess_count <- excess_count + 1
        }
        if (rain_cp$s2[tl] > q) {
          num_cond_excesses <- num_cond_excesses + 1
        }
      }

      # store the number of excesses
      N_vect <- c(n_vect, num_cond_excesses)
      n_vect <- c(N_vect, excess_count)
    }
  }
  excesses <- list(n_vect = n_vect, N_vect = N_vect)
  return(excesses)
}


q <- 0.7
excesses <- empirical_excesses(BR_df, q, h_vect)

plot(excesses$n_vect)

res <- mle2(neg_ll2, start = list(beta1 = true_param[1],
                                 beta2 = true_param[2],
                                 alpha1 = true_param[3],
                                 alpha2 = true_param[4]),
                 data = list(simu = BR_df,
                        quantile = q,
                        excesses = excesses,
                        h_vect = h_vect, tau = 0:10,
                        locations = sites_coords),
                  control = list(maxit = 10000),
                  fixed = list())


# empirical_excesses_st <- function(data_rain, quantile, h_vect) {
#   q <- quantile # quantile
#   excesses <- h_vect
#   unique_tau <- unique(h_vect$tau)
#   unique_h <- unique(h_vect$hnorm)

#   for (h in unique_h){
#     # station number inside h lag
#     ind_h <- which(h_vect$hnorm == h)
#     indices <- h_vect[ind_h, c("s1", "s2")]
#     print(paste0("h = ", h))
#     nb_pairs <- length(indices$s1)
#     # get index pairs
#     ind_s1 <- indices$s1
#     ind_s2 <- indices$s2
#     chi_val_h <- c()
#     for (i in 1:nb_pairs){
#         rain_cp <- drop_na(data_rain[, c(ind_s1[i], ind_s2[i])])
#         colnames(rain_cp) <- c("s1", "s2")
#         Tmax <- nrow(rain_cp)
#         for (j in 1:nmax){
#             rain_nolag <- rain_cp$s1[j:(Tmax - t)]
#             rain_lag <- rain_cp$s2[(j + t):Tmax]
#             data <- cbind(rains1, rains2) # get couple
#             chi_val_h <- c(chi_val_h, get_chiq(data, q))
#         }
#       }
#       chi_slag <- c(chi_slag, mean(na.omit(chi_val_h)))
#     }
#   }
#   return(excesses)
# }


empirical_excesses <- function(data_rain, quantile, h_vect) {
    q <- quantile # quantile
    n_vect <- c()
    N_vect <- c()
    unique_tau <- unique(h_vect$tau)

    for (t in unique_tau) {
        df_h_t <- h_vect[h_vect$tau == t, ]

        for (i in seq_len(nrow(df_h_t))) {
            ind_s2 <- as.numeric(as.character(df_h_t$s2[i]))
            ind_s1 <- df_h_t$s1[i]

            rain_cp <- data_rain[, c(ind_s1, ind_s2), drop = FALSE]
            rain_cp <- na.omit(rain_cp)
            colnames(rain_cp) <- c("s1", "s2")
            tmax <- nrow(rain_cp)
            for (j in 1:50){ # loop over time for each lag
                rain_nolag <- rain_cp$s1[j:(tmax - t)]
                rain_lag <- rain_cp$s2[(j + t):tmax]
                n <- length(rain_nolag)
                rain_unif <- cbind(rank(rain_nolag) / (n + 1),
                                   rank(rain_lag) / (n + 1))
                cp_cond <- rain_unif[rain_unif[, 2] > q, , drop = FALSE]

                excess_count <- sum(cp_cond[, 1] > q)
                num_cond_excesses <- nrow(cp_cond)
                N_vect <- c(N_vect, num_cond_excesses)
                n_vect <- c(n_vect, excess_count)
            }
        }
    }
  return(h_vect)
}






empirical_excesses <- function(data_rain, quantile, h_vect) {
  Tmax <- nrow(data_rain) # number of time steps
  q <- quantile # quantile

  unique_tau <- unique(h_vect$tau)

  for (t in unique_tau) {
    df_h_t <- h_vect[h_vect$tau == t, ]

    for (i in seq_len(nrow(df_h_t))) {
      ind_s2 <- as.numeric(as.character(df_h_t$s2[i]))
      ind_s1 <- df_h_t$s1[i]

      rain_cp <- data_rain[, c(ind_s1, ind_s2), drop = FALSE]
      rain_cp <- na.omit(rain_cp)
      colnames(rain_cp) <- c("s1", "s2")

      rain_nolag <- rain_cp$s1[1:(Tmax - t)]
      rain_lag <- rain_cp$s2[(1 + t):Tmax]

      n <- length(rain_nolag)
      rain_unif <- cbind(rank(rain_nolag) / (n + 1), rank(rain_lag) / (n + 1))
      cp_cond <- rain_unif[rain_unif[, 2] > q, , drop = FALSE]

      excess_count <- sum(cp_cond[, 1] > q)
      num_cond_excesses <- nrow(cp_cond)
      h_vect$N_vect[h_vect$s1 == ind_s1 & h_vect$s2 == ind_s2 & h_vect$tau == t] <- num_cond_excesses
      h_vect$n_vect[h_vect$s1 == ind_s1 & h_vect$s2 == ind_s2 & h_vect$tau == t] <- excess_count
    }
  }
  return(h_vect)
}




empirical_excesses <- function(data_rain, quantile, h_vect) {
  q <- quantile # quantile
  unique_tau <- unique(h_vect$tau) # unique temporal lags
  
  # Initialise une liste pour stocker les résultats
  results_list <- list()
  
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
      
      for (tk in 1:10) {
        tl <- tk + t
        rain_nolag <- rain_cp$s1[tk:(Tmax - t)]
        rain_lag <- rain_cp$s2[tl:Tmax]

        n <- length(rain_nolag)
        rain_unif <- cbind(rank(rain_nolag) / (n + 1), rank(rain_lag) / (n + 1))
        cp_cond <- rain_unif[rain_unif[, 2] > q, , drop = FALSE]

        excess_joint <- sum(cp_cond[, 1] > q)
        excess_marginal <- nrow(cp_cond)

        # Store the results in the list
        results_list <- append(results_list, list(data.frame(
          s1 = ind_s1,
          s2 = ind_s2,
          tau = t,
          tk = tk,
          excess_joint = excess_joint,
          excess_marginal = excess_marginal
        )))
      }
    }
  }
  
  # Combine the list into a single data frame
  results_df <- do.call(rbind, results_list)
  results_df <- h_vect %>%
    left_join(results_df, by = c("s1", "s2", "tau"))
  return(results_df)
}

q <- 0.9
excesses <- empirical_excesses(BR_df, q, h_vect)
excesses_new <- excesses[excesses$excess_joint > 0, ]
chi <- theorical_chi(true_param, excesses_new)
chi$chiemp <- chi$excess_joint / chi$excess_marginal
head(chi)
# plot density
plot(density(excesses$excess_joint))
plot(density(excesses$n_vect))
# Add a binomial distribution
n <- length(excesses$n_vect)
p <- mean(excesses$n_vect) / n
x <- seq(0, n, by = 1)
y <- dbinom(x, size = n, prob = p)
lines(x, y, col = "red")



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

theorical_chi <- function(params, excesses) {
  chi_df <- excesses
  # tau <- unique(h_vect$tau)
  # h_vectors <- unique(h_vect$hnorm)

  chi_df$chi <- theorical_chi_ind(params, excesses$hnorm, excesses$tau)

  return(chi_df)
}
-
neg_ll2 <- function(beta1, beta2, alpha1, alpha2, simu, h_vect, tau, locations,
                  latlon = FALSE, quantile = 0.9, nmin = 5,
                  simu_exp = FALSE, excesses = NULL, adv1 = 0, adv2 = 0) {
  params <- c(beta1, beta2, alpha1, alpha2)
  # hmax <- max(h_vect$hnorm)

  print(params)
  if (is.null(excesses)) {
    excesses <- empirical_excesses(simu, quantile, tau, h_vect,
                                  nmin)
  }

  lower.bound <- c(1e-6, 1e-6, 1e-6, 1e-6)
  upper.bound <- c(Inf, Inf, 1.999, 1.999)

  # Check if the parameters are in the bounds
  if (any(params < lower.bound) || any(params > upper.bound)) {
    print("out of bounds")
    return(1e8)
  }

  N_vect <- excesses$excess_marginal # number of observations
  n_vect <- excesses$excess_joint # number of excesses
  chi <- theorical_chi(params, excesses) # get chi matrix
  # transform in chi vector
  chi_vect <- as.vector(chi$chi)
  chi_vect <- ifelse(chi_vect <= 0, 0.000001, chi_vect)

  non_excesses <- N_vect - n_vect # number of non-excesses
  # log-likelihood vector
  ll_vect <- n_vect * log(chi_vect) + non_excesses * log(1 - chi_vect)

  # negative log-likelihood
  nll <- -sum(ll_vect, na.rm = TRUE)
  return(nll)
}

q <- 0.9
BR_df <- list_BR[[1]]
excesses <- empirical_excesses(BR_df, q, h_vect)

res <- mle2(neg_ll2, start = list(beta1 = true_param[1],
                                 beta2 = true_param[2],
                                 alpha1 = true_param[3],
                                 alpha2 = true_param[4]),
                 data = list(simu = BR_df,
                        quantile = q,
                        excesses = excesses,
                        h_vect = h_vect, tau = 1:10,
                        locations = sites_coords),
                  control = list(maxit = 10000),
                  fixed = list(alpha2 = 1))


### 12 august 2024


get_lag_vectors <- function(df_coords, params, tau_vect, hmax = NA) {
  alpha_spa <- params[3]
#   if (length(params) != 6) {
#     adv <- c(0, 0)
#   } else {
#     adv <- params[5:6]
#   }

  n <- nrow(df_coords)
  total_combinations <- (n * (n - 1) / 2) * length(tau_vect)

  # Preallocate the dataframe
  lags <- data.frame(
    s1 = integer(total_combinations),
    s2 = integer(total_combinations),
    h1 = numeric(total_combinations),
    h2 = numeric(total_combinations),
    hnorm = numeric(total_combinations),
    tau = numeric(total_combinations)
  )

  idx <- 1
  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      lag_latitude <- df_coords$Latitude[j] - df_coords$Latitude[i]
      lag_longitude <- df_coords$Longitude[j] - df_coords$Longitude[i]

      hnorm <- norm_Lp(lag_latitude, lag_longitude, alpha_spa)

      if (is.na(hmax) || hnorm <= hmax) {
        for (tau in tau_vect) {
          lags$s1[idx] <- i
          lags$s2[idx] <- j
          lags$h1[idx] <- lag_latitude
          lags$h2[idx] <- lag_longitude
          lags$tau[idx] <- tau
          lags$hnorm[idx] <- hnorm
          idx <- idx + 1
        }
      }
    }
  }

  # Remove the unused preallocated space
  lags <- lags[1:(idx - 1), ]

  return(lags)
}

true_param <- c(0.4, 0.2, 1.5, 1)
simu_df <- list_BR[[1]] # first simulation
nsites <- ncol(simu_df) # number of sites
plot(simu_df[, 1])

# get grid coordinates
sites_coords <- generate_grid_coords(sqrt(nsites))
tau_vect <- 0:10
df_lags <- get_lag_vectors(sites_coords, true_param, tau = tau_vect,
                          hmax = sqrt(17))


count_excesses <- function(data_rain, quantile, df_lags) {
    q <- quantile # quantile
    Tmax <- nrow(data_rain) # number of time observations
    # nsites <- ncol(data_rain)
    excesses <- df_lags
    excesses$k_ij <- c()
    excesses$n_j <- c()
    for (i in seq_len(nrow(df_lags))) {
        t <- df_lags$tau[i]
        ind_s2 <- as.numeric(as.character(df_lags$s2[i]))
        ind_s1 <- df_lags$s1[i]
        # get the data for the pair of sites
        rain_cp <- data_rain[, c(ind_s1, ind_s2), drop = FALSE]
        rain_cp <- na.omit(rain_cp)
        colnames(rain_cp) <- c("s1", "s2")

        rain_nolag <- rain_cp$s1[1:(Tmax - t)]
        rain_lag <- rain_cp$s2[(1 + t):Tmax]
        # transform in uniform distribution
        n <- length(rain_nolag)
        rain_unif <- cbind(rank(rain_nolag) / (n + 1), rank(rain_lag) / (n + 1))
        cp_cond <- rain_unif[rain_unif[, 2] > q, , drop = FALSE]

        excesses_joint <- sum(cp_cond[, 1] > q)
        excesses_marginal <- nrow(cp_cond)

        excesses$k_ij[i] <- excesses_joint
        excesses$n_j[i] <- excesses_marginal
    }
    return(excesses)
}


nll <- function(params, simu, df_lags, locations,
                  latlon = FALSE, quantile = 0.9, excesses = NULL) {
  print(params)
  if (is.null(excesses)) {
    excesses <- empirical_excesses(simu, quantile, df_lags)
  }

  lower.bound <- c(1e-6, 1e-6)
  upper.bound <- c(Inf, 1.999)

  # Check if the parameters are in the bounds
  if (any(params < lower.bound) || any(params > upper.bound)) {
    print("out of bounds")
    return(1e8)
  }

  n_j <- excesses$n_j # number of observations
  k_ij <- excesses$k_ij # number of excesses
  chi <- theorical_chi(params, excesses) # get chi matrix
  # transform in chi vector
  chi_vect <- as.vector(chi$chi)
  chi_vect <- ifelse(chi_vect <= 0, 0.000001, chi_vect)

  non_excesses <- n_j - k_ij # number of non-excesses
  # log-likelihood vector
  ll_vect <- k_ij * log(chi_vect) + non_excesses * log(1 - chi_vect)

  # negative log-likelihood
  nll <- -sum(ll_vect, na.rm = TRUE)
  return(nll)
}




excesses <- count_excesses(BR_df, 0.9, df_lags)

result <- optim(par = true_param, fn = nll,
                        simu = BR_df,
                        quantile = 0.9,
                        excesses = excesses,
                        df_lags = df_lags,
                        locations = sites_coords,
                        method = "CG",
                        control = list(parscale = c(1, 1, 1, 1),
                                        maxit = 10000))



nll <- function(beta1, beta2, alpha1, alpha2, simu, df_lags, locations,
                  latlon = FALSE, quantile = 0.9, excesses = NULL) {
  params <- c(beta1, beta2, alpha1, alpha2)
  print(params)
  if (is.null(excesses)) {
    excesses <- empirical_excesses(simu, quantile, df_lags)
  }

  lower.bound <- c(1e-6, 1e-6, 1e-6, 1e-6)
  upper.bound <- c(Inf, Inf, 1.999, 1.999)

  # Check if the parameters are in the bounds
  if (any(params < lower.bound) || any(params > upper.bound)) {
    print("out of bounds")
    return(1e8)
  }

  n_j <- excesses$n_j # number of observations
  k_ij <- excesses$k_ij # number of excesses
  chi <- theorical_chi(params, excesses) # get chi matrix
  # transform in chi vector
  chi_vect <- as.vector(chi$chi)
  chi_vect <- ifelse(chi_vect <= 0, 0.000001, chi_vect)

  non_excesses <- n_j - k_ij # number of non-excesses
  # log-likelihood vector
  ll_vect <- k_ij * log(chi_vect) + non_excesses * log(1 - chi_vect)

  # negative log-likelihood
  nll <- -sum(ll_vect, na.rm = TRUE)
  return(nll)
}

res <- mle2(nll, start = list(beta1 = true_param[1],
                                 beta2 = true_param[2],
                                 alpha1 = true_param[3],
                                 alpha2 = true_param[4]),
                 data = list(simu = BR_df,
                        quantile = 0.9,
                        excesses = excesses,
                        df_lags = df_lags,
                        locations = sites_coords),
                  control = list(maxit = 10000),
                  fixed = list(beta1 = 0.4, alpha1 = 1.5))

## 13 aout

# Définir le quantile q
q <- 0.9

data <- BR_df

library(reshape2)


# Function to calculate excess indicators for each site
excesses_indicators <- function(data, q) {
  # Assume each row represents a moment in time
  data$time <- 1:nrow(data)

  # Convert data to long format
  data_long <- melt(data, id.vars = "time", variable.name = "site",
                    value.name = "value")
  colnames(data_long) <- c("time", "site", "value")

  # Calculate threshold for each site
  thresholds <- aggregate(value ~ site, data_long, function(x) quantile(x, q))

  # Add a column for excesses
  data_long <- merge(data_long, thresholds, by = "site",
                      suffixes = c("", "_threshold"))
  data_long$excess <- ifelse(data_long$value > data_long$value_threshold, 1, 0)

  return(data_long)
}

# Function to calculate marginal excesses for each site
calculate_marginal_excesses <- function(data_long) {
  nj <- aggregate(excess ~ site, data_long, sum)
  return(nj)
}

# Function to calculate k_ij,tau for a pair of sites
calculate_kij_tau <- function(data, site_i, site_j, tau) {
  data_i <- data[data$site == site_i, ]
  data_j <- data[data$site == site_j, ]

  if (tau >= 0) {
    data_i_n <- data_i[1:(nrow(data_i) - tau), ]
    data_j_n <- data_j[(tau + 1):nrow(data_j), ]
  } else {
    data_i_n <- data_i[(-tau + 1):nrow(data_i), ]
    data_j_n <- data_j[1:(nrow(data_j) + tau), ]
  }

  sum(data_i_n$excess & data_j_n$excess)
}

# Exemple d'utilisation pour un décalage temporel tau
site_i <- "S1"
site_j <- "S4"
tau <- 1

kij_tau <- calculate_kij_tau(data_long, site_i, site_j, tau)
print(kij_tau)

# Function to calculate k_ij,tau for different tau values
calculate_kij_tau_list <- function(data_long, taus) {
  kij_tau_list <- list()
  sites <- unique(data_long$site)

  for (tau in taus) {
    # Initialize matrix
    kij_matrix <- matrix(0, nrow = length(sites), ncol = length(sites),
                         dimnames = list(sites, sites))
    for (i in 1:length(sites)) { # for each pair of sites
      for (j in 1:length(sites)) {
        if (i != j) {
          kij_matrix[i, j] <- calculate_kij_tau(data_long, sites[i], sites[j],
                                                tau)
        }
      }
    }
    kij_tau_list[[as.character(tau)]] <- kij_matrix # for each tau
  }

  return(kij_tau_list)
}

# Function to calculate gamma_theta
gamma_theta <- function(h, tau, theta) {
  # Example of modified Brown-Resnick variogram depending on h and tau
  beta1 <- theta[1]
  beta2 <- theta[2]
  alpha1 <- theta[3]
  alpha2 <- theta[4]
  gamma <- 2 * (beta1 * h^alpha1 + beta2 * abs(tau)^alpha2)

  if (gamma < 0) {
    gamma <- 0
  }
  return(gamma)
}

# Function to calculate chi_ij with distance matrix
chi_ij <- function(site_i, site_j, tau, theta, distance_matrix) {
  h <- distance_matrix[site_i, site_j]
  gamma <- gamma_theta(h, tau, theta)
  return(2 * (1 - pnorm(sqrt(0.5 * gamma))))
}

# Function to calculate composite log-likelihood
log_likelihood <- function(theta, data, distance_matrix, taus, kij_tau_list, nj) {
  print(theta)
  ll <- 0
  sites <- unique(data$site)

  lower.bound <- c(1e-6, 1e-6, 1e-6, 1e-6)
  upper.bound <- c(Inf, Inf, 1.999, 1.999)

  # Check if the parameters are in the bounds
  if (any(theta < lower.bound) || any(theta > upper.bound)) {
    print("out of bounds")
    return(1e8)
  }

  for (tau in taus) {
    kij_matrix <- kij_tau_list[[as.character(tau)]]
    for (i in 1:length(sites)) {
      for (j in 1:length(sites)) {
        if (i != j) {
          site_i <- sites[i]
          site_j <- sites[j]
          nij <- nj$excess[nj$site == site_j]
          kij <- kij_matrix[site_i, site_j]

          chi <- chi_ij(site_i, site_j, tau, theta, distance_matrix)

          ll <- ll + kij * log(chi) + (nij - kij) * log(1 - chi)
        }
      }
    }
  }
    
  return(-ll)  # Negative log likelihood
}

distance_matrix <- as.matrix(dist(sites_coords))
rownames(distance_matrix) <- colnames(distance_matrix) <- paste0("S", 
                                                          1:nrow(sites_coords))

print(distance_matrix)

q <- 0.9
data_long <- excesses_indicators(BR_df, q)

# Initialize theta parameters
theta_init <- c(0.4, 0.2, 1.5, 1)
taus <- 1:10
nj <- aggregate(excess ~ site, data_long, sum) # marginal excesses
kij_tau_list <- calculate_kij_tau_list(data_long, taus) # conditional excesses
# Optimize composite log-likelihood
optim_result <- optim(par = theta_init,
                      fn  = log_likelihood,
                      data = data_long,
                      distance_matrix = distance_matrix,
                      taus = taus,
                      kij_tau_list = kij_tau_list,
                      nj = nj,
                      method = "BFGS",
                      control = list(maxit = 5000))

theta_opt <- optim_result$par

# Display results
print(theta_opt)

optim_result <- optim(
  par = theta_init,
  fn = log_likelihood,
  data = data_long,
  distance_matrix = distance_matrix,
  kij_tau_list = kij_tau_list,
  taus = taus,
  nj = nj,
  method = "L-BFGS-B",
  lower = c(0, 0, 0, 0),
  upper = c(Inf, Inf, 1.99, 1.99),
  control = list(maxit = 5000)
)

library(ggplot2)

# Calculate chi_ij for different distances and time lags
calculate_chi_values <- function(distance_matrix, taus, theta) {
  chi_values <- list()
  
  for (tau in taus) {
    chi_matrix <- matrix(0, nrow = nrow(distance_matrix), ncol = ncol(distance_matrix))
    rownames(chi_matrix) <- rownames(distance_matrix)
    colnames(chi_matrix) <- colnames(distance_matrix)
    
    for (i in 1:nrow(distance_matrix)) {
      for (j in 1:ncol(distance_matrix)) {
        if (i != j) {
          chi_matrix[i, j] <- chi_ij(rownames(distance_matrix)[i], colnames(distance_matrix)[j], tau, theta, distance_matrix)
        }
      }
    }
    
    chi_values[[as.character(tau)]] <- chi_matrix
  }
  
  return(chi_values)
}

# Calculate chi_ij for optimized parameters
taus <- 1:20
chi_values <- calculate_chi_values(distance_matrix, taus, theta_init)

# Convert results to data frame for visualization
chi_df <- do.call(rbind, lapply(names(chi_values), function(tau) {
  chi_matrix <- chi_values[[tau]]
  data.frame(
    site_i = rep(rownames(chi_matrix), each = ncol(chi_matrix)),
    site_j = rep(colnames(chi_matrix), nrow(chi_matrix)),
    chi = as.vector(chi_matrix),
    tau = as.numeric(tau)
  )
}))

# Visualize results
ggplot(chi_df, aes(x = tau, y = chi, color = interaction(site_i, site_j))) +
  geom_line() +
  labs(title = "Variation of chi_ij with tau", x = "Time Lag (tau)", y = "chi_ij") +
  theme_minimal() +
  theme(legend.position = "none")


# Calculer les probabilités empiriques p_ij et les valeurs théoriques chi_ij
calculate_probabilities <- function(data, distance_matrix, taus, theta, 
                                    kij_tau_list, nj) {
  results <- list()
  sites <- unique(data$site)
  nsites <- length(sites)
  for (tau in taus) {
    kij_matrix <- kij_tau_list[[as.character(tau)]]

    for (i in 1:(nsites - 1)) {
      for (j in (i + 1):nsites) {
        site_i <- sites[i]
        site_j <- sites[j]
        kij <- kij_matrix[site_i, site_j]
        nij <- nj$excess[nj$site == site_j]
        p_ij <- kij / nij
        chi_ij_value <- chi_ij(site_i, site_j, tau, theta, distance_matrix)
        results[[paste(site_i, site_j, tau, sep = "_")]] <- list(
          p_ij = p_ij,
          chi_ij = chi_ij_value
        )
      }
    }
  }
  
  return(results)
}

q <- 0.9
data_long <- excesses_indicators(BR_df, q)

# Initialize theta parameters
theta_init <- c(0.4, 0.2, 1.5, 1)
taus <- 1:10
nj <- aggregate(excess ~ site, data_long, sum) # marginal excesses
kij_tau_list <- calculate_kij_tau_list(data_long, taus) # conditional excesses

# Calculer les probabilités empiriques et théoriques
results <- calculate_probabilities(data_long, distance_matrix, taus, theta_init, 
            kij_tau_list, nj)

# Convertir les résultats en data frame pour comparaison
results_df <- do.call(rbind, lapply(names(results), function(key) {
  data.frame(
    pair_tau = key,
    p_ij = results[[key]]$p_ij,
    chi_ij = results[[key]]$chi_ij
  )
}))

# Visualiser la comparaison entre p_ij et chi_ij
ggplot(results_df, aes(x = p_ij, y = chi_ij)) +
  geom_point() +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
  labs(x = "Empirical probabilities p_ij",
       y = "Theorical probabilities chi_ij") +
  theme_minimal()



# Function to calculate marginal excesses for each site and each time lag
calculate_marginal_excesses <- function(data_long, taus, q) {
  results <- list()

  for (tau in taus) {
    shifted_data <- data_long

    # Adjust the time index for the given tau
    shifted_data$time <- shifted_data$time + tau

    # Remove any rows where the adjusted time index exceeds the original time range
    shifted_data <- shifted_data[shifted_data$time <= max(data_long$time), ]

    # Calculate the excesses for the shifted data
    shifted_data$excess <- shifted_data$value > q

    # Aggregate to calculate the sum of excesses for each site
    nj_tau <- aggregate(excess ~ site, shifted_data, sum)

    # Store the results in a list
    results[[as.character(tau)]] <- nj_tau
  }
  
  return(results)
}

calculate_marginal_excesses <- function(data_long, taus, q) {
  results <- list()
  
  for (tau in taus) {
    shifted_data <- data_long
    
    # Create a new column for shifted values
    shifted_data$value_shifted <- NA
    
    # Shift the values for the given tau
    for (site in unique(data_long$site)) {
      site_data <- data_long[data_long$site == site, ]
      if (tau == 0) {
        shifted_data$value_shifted[shifted_data$site == site] <- site_data$value
      } else {
        shifted_data$value_shifted[shifted_data$site == site] <- c(rep(NA, tau), site_data$value[1:(nrow(site_data) - tau)])
      }
    }
    
    # Calculate the excesses for the shifted values
    shifted_data$excess <- shifted_data$value_shifted > q
    
    # Aggregate to calculate the sum of excesses for each site
    nj_tau <- aggregate(excess ~ site, shifted_data, sum, na.rm = TRUE)
    
    # Store the results in a list
    results[[as.character(tau)]] <- nj_tau
  }
  
  return(results)
}

# Exemple d'utilisation
# data_long doit être un data frame avec des colonnes: 'time', 'site', 'value'
data_long <- data.frame(
  time = rep(1:nrow(data), times = ncol(data)),
  site = rep(colnames(data), each = nrow(data)),
  value = as.vector(t(data))
)

# Liste des décalages temporels (taus)
taus <- 0:10

# Calculer les excès marginaux pour chaque site et chaque décalage temporel
nj_tau_list <- calculate_marginal_excesses(data_long, taus, q = 0.95)  # Remplacez 0.95 par le quantile approprié
