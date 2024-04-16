#### libraries ####
# remove objects from the environment
rm(list = ls())

# clear console
cat("\014")

# change working directory
setwd("./script")

# load libraries
source("load_libraries.R")
# load global functions
source("../R/utils.R")
source("../R/distances.R")
source("../R/optimization.R")

################################################################################
# Test on exponential ----------------------------------------------------------
################################################################################

# Simulate expo ~ chi ----------------------------------------------------------

# create a distance matrix
nsites <- 25
df_dist <- distances_regular_grid(nsites) # distance matrix
npairs <- nrow(df_dist) # number of pairs
h_vect <- get_h_vect(df_dist, sqrt(17)) # spatial lags
tau_vect <- 1:10 # temporal lags
nconfig <- npairs * length(tau_vect) # number of configurations
Tmax <- 100 # number of time steps

chi <- theorical_chi_mat(c(0.4, 0.2, 1.5, 1), h_vect, tau_vect) # chi matrix

# simulate expo and optimize
n_res <- 100 # number of simulations
df_result <- evaluate_optim_simuExp(n_res, Tmax, tau_vect, h_vect, chi, df_dist,
                                    nconfig)

true_param <- c(0.4, 0.2, 1.5, 1) # true parameters
df_valid_exp <- get_criterion(df_result, true_param) # get RMSE, MAE, Mean


################################################################################
# BR SIMU ----------------------------------------------------------------------
################################################################################

# simu with 25 sites and 300 times ---------------------------------------------

# true parameters
# beta1, beta2, alpha1, alpha2
true_param <- c(0.4, 0.2, 1.5, 1)

# Iterate through files from 1 to 100
list_BR <- list()
for (i in 1:100) {
  file_path <- paste0("../data/simulations_BR/sim_25s_300t/rainBR_", i, ".csv")
  df <- read.csv(file_path)
  list_BR[[i]] <- df
}

BR_df <- list_BR[[1]] # first simulation
nsites <- ncol(BR_df) # number of sites
df_dist <- distances_regular_grid(nsites) # distance matrix
h_vect <- get_h_vect(df_dist, sqrt(17)) # spatial lags

df_result <- evaluate_optim(list_BR, quantile = 0.9, true_param = true_param,
                            tau = 1:10, df_dist = df_dist)

df_valid_1 <- get_criterion(df_result, param)

# simu with more sites and less time -------------------------------------------

# true parameters
# beta1, beta2, alpha1, alpha2
true_param <- c(0.4, 0.2, 1.5, 1)

# Iterate through files from 1 to 100
list_BR <- list()
for (i in 1:100) {
  file_path <- paste0("../data/simulations_BR/sim_49s_100t/rainBR_", i, ".csv")
  df <- read.csv(file_path)
  list_BR[[i]] <- df
}

BR_df <- list_BR[[1]] # first simulation
nsites <- ncol(BR_df) # number of sites
df_dist <- distances_regular_grid(nsites) # distance matrix
h_vect <- get_h_vect(df_dist, sqrt(17)) # spatial lags

# for one simulation
excesses <- empirical_excesses(BR_df, 0.8, 1:10, h_vect, df_dist, nmin = 10)
result <- optim(par = c(0.4, 0.2, 1.5, 1), fn = neg_ll, excesses = excesses,
                h_vect = h_vect, tau = 1:10, df_dist = df_dist,
                method = "CG")
result$par # estimated parameters

# for all simulations
df_result <- evaluate_optim(list_BR, quantile = 0.8, true_param = true_param,
                            tau = 1:10, df_dist = df_dist, method = "CG",
                            nmin = 10)
# get RMSE, MAE, Mean
df_valid_2 <- get_criterion(df_result, param)
