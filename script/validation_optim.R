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
# Test each functions ----------------------------------------------------------
################################################################################

nsites <- 25
# npairs <- choose(nsites, 2) + nsites # or n(n+1)/2
df_dist <- distances_regular_grid(nsites)
npairs <- nrow(df_dist)

# Simulate expo ~ chi ----------------------------------------------------------
h_vect <- get_h_vect(df_dist, sqrt(17))
tau_vect <- 1:10
nconfig <- npairs * length(tau_vect)
Tmax <- 100
chi <- theorical_chi_mat(c(0.4, 0.2, 1.5, 1), h_vect, tau_vect)
n_res <- 10

df_result <- evaluate_optim_simuExp(n_res, Tmax, tau_vect, h_vect, chi, df_dist,
                                    nconfig)

true_param <- c(0.4, 0.2, 1.5, 1)

df_valid <- get_criterion(df_result, true_param)


################################################################################
# BR SIMU ----------------------------------------------------------------------
################################################################################

# true parameters
# beta1, beta2, alpha1, alpha2
param <- c(0.4, 0.2, 1.5, 1)

list_BR <- list()

# Iterate through files from 1 to 100
for (i in 1:100) {
  file_path <- paste0("../data/simulations_BR/sim_25s_300t/rainBR_", i, ".csv")
  df <- read.csv(file_path)
  list_BR[[i]] <- df
}

BR_df <- list_BR[[1]]
nsites <- ncol(BR_df)
df_dist <- distances_regular_grid(nsites)

# compute the excess ind matrix and the sum of the excesses on simulation
data_rain <- BR_df

h_vect <- get_h_vect(df_dist, sqrt(17))

excesses <- empirical_excesses(BR_df, 0.9, 1:10, h_vect, df_dist)

result <- optim(par = c(0.4, 0.2, 1.5, 1), fn = neg_ll, excesses = excesses,
                h_vect = h_vect, tau = 1:10, df_dist = df_dist,
                method = "CG")
result$par


df_result <- evaluate_optim(list_BR, quantile = 0.9, true_param = param,
                            tau = 1:10, df_dist = df_dist)

# simu 2

list_BR <- list()

# Iterate through files from 1 to 100
for (i in 1:100) {
  file_path <- paste0("../data/simulations_BR/sim_49s_100t/rainBR_", i, ".csv")
  df <- read.csv(file_path)
  list_BR[[i]] <- df
}

BR_df <- list_BR[[1]]
nsites <- ncol(BR_df)
df_dist <- distances_regular_grid(nsites)

# compute the excess ind matrix and the sum of the excesses on simulation
data_rain <- BR_df

h_vect <- get_h_vect(df_dist, sqrt(17))

excesses <- empirical_excesses(BR_df, 0.8, 1:10, h_vect, df_dist)

result <- optim(par = c(0.4, 0.2, 1.5, 1), fn = neg_ll, excesses = excesses,
                h_vect = h_vect, tau = 1:10, df_dist = df_dist,
                method = "CG")
result$par


df_result <- evaluate_optim(list_BR, quantile = 0.9, true_param = param,
                            tau = 1:10, df_dist = df_dist)
