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

