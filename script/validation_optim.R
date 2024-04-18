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
library(generain)
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
n_res <- 10 # number of simulations
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
                            tau = 1:10, df_dist = df_dist,
                            parscale = c(0.01, 0.01, 1, 1))

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
excesses <- empirical_excesses(BR_df, 0.9, 1:10, h_vect, df_dist, nmin = 5)
result <- optim(par = c(0.4, 0.2, 1.5, 1), fn = neg_ll, excesses = excesses,
                h_vect = h_vect, tau = 1:10, df_dist = df_dist)

result$par # estimated parameters

# for all simulations
df_result <- evaluate_optim(list_BR, quantile = 0.9, true_param = true_param,
                            tau = 1:10, df_dist = df_dist, method = "CG",
                            nmin = 5, parscale = c(0.01, 0.01, 1, 1))
install# get RMSE, MAE, Mean
df_valid_2 <- get_criterion(df_result, param)

# Save df_result and df_valid in a file
write.table(df_result, file = "results.txt", sep = "\t", row.names = FALSE)
write.table(df_valid, file = "results.txt", sep = "\t", row.names = FALSE, append = TRUE)



# simu with advection ----------------------------------------------------------

# Simulation with 9 sites and 50 times and advection
list_BR_adv_9_50 <- list()
for (i in 1:10) {
  file_path <- paste0("../data/simulations_BR/sim_adv_9s_50t_10/rainBR_", i, ".csv")
  df <- read.csv(file_path)
  list_BR_adv_9_50[[i]] <- df
}

BR_df <- list_BR_adv_9_50[[1]] # first simulation
nsites <- ncol(BR_df) # number of sites
df_dist <- distances_regular_grid(nsites) # distance matrix
h_vect <- get_h_vect(df_dist, sqrt(9)) # spatial lags

# for one simulation
excesses <- empirical_excesses(BR_df, 0.5, 1:10, h_vect, df_dist, nmin = 2)
result <- optim(par = c(0.4, 0.2, 1.5, 1), fn = neg_ll, excesses = excesses,
                h_vect = h_vect, tau = 1:10, df_dist = df_dist,
                method = "CG", control=list(parscale=c(0.1, 0.1, 1, 1)))

# Plot the results -------------------------------------------------------------

# Load the ggplot2 library
library(ggplot2)

# Create the boxplot using ggplot
p <- ggplot(df_result) +
    geom_boxplot(aes(x = factor(1), y = beta1), fill = btfgreen,
                                                color = "black") +
    geom_boxplot(aes(x = factor(2), y = beta2), fill = btfgreen,
                                                    color = "black") +
    geom_boxplot(aes(x = factor(3), y = alpha1), fill = btfgreen,
                                                color = "black") +
    geom_boxplot(aes(x = factor(4), y = alpha2), fill = btfgreen,
                                                color = "black") +
    btf_theme +
    ylab("Estimates") +
    xlab("") +
    scale_x_discrete(labels = c(TeX("$\\widehat{\\beta}_1$"),
                                TeX("$\\widehat{\\beta}_2$"),
                                TeX("$\\widehat{\\alpha}_1$"),
                                TeX("$\\widehat{\\alpha}_2$")))

# Save the plot
phd_extremes_im <- "../../phd_extremes/thesis/resources/images"

filename <- paste0(phd_extremes_im, "/loglike_optim/boxplot_BR_49s_100t.png")
ggsave(plot = p, filename = filename, width = 10, height = 10)
