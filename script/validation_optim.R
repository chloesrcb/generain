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

################################################################################
# functions --------------------------------------------------------------------

boxplot_optim <- function(df_result, true_param, filename) {
  # plot the results
  p <- ggplot(df_result) +
    geom_boxplot(aes(x = factor(1), y = true_param[1]), fill = btfgreen,
                 color = "black") +
    geom_boxplot(aes(x = factor(2), y = true_param[2]), fill = btfgreen,
                 color = "black") +
    geom_boxplot(aes(x = factor(3), y = true_param[3]), fill = btfgreen,
                 color = "black") +
    geom_boxplot(aes(x = factor(4), y = true_param[4]), fill = btfgreen,
                 color = "black") +
    btf_theme +
    ylab("Estimates") +
    xlab("") +
    scale_x_discrete(labels = c(TeX("$\\widehat{\\beta}_1$"),
                                TeX("$\\widehat{\\beta}_2$"),
                                TeX("$\\widehat{\\alpha}_1$"),
                                TeX("$\\widehat{\\alpha}_2$"))) +
    geom_point(aes(x = factor(1), y = true_param[1]), color = "darkred",
                shape = 4, size = 3) +
    geom_text(aes(x = factor(1), y = true_param[1], label = true_param[1]),
              color = "darkred", hjust = -0.8) +
    geom_point(aes(x = factor(2), y = true_param[2]), color = "darkred",
                shape = 4, size = 3) +
    geom_text(aes(x = factor(2), y = true_param[2], label = true_param[2]),
              color = "darkred", hjust = -0.8) +
    geom_point(aes(x = factor(3), y = true_param[3]), color = "darkred",
                shape = 4, size = 3) +
    geom_text(aes(x = factor(3), y = true_param[3], label = true_param[3]),
              color = "darkred", hjust = -0.8) +
    geom_point(aes(x = factor(4), y = true_param[4]), color = "darkred",
                shape = 4, size = 3) +
    geom_text(aes(x = factor(4), y = true_param[4], label = true_param[4]),
            color = "darkred", hjust = -2)

  ggsave(plot = p, filename = filename, width = 10, height = 10)
}

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
for (i in 1:2) {
  file_path <- paste0("../data/simulations_BR/sim_25s_300t/new_rain_br_", i,
                      ".csv")
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

df_valid_1 <- get_criterion(df_result, true_param)


# simu with more sites and less time -------------------------------------------

# true parameters
# beta1, beta2, alpha1, alpha2

true_param <- c(0.4, 0.2, 1.5, 1)

# Iterate through files from 1 to 100
list_BR <- list()
for (i in 1:2) {
  file_path <- paste0("../data/simulations_BR/sim_49s_300t/rain_br_", i, ".csv")
  df <- read.csv(file_path)
  list_BR[[i]] <- df
}

BR_df <- list_BR[[1]] # first simulation
nsites <- ncol(BR_df) # number of sites
df_dist <- distances_regular_grid(nsites) # distance matrix
h_vect <- get_h_vect(df_dist, sqrt(17)) # spatial lags

# # for one simulation
excesses <- empirical_excesses(BR_df, 0.9, 1:10, h_vect, df_dist, nmin = 5)
result <- optim(par = c(0.4, 0.2, 1.5, 1), fn = neg_ll, excesses = excesses,
                h_vect = h_vect, df_dist = df_dist, quantile = 0.9, tau = 1:10,
                nmin = 5, method = "CG",
                control = list(parscale = c(1, 1, 1, 1)))

# result$par # estimated parameters

# for all simulations
start_time <- Sys.time()
df_result <- evaluate_optim(list_BR, quantile = 0.9, true_param = true_param,
                            tau = 1:10, df_dist = df_dist, method = "CG",
                            nmin = 5, parscale = c(0.01, 0.01, 1, 1))
end_time <- Sys.time()
print(end_time - start_time)
# get RMSE, MAE, Mean
df_valid_2 <- get_criterion(df_result, param)

# simu with advection ----------------------------------------------------------

# Simulation with 9 sites and 50 times and advection
list_BR_adv_9_50 <- list()
for (i in 1:10) {
  file_path <- paste0("../data/simulations_BR/sim_adv_9s_50t_10/rainBR_", i,
                      ".csv")
  df <- read.csv(file_path)
  list_BR_adv_9_50[[i]] <- df
}
adv <- c(0.05, 0.05) # advection
tau <- 1:5 # temporal lags
BR_df <- list_BR_adv_9_50[[1]] # first simulation
nsites <- ncol(BR_df) # number of sites
df_dist <- distances_regular_grid(nsites) # distance matrix
h_vect <- get_h_vect(df_dist, sqrt(17)) # spatial lags

# for one simulation
excesses <- empirical_excesses(BR_df, 0.6, tau, h_vect, df_dist, nmin = 5)
result <- optim(par = c(0.4, 0.2, 1.5, 1, adv), fn = neg_ll,
                excesses = excesses, simu = BR_df, quantile = 0.6,
                h_vect = h_vect, tau = tau, df_dist = df_dist,
                method = "CG",
                control = list(parscale = c(0.1, 0.1, 1, 1, 0.01, 0.01)))


# Plot the results -------------------------------------------------------------

beta1 <- 0.4
beta2 <- 0.2
alpha1 <- seq(0, 2, 0.01)
alpha2 <- 1

nll <- numeric(length(alpha1))

for (i in seq_along(nll)) {
  nll[i] <- neg_ll(c(beta1, beta2, alpha1[i], alpha2),
                   excesses = excesses, h_vect = h_vect, tau = tau,
                   df_dist = df_dist, quantile = 0.9)
}

p <- ggplot() +
  geom_point(data = data.frame(alpha1 = alpha1, nll = nll),
             aes(x = alpha1, y = nll)) +
  xlab(TeX("$\\alpha_1$")) +
  ylab("Negative Log-Likelihood") +
  ylim(10000, 25000) +
  btf_theme 

p

nll_beta_alpha <- matrix(NA, nrow = length(beta1), ncol = length(alpha1))

for (i in seq_along(beta1)) {
  for (j in seq_along(alpha1)) {
    nll_beta_alpha[i, j] <- neg_ll(c(beta1[i], beta2, alpha1[j], alpha2),
                                excesses = excesses, h_vect = h_vect, tau = tau,
                                df_dist = df_dist, quantile = 0.9)
  }
}

library(reshape2)
nll_data <- melt(nll_beta_alpha)

p <- ggplot(nll_data) +
  geom_point(aes(x = X1/100, y = value, color = X2/10)) +
  ylim(25000, 150000) +
  labs(title = "Negative log-likelihood as a function of spatial parameters") +
  btf_theme +
  xlab(TeX("$\\beta_1$")) +
  ylab("Negative Log-Likelihood") +
  scale_color_continuous(name = TeX("$\\alpha_1$"))

p
# save plot
file_path <- paste0("../../phd_extremes/thesis/resources/images/loglike_optim/",
                    "nll_04_02_alpha1_1_25s_300t.png")
ggsave(file_path, plot = p, width = 7, height = 7, dpi = 300)


# Plot the results -------------------------------------------------------------

# Load the ggplot2 library
library(ggplot2)
df_result <- df_result_400
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

filename <- paste0(phd_extremes_im, "/validation/boxplot_depmodel_400s_50t.png")
ggsave(plot = p, filename = filename, width = 10, height = 10)


################################################################################
# new br process
################################################################################
file_path <- paste0("./data/simulations_BR/sim_25s_300t/rain_br_1.csv")
simulation_data <- read.csv(file_path)
ngrid <- sqrt(ncol(simulation_data))
adv <- c(0, 0) # advection
tau <- 1:10 # temporal lags max
BR_df <-  simulation_data # first simulation
nsites <- ncol(BR_df) # number of sites
df_dist <- distances_regular_grid(nsites) # distance matrix
h_vect <- get_h_vect(df_dist, sqrt(17)) # spatial lags

# for one simulation
excesses <- empirical_excesses(BR_df, 0.9, tau, h_vect, df_dist, nmin = 5)
result <- optim(par = c(0.4, 0.2, 1.5, 1), fn = neg_ll,
                excesses = excesses, simu = BR_df, quantile = 0.9,
                h_vect = h_vect, tau = tau, df_dist = df_dist,
                method = "CG",
                control = list(parscale = c(0.1, 0.1, 1, 1)))


result <- optim(par = c(0.4, 0.2, 1.5, 1), fn = neg_ll,
                excesses = excesses, simu = BR_df, quantile = 0.9,
                h_vect = h_vect, tau = tau, df_dist = df_dist,
                method = "L-BFGS-B", lower = c(0, 0, 1e-6, 1e-6),
                upper = c(Inf, Inf, 2, 2),
                control = list(parscale = c(0.1, 0.1, 1, 1)))
