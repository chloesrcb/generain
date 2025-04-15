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


################################################################################
# functions --------------------------------------------------------------------

boxplot_optim <- function(df_result, true_param, filename) {
  # plot the results
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

parscale.parameters <- function(par, scale, fix = 1) {
  # check if length of scale is equal to par
  if (length(par) !=  length(scale)) {
    stop("parscale.parameters has parameter and scaling vectors of different sizes.")
  }

  if (any(scale == 0)) {
    scale[scale == 0] <- 1
  }
  # parscale fixes the larged par/parscale value to deviate only 10 percent,
  # others can then vary get fixed and maximal value
  fix.value <- par[fix]
  if(fix.value == 0){
    fix.value <- 1
  }

  max.value <- abs(par[which.max(abs(par))[1]])
  if (max.value == 0){
    print("Warning: Vector contains only zeroes. Scaling is set to a default of 1.")
    max.value <- 1
  }

  #fill in scaling vector
  par.scale <- scale / (0.1 * fix.value * max.value)
  par.scale[fix] <- 1 / max.value

  return(abs(par.scale))
}

################################################################################
# Test on exponential ----------------------------------------------------------
################################################################################

# Simulate expo ~ chi ----------------------------------------------------------

# create a distance matrix
nsites <- 25
df_dist <- distances_regular_grid(nsites) # distance matrix
sites_coords <- generate_grid_coords(sqrt(nsites)) # grid coordinates
npairs <- nrow(df_dist) # number of pairs
true_param <- c(0.4, 0.2, 1.5, 1) # true parameters
tau_vect <- 1:10 # temporal lags
nconfig <- npairs * length(tau_vect) # number of configurations
Tmax <- 100 # number of time steps

h_vect <- get_lag_vectors(sites_coords, true_param,
                          hmax = sqrt(17), tau_vect = 1:10) # lag vectors

# chi <- theorical_chi_mat(c(0.4, 0.2, 1.5, 1), h_vect, tau_vect) # chi matrix
chi <- theoretical_chi(true_param, h_vect) # chi matrix

# simulate expo and optimize
n_res <- 100 # number of simulations
df_result <- evaluate_optim_simuExp(n_res, Tmax, tau_vect, h_vect, chi,
                                    sites_coords, nconfig)

df_valid_exp <- get_criterion(df_result, true_param) # get RMSE, MAE, Mean

################################################################################
# BR SIMU ----------------------------------------------------------------------
################################################################################

# simu with 25 sites and 300 times ---------------------------------------------

# true parameters
# beta1, beta2, alpha1, alpha2
true_param <- c(0.4, 0.2, 1.5, 1)

nsites <- 25
ntimes <- 300

file_name <- paste0("../data/simulations_BR/sim_", nsites,
                    "s_", ntimes, "t/br_", nsites, "s_", ntimes, "t_")

# Iterate through files from 1 to 100
list_BR <- list()
for (i in 1:5) {
  file_path <- paste0(
          file_name, i, ".csv")
  df <- read.csv(file_path)
  list_BR[[i]] <- df
}


BR_df <- list_BR[[1]] # first simulation
plot(BR_df$S2)
nsites <- ncol(BR_df) # number of sites

sites_coords <- generate_grid_coords(sqrt(nsites))
# dist_mat <- get_dist_mat(sites_coords, adv = adv, tau = 1:10,
#                          latlon = FALSE) # distance matrix
# df_dist <- reshape_distances(dist_mat) # reshape the distance matrix
# h_vect <- get_h_vect(df_dist, sqrt(17)) # spatial lags

h_vect <- get_lag_vectors(sites_coords, true_param,
                          hmax = sqrt(17), tau_vect = 1:10) # lag vectors

start_time <- Sys.time()
# for all simulations
df_result <- evaluate_optim(list_BR, quantile = 0.9, true_param = true_param,
                            tau_vect = 1:10, hmax = sqrt(17), nmin = 5,
                            locations = sites_coords,
                            parscale = c(1, 1, 1, 1))
end_time <- Sys.time()
print(end_time - start_time)

q <- 0.9
BR_df <- list_BR[[1]]
excesses <- empirical_excesses(BR_df, q, h_vect)

excesses_filtered <- excesses$n_vect[excesses$n_vect > 0]
density_plot <- ggplot(data.frame(x = excesses_filtered), aes(x)) +
  geom_density() +
  labs(
    title = "Density plot of the number of excesses",
    x = "Number of excesses",
    y = "Density"
  ) +
  theme_minimal()

# Display the plot
print(density_plot)

# Save the plot to a file
ggsave("../images/optim/density_plot_25s_100t_90.png", plot = density_plot, width = 8,
       height = 6)



beta1 <- 0.4
beta2 <- 0.2
alpha1 <- seq(0.2, 1.99, 0.01)
alpha2 <- 1


q <- 0.8
BR_df <- list_BR[[1]]
excesses <- empirical_excesses(BR_df, q, h_vect)

nll <- c()
for (i in seq_along(alpha1)) {
  nll[i] <- neg_ll(c(beta1, beta2, alpha1[i], alpha2), simu = BR_df,
                   excesses = excesses, h_vect = h_vect, tau = tau,
                   locations = sites_coords, quantile = q)
}
# which.min(nll)
# alpha1[which.min(nll)]

par(mfrow = c(1, 1))
png("../images/optim/nll_25s_100t_alpha1.png", width = 800, height = 600)
plot(alpha1, nll, type = "l")
dev.off()

true_param <- c(0.4, 0.2, 1.5, 1)
result <- optim(par = c(true_param), fn = neg_ll,
                        simu = BR_df,
                        quantile = q,
                        excesses = excesses,
                        h_vect = h_vect, tau = tau_vect,
                        locations = sites_coords,
                        method = "CG",
                        control = list(parscale = c(1, 1, 1, 1),
                                        maxit = 10000))



library(bbmle)
q <- 0.75
excesses <- empirical_excesses(BR_df, q, 1:10, h_vect, nmin = 5)

res <- mle2(neg_ll, start = list(beta1 = true_param[1],
                                 beta2 = true_param[2],
                                 alpha1 = true_param[3],
                                 alpha2 = true_param[4]),
                 data = list(simu = BR_df,
                        quantile = q,
                        excesses = excesses,
                        h_vect = h_vect, tau = 1:10,
                        locations = sites_coords),
                  control = list(maxit = 10000),
                  fixed = list(beta1 = true_param[1]))

# (fit1 <- mle2(LL, method="L-BFGS-B", lower=c(ymax=0, xhalf=0)))
p1 <- profile(res)

plot(p1, absVal=FALSE)







df_res <- na.omit(df_result)
df_valid_1 <- get_criterion(df_res, true_param)

grad_test <- gHgenb(par = c(0.4, 0.2, 1.5, 1),
          fn = function(par) {
            neg_ll(par, excesses = excesses, simu = BR_df, quantile = 0.9,
                    h_vect = h_vect, tau = tau, df_dist = df_dist)
            }, lower = c(1e-6, 1e-6, 1e-6, 1e-6),
                upper = c(Inf, Inf, 1.999, 1.999))



excesses <- empirical_excesses(BR_df, 0.9, 1:10, h_vect, df_dist, nmin = 5)
# conjugate gradient algorithm with the Dai and Yuan update (2001)

# scaling <- c(0.4, 0.2, 1.5, 1)
# p.scale <- parscale.parameters(true_param, scaling)

param <- c(0.5, 0.1, 1.4, 1.1)
result <- optimr(par = list(beta1 = true_param[1],
                                 beta2 = true_param[2],
                                 alpha1 = true_param[3],
                                 alpha2 = true_param[4]),
                  method = "CG", 
          fn = function(par) {
            neg_ll(par, excesses = excesses, quantile = 0.96,
                    locations = sites_coords,
                    h_vect = h_vect, tau = tau, simu = BR_df)
            }, 
                control = list(parscale = c(1, 1, 0.1, 0.1),
                                maxit = 1000))


# compute the Hessian
optimHess(result$par, fn = function(par) {
            neg_ll(par, excesses = excesses, quantile = 0.9,
                    h_vect = h_vect, tau = tau, df_dist = df_dist)
            })


# Plot the results -------------------------------------------------------------

beta1 <- 0.4
beta2 <- 0.2
alpha1 <- seq(0, 2, 0.01)
alpha2 <- 0.8

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
                    "boxplot_50_optimr_49s_300t.png")
ggsave(file_path, plot = p, width = 7, height = 7, dpi = 300)


boxplot_optim(df_res, true_param, filename = file_path)

# Plot the results -------------------------------------------------------------

# Load the ggplot2 library
library(ggplot2)
df_result <- df_res
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
                control = list(parscale = c(0.1, 0.1, 1.5, 1)))


result <- optim(par = c(0.4, 0.2, 1.5, 1), fn = neg_ll,
                excesses = excesses, simu = BR_df, quantile = 0.9,
                h_vect = h_vect, tau = tau, df_dist = df_dist,
                method = "L-BFGS-B", lower = c(0, 0, 1e-6, 1e-6),
                upper = c(Inf, Inf, 2, 2),
                control = list(parscale = c(0.1, 0.1, 1, 1)))

param <- c(0.4, 0.2, 1.5, 1)
result <- optimr(
  par = true_param,
  method = "Rcgmin",
  gr = "grfwd",
  fn = function(par) {
    neg_ll(par, simu = simu_df, quantile = quantile, h_vect = h_vect,
           tau = tau, locations = sites_coords, df_dist = df_dist,
           excesses = excesses)
  },
  control = list(parscale = c(1, 1, 1, 1),
                 maxit = 10000,
                 trace = 1
                 )
)


################################################################################
# simu rpareto -----------------------------------------------------------------
################################################################################

# true parameters
ngrid <- 5
temp <- 1:300
param <- c(0.2, 0.2, 1.5, 0.8) # true parameters for the variogram
beta1 <- param[1] / 2
beta2 <- param[2] / 2
alpha1 <- param[3]
alpha2 <- param[4]
adv <- c(0.1, 0.1)

file_path <- paste0("../data/simulations_rpar/rpar_", ngrid^2, "s_",
                                length(temp), "t_1.csv")
simu_df <- read.csv(file_path)

params <- c(0.4, 0.2, 1.5, 1, 0.1, 0.1)
nsites <- ncol(simu_df) # number of sites
# get grid coordinates
sites_coords <- generate_grid_coords(sqrt(nsites))
# sites_coords$id <- 1:nsites
h_vect <- get_lag_vectors(sites_coords, params,
                          hmax = sqrt(17), tau_vect = 1:10) # lag vectors

tau <- 1:10 # temporal lags
quantile <- 0.9
nmin <- 5

simu_df <- list_BR[[1]] # first simulation
# get the empirical excesses
excesses <- empirical_excesses(simu_df, quantile, tau, h_vect,
                              nmin)

parscale <- c(0.1, 0.1, 1, 1, 0.1, 0.1)  # scale parameters
lower.bound <- c(1e-6, 1e-6, 1e-6, 1e-6, -Inf, -Inf)
upper.bound <- c(Inf, Inf, 1.999, 1.999, Inf, Inf)

result <- optimr(par = true_param, method = "CG", fn = function(par) {
                  neg_ll(par, simu = simu_df, quantile = quantile,
                        h_vect = h_vect, tau = tau,
                        locations = sites_coords)
                  }, 
                  control = list(parscale = parscale,
                                 maxit = 10000))

# Conjugate gradient method
result <- optim(par = true_param, fn = neg_ll, quantile = 0.9,
                h_vect = h_vect, tau = tau,
                locations = sites_coords,
                simu = simu_df, method = "CG",
                control = list(parscale = parscale))

# Fix all parameters except adv2
param = params
fixed_params <- c(beta1 = param[1], beta2 = param[2], alpha1 = param[3],
                alpha2 = param[4], adv1 = param[5], adv2 = param[6])

# Function to plot neg_ll with fixed parameters
values <- seq(0.05, 1, 0.05)
neg_ll_vector <- numeric(0)
for (i in values) {
  fixed_params[2] <- i
  neg_ll_values <- neg_ll(fixed_params, quantile = quantile, h_vect = h_vect,
                          tau = tau, locations = sites_coords, simu = simu_df)
  neg_ll_vector <- c(neg_ll_vector, neg_ll_values)
}
plot(values, neg_ll_vector, type = "l", xlab = "adv2", ylab = "neg_ll")


beta1 <- 0.4
beta2 <- 0.2
alpha1 <- seq(1, 1.95, 0.05)
alpha2 <- 0.8
adv1 <- 0.1
adv2 <- 0.1

nll <- numeric(length(alpha1))

for (i in seq_along(nll)) {
  nll[i] <- neg_ll(c(beta1, beta2, alpha1[i], alpha2, adv1, adv2),
                   h_vect = h_vect, tau = tau,
                   locations = sites_coords, simu = simu_df, quantile = 0.9)
}


p <- ggplot() +
  geom_point(data = data.frame(alpha1 = alpha1, nll = nll),
             aes(x = alpha1, y = nll)) +
  xlab(TeX("$\\alpha_1$")) +
  ylab("Negative Log-Likelihood") +
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