
library(generain)
library(ggplot2)
library(reshape2)
library(animation)


################################################################################
# Simulate data using the rpareto model
ngrid <- 5
spa <- 1:ngrid
temp <- 1:50
n.res <- 100
param <- c(0.4, 0.4, 1.5, 1) # true parameters for the variogram
beta1 <- param[1]
beta2 <- param[2]
alpha1 <- param[3]
alpha2 <- param[4]
adv <- c(0.05, 0.02)
s0 <- c(1, 1)
t0 <- 1

# Simulate spatio-temporal r-Pareto process
simu_rpar <- sim_rpareto(param[1], param[2], param[3], param[4], spa, spa, temp,
                          adv, s0, t0, n.res)

# Save the simulations
save_simulations(simu_rpar, ngrid, n.res,
                 folder = "./data/simulations_rpar/",
                 file = paste0("rpar_", ngrid^2, "s_",
                                length(temp), "t"))

file_path <- paste0("./data/simulations_rpar/rpar_", ngrid^2, "s_",
                                length(temp), "t_1.csv")
simulation_data <- read.csv(file_path)

create_simu_gif(simulation_data, c(param, adv), type = "rpar", forcedtemp = 30)

################################################################################
# Simulation
list_rpar <- list()
for (i in 1:n.res) {
  file_path <- paste0("./data/simulations_rpar/rpar_", ngrid^2, "s_",
                                length(temp), "t_", i, ".csv")
  df <- read.csv(file_path)
  list_rpar[[i]] <- df
}

# dependece buhl
simu_df <- list_rpar[[1]] # first simulation
nsites <- ncol(simu_df) # number of sites
par(mfrow = c(1, 1))
plot(simu_df[, 1], main = "rpareto simulation")
# get grid coordinates
sites_coords <- generate_grid_coords(sqrt(nsites))

params <- c(param, adv)

# get the lag vectors
df_lags <- get_conditional_lag_vectors(sites_coords, params, s0, t0,
                          hmax = sqrt(17), tau_vect = 0:10)

chi_theorical <- theorical_chi(params, df_lags)
chi <- unique(chi_theorical$chi)
plot(chi)
tau <- 0:10 # temporal lags
quantile <- 0.8

# get the empirical excesses
excesses <- empirical_excesses(simu_df, quantile, df_lags)
# plot(density(excesses$kij), main = "Excesses")


result <- optim(par = c(params), fn = neg_ll,
                  data = simu_df,
                  quantile = quantile,
                  df_lags = df_lags,
                  excesses = excesses,
                  locations = sites_coords,
                  hmax = sqrt(17),
                  s0 = s0,
                  t0 = t0,
                  method = "BFGS",
                  control = list(parscale = c(1, 1, 1, 1, 1, 1),
                                 maxit = 10000))

rmse_optim <- round(sqrt((result$par - params)^2), 5)
print(rmse_optim)

################################################################################
# All simulations together
################################################################################
# get combination of all simulations in list_rpar
simu_all <- do.call(rbind, list_rpar)
dim(simu_all)
# get the empirical excesses
df_lags <- get_conditional_lag_vectors(sites_coords, params, s0, t0,
                          hmax = sqrt(17), tau_vect = 0:10)
quantile <- 0.8
excesses_all <- empirical_excesses(simu_all, quantile, df_lags)
dim(excesses_all)

neg_ll_composite <- function(params, list_simu, df_lags, locations, quantile,
                    list_excesses, latlon = FALSE, s0 = NA, t0 = NA, hmax = NA) {

  nll_composite <- 0
  # number of simulations  in list_simu
  nsim <- length(list_simu)
  for (i in 1:nsim) {
    simu <- list_simu[[i]]
    excesses <- list_excesses[[i]]
    nll_i <- neg_ll(params, simu, df_lags, locations, quantile,
                    latlon = latlon, excesses = excesses, hmax = hmax, s0 = s0,
                    t0 = t0)

    nll_composite <- nll_composite + nll_i
  }

  return(nll_composite)
}

nfiles <- 100

quantile <- 0.7
list_excesses <- list()
for (i in 1:nfiles) {
  excesses <- empirical_excesses(list_rpar[[i]], quantile = quantile,
                                df_lags = df_lags)
  list_excesses[[i]] <- excesses
}

result <- optim(par = c(params), fn = neg_ll_composite,
                  list_simu = list_rpar,
                  df_lags = df_lags,
                  quantile = quantile,
                  list_excesses = list_excesses,
                  locations = sites_coords,
                  hmax = sqrt(17),
                  s0 = s0,
                  t0 = t0,
                  method = "BFGS",
                  control = list(parscale = c(1, 1, 1, 1, 1, 1),
                                 maxit = 10000))


# plot negative log-likelihood fixing every parameter but beta2
neg_ll_beta2 <- function(beta2, beta1, alpha1, alpha2, adv1, adv2) {
  params <- c(beta1, beta2, alpha1, alpha2, adv1, adv2)
  neg_ll(params, simu_all, df_lags, sites_coords, quantile, excesses_all,
         hmax = sqrt(17), s0 = s0, t0 = t0)
}

beta2 <- seq(0.01, 0.5, 0.01)
nll_beta2 <- sapply(beta2, neg_ll_beta2)

plot(beta2, nll_beta2, type = "l", main = "Negative log-likelihood",
     xlab = "beta2", ylab = "Negative log-likelihood")

library(graphics)

# Définir les paramètres fixes
fixed_params <- c(0.4, 0.2, 1.5, 1)

# Créer une grille de valeurs pour deux paramètres (par exemple, paramètre 1 et paramètre 2)
param1_values <- seq(0.00001, 0.9, length.out = 50)
param2_values <- seq(0.00001, 0.9, length.out = 50)

# Matrice pour stocker les résultats de la log-vraisemblance négative
nll_matrix <- matrix(NA, nrow = length(param1_values), ncol = length(param2_values))

# Calculer la log-vraisemblance pour chaque combinaison de param1 et param2
for (i in 1:length(param1_values)) {
  for (j in 1:length(param2_values)) {
    params <- c(param1_values[i], param2_values[j], fixed_params[3:4])
    nll_matrix[i, j] <- neg_ll(params, simu_df, df_lags, locations, quantile, 
                    excesses, hmax = sqrt(17), s0 = s0, t0 = t0)
  }
}

par(mfrow=c(1,1))

png("../images/optim/contour_25s_300t_betas_rpar.png")
contour(param1_values, param2_values, nll_matrix, nlevels = 20,
        xlab = "Beta1", ylab = "Beta2",
        main = "Contour plot of Negative Log-Likelihood")

# Save plot as PNG file
dev.off()

# only beta2 varies
beta2_values <- seq(0.000001, 0.5, length.out = 50)
q <- 0.8
nll <- numeric(length(beta2_values))
excesses <- empirical_excesses(simu, q, df_lags)
for (i in 1:length(beta2_values)) {
  params <- c(fixed_params[1], beta2_values[i], fixed_params[3:4])
  nll[i] <- neg_ll(params, simu_df, df_lags, locations, q, excesses,
                    hmax = sqrt(17), s0 = s0, t0 = t0)
}

plot(beta2_values, nll, type = "l", xlab = "Beta2", ylab = "Negative Log-Likelihood",
     main = "Negative Log-Likelihood as a function of Beta2")


neg_ll_composite <- function(params, data, df_lags, locations, quantile,
                    excesses, latlon = FALSE, hmax = NA, nsample = 1, s0 = NA,
                    t0 = NA) {
  ntemp <- nrow(data) / nsample # number of time steps in each simulation

  nll_composite <- 0 # composite negative log-likelihood
  for (i in 1:nsample) {
    # extract simulation data from i-th simulation
    simu <- data[((i - 1) * ntemp + 1):(i * ntemp), ]

    # excesses <- list_excesses[[i]]
    nll_i <- neg_ll(params, simu, df_lags, locations, quantile,
                    latlon = latlon, hmax = hmax,
                     excesses = excesses, s0 = s0, t0 = t0)

    nll_composite <- nll_composite + nll_i
  }

  return(nll_composite)
}

excesses_all <- empirical_excesses(simu_all, quantile, df_lags)

result <- optim(par = c(params), fn = neg_ll_composite,
                  data = simu_all,
                  quantile = quantile,
                  df_lags = df_lags,
                  excesses = excesses_all,
                  locations = sites_coords,
                  hmax = sqrt(17),
                  s0 = s0,
                  t0 = t0,
                  nsample = n.res,
                  method = "BFGS",
                  control = list(parscale = c(1, 1, 1, 1, 1, 1),
                                 maxit = 10000))


################################################################################
# Verification
################################################################################

library(extRemes)
library(ismev)
library(POT)

rpar <- simulation_data$S1
threshold <- quantile(rpar, probs = 0.98)
rpar_exc <- rpar[rpar > threshold]
fit_gpd <- gpd.fit(rpar, threshold)
sigma <- fit_gpd$mle[1]
xi <- fit_gpd$mle[2]

theorical_qgpd <- qgpd(ppoints(rpar_exc), loc=min(rpar_exc),
                       shape=xi, scale=sigma)

qqplot(rpar_exc, theorical_qgpd, main = "GPD Q-Q plot",
  xlab = "Empirical quantiles",
  ylab = "Theoretical quantiles")
