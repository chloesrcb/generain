
library(generain)
library(ggplot2)
library(reshape2)
library(animation)
library(RandomFields)
library(RandomFieldsUtils)
source("./script/load_libraries.R")

################################################################################
# Simulate data using the rpareto model
ngrid <- 5
spa <- 1:ngrid
temp <- 0:20
m <- 500 # number of episodes
M <- 1 # number of simulations
n.res <- m * M
param <- c(0.4, 0.2, 1.5, 1) # true parameters for the variogram
beta1 <- param[1]
beta2 <- param[2]
alpha1 <- param[3]
alpha2 <- param[4]
adv <- c(0.1, 0.2)
true_param <- c(beta1, beta2, alpha1, alpha2, adv[1], adv[2])
s0 <- c(1, 1)
t0 <- 0

# Simulate spatio-temporal r-Pareto process
simu_rpar <- sim_rpareto(param[1], param[2], param[3], param[4], spa, spa, temp,
                          adv, s0, t0, n.res)


# Apply formatting
param_str <- format_value(true_param)
adv_str <- format_value(adv)
s0_str <- format_value(s0)
t0_str <- format_value(t0)

# Save the data
foldername <- paste0(data_folder, "simulations/simulations_rpar/rpar_",
                    param_str, "/sim_", ngrid^2, "s_", length(temp), "t_s0_",
                    s0_str, "_t0_", t0_str, "/")

if (!dir.exists(foldername)) {
  dir.create(foldername, recursive = TRUE)
}

save_simulations(simu_rpar, ngrid, n.res,
                 folder = foldername,
                 file = paste0("rpar_", ngrid^2, "s_",
                                length(temp), "t"))



file_path <- paste0(foldername, "rpar_", ngrid^2, "s_",
                    length(temp), "t_1.csv")
simulation_data <- read.csv(file_path)

# create_simu_gif(simulation_data, c(param, adv), type = "rpar", forcedtemp = 30)

################################################################################
# Simulation
files <- list.files(foldername, full.names = TRUE)
length(files)
list_rpar <- list()
for (i in 1:n.res) {
  file_name <- files[i]
  list_rpar[[i]]<- read.csv(file_name)
}


# files <- list.files(foldername, full.names = TRUE)
# length(files)
# list_simuM <- list()
# m <- 500
# nres <- m*M
# for (i in 1:nres) {
#   file_name <- files[i]
#   list_simuM[[i]] <- read.csv(file_name)
# }

# dependece buhl
simu_df <- list_rpar[[1]] # first simulation
nsites <- ncol(simu_df) # number of sites
par(mfrow = c(1, 1))
plot(simu_df[, 1], main = "rpareto simulation")
# get grid coordinates
sites_coords <- generate_grid_coords(sqrt(nsites))

params <- c(param, adv)

# get the lag vectors
df_lags <- get_conditional_lag_vectors(sites_coords, s0, t0,
                                        tau_vect = 0:10)

chi_theorical <- theoretical_chi(params, df_lags, latlon = FALSE,
                                            directional = FALSE)
chi <- unique(chi_theorical$chi)
plot(chi)
tau <- 0:10 # temporal lags
quantile <- 0.8

# get the empirical excesses
excesses <- empirical_excesses(simu_df, quantile, df_lags)
# plot(density(excesses$kij), main = "Excesses")

u = 1
i = 1
# Get the m corresponding simulations from list_simu inside a list
list_episodes <- list_rpar[((i - 1) * m + 1):(i * m)]
# Compute excesses
list_excesses <- lapply(list_episodes, function(episode) {
empirical_excesses_rpar(episode, u, df_lags, threshold = TRUE, t0 = t0)
})

library(parallel)
num_cores <- detectCores() - 1
true_param <- c(0.4, 0.2, 1.5, 1, 0.1, 0.2)
init_params <- c(0.4, 0.2, 1.5, 1, 1.1, 0.8)
result_list <- mclapply(1:M, process_simulation, M = M, m = m,
                        list_simuM = list_rpar, u = u, df_lags = df_lags,
                        t0 = t0, true_param = init_params,
                        wind_df = adv, hmax = sqrt(17),
                        mc.cores = num_cores)
result_list

generate_init_close <- function(true_param) {
  p1 <- abs(true_param[1] + rnorm(1, 0, sd = 0.1))         # > 0
  p2 <- abs(true_param[2] + rnorm(1, 0, sd = 0.1))         # > 0

  p3 <- true_param[3] + rnorm(1, 0, sd = 0.1)
  p4 <- true_param[4] + rnorm(1, 0, sd = 0.1)

  p5 <- true_param[5] + rnorm(1, 0, sd=0.1)
  p6 <- true_param[6] + rnorm(1, 0, sd=0.1)

  return(c(p1, p2, p3, p4, p5, p6))
}
generate_init_close <- function(true_param) {
  p1 <- true_param[1]
  p2 <- true_param[2]
  p3 <- true_param[3]
  p4 <- true_param[4]

  p5 <- true_param[5] + rnorm(1, 0, sd = 0.2)
  p6 <- true_param[6] + rnorm(1, 0, sd = 0.2)

  return(c(p1, p2, p3, p4, p5, p6))
}


set.seed(123)
n_inits <- 10
true_param_eta <- c(0.4, 0.2, 1.5, 1, 1, 1)
init_list <- replicate(n_inits, generate_init_close(true_param_eta), 
                      simplify = FALSE)

results_all <- lapply(init_list, function(params) {
  mclapply(1:M, process_simulation, M = M, m = m,
           list_simuM = list_rpar, u = u, df_lags = df_lags,
           t0 = t0, true_param = params,
           wind_df = adv, hmax = sqrt(17),
           mc.cores = num_cores)
})

length(results_all)

final_params <- do.call(rbind, lapply(results_all, function(res) do.call(rbind, res)))

colnames(final_params) <- paste0("param", 1:6)

# # Moyenne, écart-type
# summary_stats <- as.data.frame(final_params) %>%
#   summarise(across(everything(), list(mean = mean, sd = sd)))

# print(summary_stats)

# # Boxplot + true_param en rouge pointillé
# boxplot(final_params, main = "")
# abline(h = true_param, col = "red", lty = 2)


# df <- data.frame(value = as.vector(final_params),
#                  param = rep(colnames(final_params), each = nrow(final_params)))

# Define true_param (vector of true values for each parameter)
# Example: true_param <- c(a = 1, b = 2, c = 3)

# Convert true_param into a data frame for plotting

# Parameter names in latex
param_names <- c("beta1", "beta2", "alpha1", "alpha2", "adv1", "adv2")
param_str <- c(expression(beta[1]), expression(beta[2]),
               expression(alpha[1]), expression(alpha[2]),
               expression(v[x]), expression(v[y]))
true_df <- data.frame(param = colnames(final_params), value = true_param)

# Create the boxplot
ggplot(df, aes(x = param, y = value)) +
  geom_boxplot() +
  geom_point(data = true_df, aes(x = param, y = value), 
             color = "red", shape = 4, size = 3, stroke = 1.5) +
  scale_x_discrete(labels = param_str) +
  theme_minimal()

folder_save <- paste0(im_folder, "optim/rpar/", param_str, "/")
if (!dir.exists(folder_save)) {
  dir.create(folder_save, recursive = TRUE)
}
filename <- paste0(folder_save, "boxplot_rpar_", n_inits, 
                                      "_different_init.png")
ggsave(filename, width = 10, height = 6)







library(ggplot2)
library(dplyr)
library(tidyr)

# Parameter names in latex
param_names <- c("beta1", "beta2", "alpha1", "alpha2", "adv1", "adv2")
param_str <- c(expression(beta[1]), expression(beta[2]),
               expression(alpha[1]), expression(alpha[2]),
               expression(v[x]), expression(v[y]))

init_params <- do.call(rbind, init_list)
colnames(init_params) <- colnames(final_params)

df_init <- data.frame(value = as.vector(init_params),
                      param = rep(param_names,
                              each = nrow(init_params)),
                      type = "Initial")

df_final <- data.frame(value = as.vector(final_params),
                       param = rep(param_names,
                              each = nrow(final_params)),
                       type = "Final")




df_all <- rbind(df_init, df_final)

# Plot
ggplot(df_all, aes(x = param, y = value, fill = type)) +
  geom_boxplot(position = position_dodge(width = 0.75)) +
  geom_point(data = true_df, aes(x = param_names, y = value),
             inherit.aes = FALSE,
             color = "red", shape = 4, size = 3, stroke = 1.5) +
  scale_fill_manual(values = c("Initial" = "skyblue", "Final" = "orange")) +
  theme_minimal() +
  labs(title = "Boxplot of Parameters: Initial vs Final",
       y = "Parameter Value", x = "Parameter",
       fill = "Type")




df_final[df_final$param == "beta2", ]$value





# rmse_optim <- round(sqrt((result$par - params)^2), 5)
# print(rmse_optim)

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
