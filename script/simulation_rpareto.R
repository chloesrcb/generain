
library(generain)
library(ggplot2)
library(reshape2)
library(animation)
library(RandomFields)
library(RandomFieldsUtils)
source("./script/load_libraries.R")

################################################################################
# Simulate data using the rpareto model
ngrid <- 7
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
adv <- c(1, 2)
true_param <- c(beta1, beta2, alpha1, alpha2, adv[1], adv[2])
s0 <- c(1, 1)
t0 <- 0
random_s0 <- FALSE

# Simulate spatio-temporal r-Pareto process
simu_rpar <- sim_rpareto(param[1], param[2], param[3], param[4], spa, spa, temp,
                          adv, t0, n.res, random_s0, s0)


# Apply formatting
param_str <- format_value(true_param)
adv_str <- format_value(adv)
s0_str <- format_value(s0)
t0_str <- format_value(t0)

# Save the data
if (random_s0) {
  s0_str <- "random_s0"
} else {
  s0_str <- paste0("s0_", s0_str)
}
foldername <- paste0(data_folder, "simulations/simulations_rpar/rpar_",
                    param_str, "/sim_", ngrid^2, "s_", length(temp), "t_",
                    s0_str, "_t0_", t0_str, "/")

if (!dir.exists(foldername)) {
  print("Creating folder")
  dir.create(foldername, recursive = TRUE)
}

save_simulations(simu_rpar, ngrid, n.res,
                 folder = foldername,
                 file = paste0("rpar_", ngrid^2, "s_", length(temp), "t"))



file_path <- paste0(foldername, "rpar_", ngrid^2, "s_",
                    length(temp), "t_1.csv")
simulation_data <- read.csv(file_path)


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
simu_df <- list_rpar[[100]] # first simulation
nsites <- ncol(simu_df) # number of sites
par(mfrow = c(1, 1))
plot(simu_df[, 20], main = "rpareto simulation")
# get grid coordinates
sites_coords <- generate_grid_coords(sqrt(nsites))
gif_folder <- "/user/cserreco/home/Documents/These/phd_extremes/thesis/resources/images/simulation/rpar"
create_simu_gif(simu_df, site_coords, c(param, adv), type = "rpar",
                foldername = gif_folder)

length(list_rpar)
# concat all simulations together
simu_all <- do.call(rbind, list_rpar)
gif_folder <- "/user/cserreco/home/Documents/These/phd_extremes/thesis/resources/images/simulation/rpar"
create_simu_gif(simu_all, site_coords, c(param, adv), type = "rpar",
                foldername = gif_folder, forcedtemp = 200)


params <- c(param, adv)

# get the lag vectors

# Spatio-temporal neighborhood parameters
min_spatial_dist <- 3 # in km
delta <- 10 # in hours
episode_size <- 2 * delta # size of the episode
threshold <- 1

# Get the selected episodes
selected_points <- select_extreme_episodes(sites_coords, simu_all,
                                        min_spatial_dist,
                                        threshold = threshold,
                                        episode_size = episode_size,
                                        n_max_episodes = 10000,
                                        time_ext = 0, latlon = FALSE)

n_episodes <- length(selected_points$s0)
print(n_episodes)
length(unique(selected_points$s0)) # can be same s0
length(unique(selected_points$t0)) # never same t0?
print(min(selected_points$u_s0)) # min threshold
t0_list <- selected_points$t0
s0_list <- selected_points$s0
list_episodes_points <- get_extreme_episodes(selected_points, simu_all,
                                     episode_size = episode_size, unif = FALSE,
                                     beta = 0)

list_episodes <- list_episodes_points$episodes
episode <- list_episodes[[1]]
# check the episode
head(episode)
library(ggplot2)
library(reshape2)  # for melting wide data to long format

# Convert matrix to data frame
index <- 3
sort(t0_list)
which(t0_list == t0_list[index])
episode_test <- list_episodes[[index]]
df_episode <- as.data.frame(episode_test)
df_episode$Time <- 0:(nrow(df_episode) - 1)  # Add a time column
s0_list[index]
u_episode <- selected_points$u_s0[index]
# t0_episode <- t0_list[index]
# Convert from wide to long format
df_long <- melt(df_episode, id.vars = "Time")
head(df_long)
colnames(df_long) <- c("Time", "Pixel", "Value")
ggplot(df_long, aes(x = Time, y = Value, group = Pixel)) +
  geom_line(color = btfgreen) +
  geom_hline(yintercept = u_episode, color = "red", linetype = "dashed") +
  # labs(title = "Extreme Episode", x = "Time", y = "Value") +
  # annotate("text", x = t0_episode, y = 0, label = expression(t[0]),
  #          color = "red", vjust = 0, size = 7) +
  theme_minimal()

ggsave("image.png", width = 8, height = 6)

filename <- paste(im_folder, "optim/simulation/rpar_",
                  param_str, "/episodes_",
                  ngrid^2, "s_", length(temp), "t_",
                  "t0_", t0_str, "_",
                  s0_list[index], ".png", sep = "")
# filename <- "test.png"
ggsave(filename, width = 20, height = 15, units = "cm")

library(ggplot2)
library(reshape2)

unique_s0_list <- unique(s0_list)

for (target_s0 in unique_s0_list) {
  matching_indices <- which(s0_list == target_s0)

  df_list <- list()

  for (i in matching_indices) {
    episode_test <- list_episodes[[i]]
    df_episode <- as.data.frame(episode_test)
    df_episode$Time <- 0:(nrow(df_episode) - 1)

    df_long <- melt(df_episode, id.vars = "Time", variable.name = "Pixel",
                    value.name = "Value")
    colnames(df_long) <- c("Time", "Pixel", "Value")
    df_filtered <- subset(df_long, Pixel == target_s0)
    df_filtered$Episode <- as.factor(i)
    df_list[[i]] <- df_filtered
  }

  df_all <- do.call(rbind, df_list)
  u_episode <- selected_points$u_s0[matching_indices[1]]
  p <- ggplot(df_all, aes(x = Time, y = Value, group = Episode)) +
    geom_line(color = btfgreen, alpha = 0.7) +
    labs(x = "Time", y = "Rainfall (mm)") +
    theme_minimal() +
    geom_hline(yintercept = u_episode, color = "red", linetype = "dashed")

  # Find the min and max values of the 'Value' for setting the position of the annotation
  y_min <- min(df_all$Value, na.rm = TRUE)
  y_max <- max(df_all$Value, na.rm = TRUE)

  # Set the annotation's y position dynamically based on the min value
  annotation_y <- y_min - 0.1 * (y_max - y_min)

  # Add vertical line at t0 and annotation for t0
  # p <- p + geom_vline(xintercept = delta, color = "#ff0000a6", linetype = "dashed") +
  #   annotate("text", x = delta + 0.5, y = annotation_y, label = expression(t[0]),
  #            color = "red", vjust = 0, size = 5)

  filename <- paste0(im_folder,
                     "optim/simulation/rpar_", param_str, "/episodes_",
                     ngrid^2, "s_", length(temp), "t_",
                     "t0_", t0_str, "_",
                     "s0_", target_s0, "_",
                     "episodes_", n_episodes, "_min",
                     min_spatial_dist, "km", "_delta_", delta, ".png")
  ggsave(filename, plot = p, width = 20, height = 15, units = "cm")

}


list_episodes_unif_points <- get_extreme_episodes(selected_points, simu_all,
                                episode_size = episode_size, unif = TRUE)

list_episodes_unif <- list_episodes_unif_points$episodes
s0_list <- selected_points$s0
t0_list <- selected_points$t0
u_list <- selected_points$u_s0


library(parallel)

rownames(sites_coords) <- paste0("S", 1:nrow(sites_coords))
tau_vect <- 0:10
u <- 1
tmax <- max(tau_vect)
sites_coords <- as.data.frame(sites_coords)
# Compute the lags and excesses for each conditional point
list_results <- mclapply(1:length(s0_list), function(i) {
  s0 <- s0_list[i]
  s0_coords <- sites_coords[s0, ]
  # t0 <- t0_list[i]
  episode <- list_episodes[[i]]
  ind_t0_ep <- 0 # index of t0 in the episode
  lags <- get_conditional_lag_vectors(sites_coords, s0_coords, ind_t0_ep,
                                tau_vect, latlon = FALSE)
  # lags$hnorm <- lags$hnorm / 1000 # convert to km
  excesses <- empirical_excesses(episode, u, lags, type = "rpareto",
                                t0 = ind_t0_ep, threshold = TRUE)
  list(lags = lags, excesses = excesses)
}, mc.cores = detectCores() - 1)

list_lags <- lapply(list_results, `[[`, "lags")
list_excesses <- lapply(list_results, `[[`, "excesses")

s0 <- s0_list[1]
s0_coords <- sites_coords[s0, ]
excesses <- list_excesses[[1]]
# sort excesses by hnorm from smallest to largest
excesses <- excesses[order(excesses$hnorm), ]
excesses$kij
excesses[1:20, ]
df_lags <- list_lags[[1]]
head(df_lags)


wind_df <- as.data.frame(t(adv))
colnames(wind_df) <- c("vx", "vy")

init_param <- c(beta1, beta2, alpha1, alpha2, 1, 1)
init_param <- params
result <- optim(par = init_param, fn = neg_ll_composite,
        list_lags = list_lags, list_episodes = list_episodes,
        list_excesses = list_excesses, hmax = 7,
        latlon = FALSE,
        directional = FALSE,
        method = "L-BFGS-B",
        lower = c(1e-08, 1e-08, 1e-08, 1e-08, 1e-08, 1e-08),
        upper = c(Inf, Inf, 1.999, 1.999, Inf, Inf),
        control = list(maxit = 10000,
                      trace = 1,
                      parscale = c(1, 1, 1, 1, 1, 1)),
        hessian = F)
result







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

true_df <- data.frame(value = true_param_eta,
                param = param_names,
                type = "True")
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

# Save the plot
ggsave("./images/optim/boxplot_25s_300t_rpar.png", width = 8, height = 6)





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
