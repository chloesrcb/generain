library(generain)
library(reshape2)
library(ggplot2)
setwd("./script")
# load libraries
source("load_libraries.R")

# spatial and temporal structures
ngrid <- 5
spa <- 1:ngrid
nsites <- ngrid^2 # if the grid is squared
temp <- 1:300

# beta1, beta2, alpha1, alpha2
param <- c(0.4, 0.2, 1.5, 1) # true parameters for the variogram
beta1 <- param[1]
beta2 <- param[2]
alpha1 <- param[3]
alpha2 <- param[4]
adv <- c(0, 0)
true_param <- c(beta1, beta2, alpha1, alpha2)

# BR <- sim_BR(param[1], param[2], param[3], param[4], spa, spa, temp, 1,
#                adv)
# save_simulations(BR, ngrid, 1,
#                   folder = paste0("../data/simulations_BR/sim_", ngrid^2, "s_",
#                                 length(temp), "t/"),
#                   file = paste0("br_", ngrid^2, "s_",
#                                 length(temp), "t"), forcedind = 1)


library(foreach)
library(doParallel)

n.BR <- 1
num_iterations <- 5

# Create a parallel backend
cl <- makeCluster(detectCores())
registerDoParallel(cl)

# Run the simulation in parallel
results <- foreach(i = 1:num_iterations, .combine = "c") %dopar% {
  BR <- sim_BR(param[1], param[2], param[3], param[4], spa, spa, temp, n.BR)
  save_simulations(BR, ngrid, n.BR,
                  folder = paste0("../data/simulations_BR/oldsim_", ngrid^2, "s_",
                                length(temp), "t/"),
                  file = paste0("br_", ngrid^2, "s_",
                                length(temp), "t"), forcedind = i)

}

# Stop the parallel backend
stopCluster(cl)


# plot(BR[5,4,,1], main = "BR simulation")


# save simulations to CSV files
# save_simulations(BR, ngrid, n.BR,
#                  folder = paste0("../data/simulations_BR/sim_", ngrid^2, "s_",
#                                 length(temp), "t/"),
#                  file = paste0("br_", ngrid^2, "s_",
#                                 length(temp), "t"))


# Simulation
num_iterations <- 1
list_BR <- list()
for (i in 1:num_iterations) {
  file_path <- paste0("../data/simulations_BR/oldsim_25s_300t/br_",
                      ngrid^2, "s_", length(temp), "t_", i, ".csv")
  df <- read.csv(file_path)
  list_BR[[i]] <- df
}

# dependece buhl
# simu_df <- BR
simu_df <- list_BR[[1]] # first simulation
# simu_df <- simulation_data
par(mfrow=c(1,1))
nsites <- ncol(simu_df) # number of sites
plot(simu_df[, 1])
# get grid coordinates
sites_coords <- generate_grid_coords(sqrt(nsites))
dist_mat <- get_dist_mat(sites_coords,
                         latlon = FALSE) # distance matrix
df_dist <- reshape_distances(dist_mat) # reshape the distance matrix


params <- true_param
nsites <- ncol(simu_df) # number of sites
# get grid coordinates
sites_coords <- generate_grid_coords(sqrt(nsites))
# sites_coords$id <- 1:nsites
h_vect <- get_lag_vectors(sites_coords, params,
                          hmax = sqrt(17), tau_vect = 0)

# Evaluate the estimates with Buhl method
spa_estim <- evaluate_vario_estimates(list_BR, 0.8,
                                  spatial = TRUE, df_dist = df_dist,
                                  hmax = sqrt(17))

temp_estim <- evaluate_vario_estimates(list_BR, 0.85,
                                          spatial = FALSE, tmax = 10)

df_result <- data.frame(beta1 =  unlist(spa_estim$beta),
                        alpha1 = unlist(spa_estim$alpha),
                        beta2 = unlist(temp_estim$beta),
                        alpha2 = unlist(temp_estim$alpha))
colnames(df_result) <- c("beta1", "alpha1", "beta2", "alpha2")

true_param <- c(beta1, beta2, alpha1, alpha2)
df_valid <- get_criterion(df_result, true_param)

# save the results in latex table
write.table(round(df_valid, 3), file = "valid_br.tex",
            row.names = c("Beta1", "Alpha1", "Beta2", "Alpha2"),
            col.names = c("Mean", "RMSE", "MAE"),
            quote = FALSE)


################################################################################
# VERIFICATION
################################################################################


# qq-plots margins
par(mfrow = c(2, 2))


library(ismev)

BR_df <- simulation_data
BR_loc <- simu_df$S5
plot(BR_loc)
BR_loc_log <- log(BR_loc)
gumbel.fit <- gum.fit(BR_loc_log)

mu <- gumbel.fit$mle[1]
sigma <- gumbel.fit$mle[2]

library(evd)

theorical_qgum <- qgumbel(ppoints(BR_loc_log), mu, sigma)

qqplot(BR_loc_log, theorical_qgum, main = "Gumbel Q-Q plot",
       xlab = "Empirical quantiles",
       ylab = "Theoretical quantiles")


# fit frechet
library(extRemes)
frechet.fit <- fevd(BR_loc, type = "Gumbel", method = "MLE")

mu <- frechet.fit$results$par[1]
sigma <- frechet.fit$results$par[2]

theorical_qfrechet <- qgev(ppoints(BR_loc), mu, sigma, 0)

qqplot(BR_loc, theorical_qfrechet, main = "Frechet Q-Q plot",
  xlab = "Empirical quantiles",
  ylab = "Theoretical quantiles")

################################################################################
# SIMULATION WITH ADVECTION
################################################################################

adv <- c(0, 0)
param <- c(0.4, 0.2, 1.5, 1, adv)
ngrid <- 5
spa <- 1:ngrid
nsites <- ngrid^2 # if the grid is squared
temp <- 1:300
n.BR <- 1

if (all(adv == c(0, 0))) {
  foldername <- paste0("../data/simulations_BR/sim_", ngrid^2, "s_",
                                length(temp), "t/")
} else {
  foldername <- paste0("../data/simulations_BR/sim_", ngrid^2, "s_",
                                length(temp), "t_adv/")
}

if (!dir.exists(foldername)) {
  dir.create(foldername, recursive = TRUE)
}

# generate the simulations
# BR_adv <- sim_BR(true_param[1] * 2, true_param[2] * 2, true_param[3],
#                     true_param[4], spa, spa, temp, n.BR, adv)


library(foreach)
library(doParallel)
library(generain)

num_iterations <- 8

# Create a parallel backend
cl <- makeCluster(detectCores())
registerDoParallel(cl)

results <- foreach(i = 1:num_iterations, .combine = rbind) %dopar% {
  library(generain)
  BR <- generain::sim_BR(param[1], param[2], param[3],
                param[4], spa, spa, temp, 1, adv)

  save_simulations(BR, ngrid, 1, folder = foldername,
          file = paste0("br_", ngrid^2, "s_", length(temp), "t"),
          forcedind = i)

}

stopCluster(cl)

# model <- modelBuhlCklu
# library(parallel)
# advected_variogram <- function(x, y, z, grid, model, adv) {
#   lx <- length(x)
#   ly <- length(y)
#   lz <- length(z)
#   N <- nrow(grid)  # spatio-temporal dimension

#   # Initialize the result array
#   gamma <- array(0, dim = c(lx, ly, lz, N)) # variogram

#   # Precompute common values that will be reused in the loop
#   adv_x <- adv[1]
#   adv_y <- adv[2]

#   # Use parallelization to distribute the work across multiple cores
#   cl <- makeCluster(detectCores() - 1)  # Use all cores except one
#   print(model)
#   # Export necessary variables and functions to each worker
#   clusterExport(cl, varlist = c("x", "y", "z", "grid", "model",
#             "adv_x", "adv_y", "lx", "ly", "lz", "RFvariogram"))

#   # Apply parallelized loop using parLapply
#   gamma_list <- parLapply(cl, seq_len(N), function(n) {
#     print(model)
#     gamma_n <- array(0, dim = c(lx, ly, lz))  # Temporary array for current 'n'
#     for (i in seq_len(lx)) {
#       for (j in seq_len(ly)) {
#         for (k in seq_len(lz)) {
#           gamma_n[i, j, k] <- RandomFields::RFvariogram(
#             model,
#             x = x[i] - grid[n, 1] - adv_x * (z[k] - grid[n, 3]),
#             y = y[j] - grid[n, 2] - adv_y * (z[k] - grid[n, 3]),
#             z = z[k] - grid[n, 3]
#           )
#         }
#       }
#     }
#     return(gamma_n)
#   })

#   # Stop the cluster after computation
#   stopCluster(cl)

#   # Combine results back into the original array
#   for (n in seq_len(N)) {
#     gamma[, , , n] <- gamma_list[[n]]
#   }

#   return(gamma)
# }

# gamma = advected_variogram(x, y, z, grid, modelBuhlCklu, adv = c(0.1, 0.2))
