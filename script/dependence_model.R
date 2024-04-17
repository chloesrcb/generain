library(generain)

################################################################################
# Validation of the dependence model -------------------------------------------
################################################################################

# Brown-Resnick simulation data ------------------------------------------------

# true parameters
# beta1, beta2, alpha1, alpha2
true_param <- c(0.4, 0.2, 1.5, 1)

# Simulation with 25 sites and 300 times
list_BR_25_300 <- list()
for (i in 1:100) {
  file_path <- paste0("../data/simulations_BR/sim_25s_300t/rainBR_", i, ".csv")
  df <- read.csv(file_path)
  list_BR_25_300[[i]] <- df
}

# Simulation with 49 sites and 100 times
list_BR_49_100 <- list()
for (i in 1:100) {
  file_path <- paste0("../data/simulations_BR/sim_49s_100t/rainBR_", i, ".csv")
  df <- read.csv(file_path)
  list_BR_49_100[[i]] <- df
}

# Simulation with 9 sites and 50 times and advection
list_BR_adv_9_50 <- list()
for (i in 1:10) {
  file_path <- paste0("../data/simulations_BR/sim_adv_9s_50t_10/rainBR_", i, ".csv")
  df <- read.csv(file_path)
  list_BR_adv_9_50[[i]] <- df
}

# Evaluation of the dependence model

# Simulation with 25 sites and 300 times ---------------------------------------
simu_df <- list_BR_25_300[[1]] # first simulation
nsites <- ncol(simu_df) # number of sites
df_dist <- distances_regular_grid(nsites) # distance matrix

# Evaluate the estimates
spa_estim_25 <- evaluate_vario_estimates(list_BR_25_300, 0.9,
                                      c(param[1], param[3]), spatial = TRUE,
                                      df_dist = df_dist, hmax = sqrt(17))

temp_estim_25 <- evaluate_vario_estimates(list_BR_25_300, 0.9,
                      c(param[2], param[4]), spatial = FALSE, tmax = 10)

# Simulation with 49 sites and 100 times ---------------------------------------

simu_df <- list_BR_49_100[[1]] # first simulation
nsites <- ncol(simu_df) # number of sites
df_dist <- distances_regular_grid(nsites) # distance matrix

# Evaluate the estimates
spa_estim_49 <- evaluate_vario_estimates(list_BR_49_100, 0.8,
                                      c(param[1], param[3]), spatial = TRUE,
                                      df_dist = df_dist, hmax = sqrt(17))

temp_estim_49 <- evaluate_vario_estimates(list_BR_49_100, 0.8,
                      c(param[2], param[4]), spatial = FALSE, tmax = 10)

# Simulation with 9 sites and 50 times and advection --------------------------

simu_df <- list_BR_adv_9_50[[1]] # first simulation
nsites <- ncol(simu_df) # number of sites
df_dist <- distances_regular_grid(nsites) # distance matrix

# Evaluate the estimates
spa_estim_adv <- evaluate_vario_estimates(list_BR_adv_9_50, 0.8,
                                      c(param[1], param[3]), spatial = TRUE,
                                      df_dist = df_dist, hmax = sqrt(17))

temp_estim_adv <- evaluate_vario_estimates(list_BR_adv_9_50, 0.8,
                      c(param[2], param[4]), spatial = FALSE, tmax = 10)