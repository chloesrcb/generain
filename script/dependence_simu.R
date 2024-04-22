library(generain)

################################################################################
# Validation of the dependence model -------------------------------------------
################################################################################

# Brown-Resnick simulation data ------------------------------------------------

# true parameters
# beta1, beta2, alpha1, alpha2
param <- c(0.4, 0.2, 1.5, 1)

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

# Simulation with 400 sites and 50 times
list_BR_400_50 <- list()
for (i in 1:100) {
  file_path <- paste0("../data/simulations_BR/sim_400s_50t/rainBR_", i, ".csv")
  df <- read.csv(file_path)
  list_BR_400_50[[i]] <- df
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
                                  true_param = c(param[1], param[3]),
                                  spatial = TRUE, df_dist = df_dist,
                                  hmax = sqrt(17))

temp_estim_25 <- evaluate_vario_estimates(list_BR_25_300, 0.9,
                      c(param[2], param[4]), spatial = FALSE, tmax = 10)

df_result <- cbind(spa_estim_25, temp_estim_25)
colnames(df_result) <- c("beta1", "alpha1", "beta2", "alpha2")

df_valid <- get_criterion(df_result, param)

# Simulation with 400 sites and 50 times ---------------------------------------

simu_df <- list_BR_400_50[[1]] # first simulation
nsites <- ncol(simu_df) # number of sites
df_dist <- distances_regular_grid(nsites) # distance matrix

# Evaluate the estimates
spa_estim_400 <- evaluate_vario_estimates(list_BR_400_50, 0.8,
                                      c(param[1], param[3]), spatial = TRUE,
                                      df_dist = df_dist, hmax = sqrt(17))

temp_estim_400 <- evaluate_vario_estimates(list_BR_400_50, 0.8,
                      c(param[2], param[4]), spatial = FALSE, tmax = 10)


df_result_400 <- cbind(spa_estim_400, temp_estim_400)
colnames(df_result_400) <- c("beta1", "alpha1", "beta2", "alpha2")

df_valid_400 <- get_criterion(df_result_400, true_param)

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


df_result_49 <- cbind(spa_estim_49, temp_estim_49)
colnames(df_result_49) <- c("beta1", "alpha1", "beta2", "alpha2")

df_valid_49 <- get_criterion(df_result_49, param)

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