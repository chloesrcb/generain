
# remove objects from the environment
rm(list = ls())

# clear console
cat("\014")

# LOAD LIBRARIES ###############################################################
muse <- FALSE

if (muse) {
  # Get the muse folder
  folder_muse <- "/home/serrec/work_rainstsimu/rpareto/"
  setwd(folder_muse)
  # Load libraries and set theme
  source("load_libraries.R")
  source("config.R")
  data_folder <- "./data/"
} else {
  # Load libraries and set theme
  source("./script/load_libraries.R")
  source("./script/optimisation/config.R")
}

# Get all files in the folder "R"
functions_folder <- "./R"
files <- list.files(functions_folder, full.names = TRUE)
# load all functions in files
invisible(lapply(files, function(f) source(f, echo = FALSE)))
library(latex2exp)

# SIMULATION ###################################################################
result_folder <- paste0(data_folder, "simulations/simulations_rpar/")
if (!dir.exists(result_folder)) {
  dir.create(result_folder, recursive = TRUE)
}

# Configuration
true_param <- c(params, adv)
beta1 <- params[1]
beta2 <- params[2]
alpha1 <- params[3]
alpha2 <- params[4]
spa <- 1:ngrid
nsites <- ngrid^2 # if the grid is squared

# Number of realizations
M <- 1 # number of simulations
m <- 500 # number of extreme episodes
nres <- M * m

# Function to format values correctly, considering decimal places
format_value <- function(x) {
  formatted_values <- sapply(x, function(val) {
    # Vérifier si val est négatif
    is_negative <- val < 0
    val <- abs(val)  # Travailler avec la valeur absolue pour le formatage
    
    # Vérifier si val est un entier
    if (val == as.integer(val)) {
      formatted_value <- sprintf("%d", val)
    } else {
      # Compter le nombre de décimales
      num_decimals <- nchar(sub("^[^.]*\\.", "", as.character(val)))
      
      if (val >= 1) {
        # Si le nombre est supérieur ou égal à 1
        if (num_decimals == 1) {
          formatted_value <- sprintf("%02d", round(val * 10))
        } else {
          formatted_value <- sprintf("%03d", round(val * 100))
        }
      } else {
        # Si le nombre est inférieur à 1
        if (num_decimals == 1) {
          formatted_value <- sprintf("%02d", round(val * 10))
        } else {
          formatted_value <- sprintf("%03d", round(val * 100))
        }
      }
    }
    
    # Ajouter "neg" devant si le nombre était négatif
    if (is_negative) {
      formatted_value <- paste0("neg", formatted_value)
    }
    
    return(formatted_value)
  })
  
  # Concaténer les valeurs avec un underscore "_"
  return(paste(formatted_values, collapse = "_"))
}

# Apply formatting
param_str <- format_value(true_param)
adv_str <- format_value(adv)
s0_str <- format_value(s0)
t0_str <- format_value(t0)

# Save the data
foldername <- paste0(data_folder, "simulations/simulations_rpar/rpar_", param_str, 
                   "/sim_", ngrid^2, "s_", length(temp), "t_s0_",
                    s0_str, "_t0_", t0_str, "/")

if (!dir.exists(foldername)) {
  dir.create(foldername, recursive = TRUE)
}

# save_simulations(simu, ngrid, nres, folder = foldername,
#         file = paste0("rpar_", ngrid^2, "s_", length(temp), "t"))
# get files in the folder
# foldername = "/home/cserreco/Documents/These/phd_extremes/data/simulations/simulations_rpar/rpar_001_02_15_1_02_01/sim_25s_30t_s0_1_1_t0_14/"
foldername = "/user/cserreco/home/Documents/These/phd_extremes/data/simulations/simulations_rpar/rpar_001_02_15_1_02_01/sim_49s_30t_s0_1_1_t0_0/"
files <- list.files(foldername, full.names = TRUE)
length(files)
list_simuM <- list()
m <- 500
nres <- m*M
for (i in 1:nres) {
  file_name <- files[i]
  list_simuM[[i]] <- read.csv(file_name)
}

adv <- c(0.2, 0.1)
true_param <- c(0.01, 0.2, 1.5, 1, adv)
# create a wind dataframe with m rows and adv at each row
wind_df <- data.frame(matrix(rep(adv, m), nrow = m, byrow = TRUE))
colnames(wind_df) <- c("vx", "vy")
### Optimization ###############################################################
library(parallel)
num_cores <- detectCores() - 1  # Reserve 1 core for the OS
ngrid <- 7
# Parallel execution
t0 = 0
tau_vect <- 0:10
sites_coords <- generate_grid_coords(ngrid)
s0 <- c(1, 1)
df_lags <- get_conditional_lag_vectors(sites_coords, s0, t0, tau_vect)
# init_param <- true_param
# init_param[5:6] <- c(1, 1)
u <- 1 # threshold corresponding to the r-pareto simulation

i=1
# Get the m corresponding simulations from list_simu inside a list
list_episodes <- list_simuM[((i - 1) * m + 1):(i * m)]
# Compute excesses
list_excesses <- lapply(list_episodes, function(episode) {
empirical_excesses_rpar(episode, u, df_lags, threshold = TRUE, t0 = t0)
})

episode <- list_episodes[[1]]
excesses <- list_excesses[[1]]
excesses$kij
hmax <- sqrt(17)
# create a list of lags with df_lags for all 1:m
list_lags <- lapply(1:m, function(i) {
  df_lags
})
list_lags[[1]]

# Optimize
result <- optim(
    par = c(0.01, 0.2, 1.5, 1, 0.5, 1),
    fn = neg_ll_composite,
    list_episodes = list_episodes,
    list_lags = list_lags,
    list_excesses = list_excesses,
    hmax = hmax,
    wind_df = adv,
    threshold = TRUE,
    latlon = FALSE,
    directional = FALSE,
    method = "L-BFGS-B",
    lower = c(1e-8, 1e-8, 1e-8, 1e-8, 1e-8,  1e-8),
    upper = c(Inf, Inf, 1.999, 1.999, Inf, Inf),
    control = list(maxit = 10000)
)
result


eta1_vals <- seq(0.05, 0.3, length.out = 20)
eta2_vals <- seq(0.05, 0.2, length.out = 20)
grid_ll <- matrix(NA, nrow = length(eta1_vals), ncol = length(eta2_vals))
for (i in seq_along(eta1_vals)) {
  for (j in seq_along(eta2_vals)) {
    p_try <- c(true_param[1:4], eta1_vals[i], eta2_vals[j])
    grid_ll[i, j] <- neg_ll_composite(p_try, list_episodes, list_excesses, list_lags,
                                      threshold = TRUE, latlon = FALSE,
                                      directional = FALSE)
  }
}

contour(eta1_vals, eta2_vals, grid_ll, xlab = "eta1", ylab = "eta2", main = "Profile likelihood")


# Optim with eta1 and eta2
result <- optim(
    par = c(0.01, 0.2, 1.5, 1, 1, 1),
    fn = neg_ll_composite,
    list_episodes = list_episodes,
    list_lags = list_lags,
    list_excesses = list_excesses,
    hmax = hmax,
    wind_df = adv,
    threshold = TRUE,
    latlon = FALSE,
    directional = TRUE,
    method = "L-BFGS-B",
    lower = c(1e-8, 1e-8, 1e-8, 1e-8, 1e-8, 1e-8),
    upper = c(Inf, Inf, 1.999, 1.999, Inf, Inf),
    control = list(maxit = 10000)
)
result

params <- c(0.01, 0.2, 1.5, 1, 0.2, 0.1)
df_chi = theoretical_chi(params, df_lags, latlon =FALSE, directional = FALSE)


result <- optim(
    par = c(0.01, 0.2, 1.5, 1, 0.3, 0.1),
    fn = neg_ll_composite_simu,
    list_simu = list_episodes,
    df_lags = df_lags,
    list_excesses = list_excesses,
    quantile = 1,
    hmax = hmax,
    adv = NA,
    threshold = TRUE,
    method = "L-BFGS-B",
    lower = c(1e-8, 1e-8, 1e-8, 1e-8, -Inf, -Inf),
    upper = c(Inf, Inf, 1.999, 1.999, Inf, Inf),
    control = list(maxit = 10000)
)







result_list <- mclapply(1:M, process_simulation, M = M, m = m,
                        list_simuM = list_simuM, u = u, df_lags = df_lags,
                        t0 = t0, true_param = init_param,
                        wind_df = wind_df, hmax = sqrt(17),
                        mc.cores = num_cores)

# Combine results into a data frame
df_result_all <- do.call(rbind, result_list)
colnames(df_result_all) <- c("beta1", "beta2", "alpha1",
                             "alpha2", "adv1", "adv2")

df_bplot <- as.data.frame(df_result_all)

# > df_result_all adv +, -10:10
#              [,1]      [,2]     [,3]      [,4]       [,5]       [,6]
#  [1,] 0.010661974 0.4029309 1.528410 1.0052373  0.2289262  0.2780111
#  [2,] 0.008637025 0.4160721 1.687900 0.9440856 -0.3333151 -0.4368077
#  [3,] 0.010263775 0.3915771 1.403014 1.0234373 -0.2016857 -0.2097121
#  [4,] 0.009787595 0.3889252 1.570193 1.0080372  0.3264038  0.3396228
#  [5,] 0.012557590 0.4144650 1.433984 0.9813914  0.2988893  0.3090406
#  [6,] 0.009846980 0.3572470 1.494960 1.0426908 -0.2784679 -0.2768801
#  [7,] 0.008069802 0.3751421 1.556749 1.0013524 -0.6856584 -0.6677612
#  [8,] 0.010075535 0.4053543 1.463154 0.9884304 -0.4487910 -0.3246396
#  [9,] 0.009995956 0.4071360 1.496522 0.9561426  0.6162234  0.6301404
# [10,] 0.010795546 0.3990308 1.499959 0.9992351  0.1001270  0.1000996

# > result_list adv -, -10:10
# [1]  0.01152073  0.40703984  1.43689895  0.99820642 -0.25572222 -0.35909073
# [1] 0.007704167 0.421337853 1.808492050 0.928379014 0.353110463 0.535483230
# [1] 0.009565259 0.391407200 1.471368674 1.021037043 0.206454886 0.238160531
# [1]  0.01052999  0.38731634  1.49663297  1.00795458 -0.36977565 -0.36776822
# [1]  0.01391657  0.41565761  1.32979420  0.98096582 -0.30642827 -0.32322678
# [1] 0.0101511 0.3565449 1.4914919 1.0433028 0.2737667 0.2993415
# [1] 0.006618769 0.373598004 1.735770212 0.992054520 0.747725312 0.716338517
# [1] 0.01000697 0.40325170 1.46085957 0.99187402 0.39389899 0.30198489
# [1]  0.01016631  0.40431457  1.49671136  0.94789362 -0.74466754 -0.79885082
# [1]  0.01217717  0.40107062  1.39982817  0.99351691 -0.28077293 -0.23116656
# chi_df <- df_lags[c("s1", "s2", "tau", "s1x", "s1y", "s2x", "s2y", "hnorm")]

# chi_df$s1xv <- chi_df$s1x
# chi_df$s1yv <- chi_df$s1y
# chi_df$s2xv <- chi_df$s2x + adv[1] * chi_df$tau
# chi_df$s2yv <- chi_df$s2y + adv[2] * chi_df$tau
# chi_df$dx <- chi_df$s2xv - chi_df$s1x  # Distance corrigée en x
# chi_df$dy <- chi_df$s2yv - chi_df$s1y  # Distance corrigée en y
# chi_df$dx_orig <- chi_df$s2x - chi_df$s1x  # Distance originale en x
# chi_df$dy_orig <- chi_df$s2y - chi_df$s1y  # Distance originale en y

# save data in csv
foldername <- "./data/optim/rpar/"
if (!dir.exists(foldername)) {
  dir.create(foldername, recursive = TRUE)
}

name_file <- paste0("optim_rpar_", M, "simu_", m, "rep_", ngrid^2,
                "s_", length(temp), "t_", param_str, "_",
                adv_str, ".csv")
write.csv(df_bplot, paste0(foldername, name_file), row.names = FALSE)


# Plot results
df_bplot <- stack(df_bplot)

ggplot(df_bplot, aes(x = ind, y = values)) +
  geom_boxplot() +
  labs(title = "",
    x = "Parameters", y = "Estimated values") +
  theme_minimal() +
  geom_point(aes(y = true_param[as.numeric(ind)]), color = "red", pch=4) 


# save plot
foldername <- "./images/optim/rpar/"
if (!dir.exists(foldername)) {
  dir.create(foldername, recursive = TRUE)
}


# Save the plot
name_file <- paste0("bp_optim_", M, "simu_", m, "rep_", ngrid^2,
                "s_", length(temp), "t_", param_str, "_",
                adv_str, ".png")

ggsave(paste0(foldername, name_file), width = 5, height = 5)
