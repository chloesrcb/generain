
# remove objects from the environment
rm(list = ls())

# clear console
cat("\014")

# LOAD LIBRARIES ###############################################################
library(parallel)
library(abind)
muse <- TRUE

if (muse) {
  # Get the muse folder
  folder_muse <- "/home/serrec/work_rainstsimu/rpareto/"
  setwd(folder_muse)
  # Load libraries and set theme
  source("load_libraries.R")
  source("config.R")
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
result_folder <- "./data/simulations_rpar/"
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
nres <- M * m

# Apply formatting
param_str <- format_value(true_param)
adv_str <- format_value(adv)
s0_str <- format_value(s0)
t0_str <- format_value(t0)

# Save the data
foldername <- paste0("./data/simulations_rpar/rpar_", param_str, "_", adv_str,
                   "/sim_", ngrid^2, "s_", length(temp), "t_s0_",
                    s0_str, "_t0_", t0_str, "/")

if (!dir.exists(foldername)) {
  dir.create(foldername, recursive = TRUE)
}


# Parallel configuration
num_cores <- detectCores() - 1
cl <- makeCluster(num_cores)

# Export variables and functions to the cluster
clusterExport(cl, varlist = c(
  "sim_rpareto", "beta1", "beta2", "alpha1", "alpha2",
  "spa", "temp", "adv", "t0", "m", "random_s0", "s0",
  "compute_gamma_point", "compute_W_s_t",
  "foldername", "generate_grid_coords", "save_simulations",
  "ngrid"
))

# Load libraries in the cluster
clusterEvalQ(cl, {
  library(RandomFields)
  library(dplyr)
  library(tidyr)
  RandomFields::RFoptions(
    spConform = FALSE,
    allow_duplicated_locations = TRUE,
    storing = FALSE,
    printlevel = 0
  )
})

# Simulate the process in parallel
result_list <- parLapply(cl, 1:M, function(i) {
  simu <- sim_rpareto(
    beta1 = beta1,
    beta2 = beta2,
    alpha1 = alpha1,
    alpha2 = alpha2,
    x = spa,
    y = spa,
    t = temp,
    adv = adv,
    t0 = t0,
    nres = m,
    random_s0 = random_s0,
    s0 = s0
  )

  # Save the simulation
  Z_rpar <- simu$Z
  s0_used_list <- simu$s0_used

  save_simulations(Z_rpar, ngrid,
                   folder = foldername,
                   file = paste0("rpar_", ngrid^2, "s_", length(temp),
                   "t_simu", i))

  return(list(Z = Z_rpar, s0_used = s0_used_list))
})

stopCluster(cl) # Stop the cluster

# Combine all Z_rpar matrices into one big 4D array
Z_list <- lapply(result_list, function(res) res$Z)

# Combine them along the 4th dimension (simulation axis)
Z_rpar <- abind::abind(Z_list, along = 4)

# Combine all s0_used lists into one big list
s0_list <- do.call(c, lapply(result_list, function(res) res$s0_used))

sites_coords <- generate_grid_coords(ngrid)

# Save the data
files <- list.files(foldername, full.names = TRUE)
length(files)
list_rpar <- list()
for (i in 1:M) {
  for (j in 1:m) {
    file_name <- paste0(foldername, "rpar_", ngrid^2, "s_",
                        length(temp), "t_simu", i, "_", j, ".csv")
    list_rpar[[((i-1) * m + j)]] <- read.csv(file_name)
  }
}

### Optimization ###############################################################
library(parallel)
num_cores <- detectCores() - 1  # Reserve 1 core for the OS

# Parallel execution
tau_vect <- 0:10
u <- 1
tmax <- max(tau_vect)
sites_coords <- as.data.frame(sites_coords)

# Compute lags and excesses for each episode inside the simulation
list_results <- mclapply(1:nres, function(i) {
  s0_x <- s0_list[[i]]$x
  s0_y <- s0_list[[i]]$y
  s0_coords <- sites_coords[sites_coords$Longitude == s0_x &
                            sites_coords$Latitude == s0_y, ]
  episode <- list_rpar[[i]]
  lags <- get_conditional_lag_vectors(sites_coords, s0_coords, t0,
                                      tau_vect, latlon = FALSE)
  excesses <- empirical_excesses(episode, u, lags, type = "rpareto",
                                 t0 = t0, threshold = TRUE)
  list(lags = lags, excesses = excesses)
}, mc.cores = num_cores)

# Extract lags and excesses from the list of results
list_lags <- lapply(list_results, `[[`, "lags")
list_excesses <- lapply(list_results, `[[`, "excesses")

result_list <- mclapply(1:M, process_simulation, M = M, m = m,
                        list_simu = list_rpar, u = u,
                        list_lags = list_lags,
                        list_excesses = list_excesses,
                        init_params = true_param,
                        hmax = 7, wind_df = NA,
                        mc.cores = num_cores)

# Combine results int a data frame
df_result_all <- do.call(rbind, result_list)
colnames(df_result_all) <- c("beta1", "beta2", "alpha1",
                             "alpha2", "adv1", "adv2")
# Convert in latex
colnames(df_result_all) <- c(
  TeX("$\\beta_1$"), TeX("$\\beta_2$"),
  TeX("$\\alpha_1$"), TeX("$\\alpha_2$"),
  TeX("$v_x$"), TeX("$v_y$")
)
df_bplot <- as.data.frame(df_result_all)

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

labels_latex <- c(
  beta1 = TeX("$\\beta_1$"),
  beta2 = TeX("$\\beta_2$"),
  alpha1 = TeX("$\\alpha_1$"),
  alpha2 = TeX("$\\alpha_2$"),
  vx = TeX("$v_x$"),
  vy = TeX("$v_y$")
)

bplot <- ggplot(df_bplot, aes(x = ind, y = values)) +
  geom_boxplot() +
  labs(title = "",
    x = "Parameters", y = "Estimated values") +
  theme_minimal() +
  geom_point(aes(y = true_param[as.numeric(ind)]), color = "red", pch=4) +
  scale_x_discrete(labels = labels_latex) 


# folder to save the plot
foldername <- "./images/optim/rpar/"
if (!dir.exists(foldername)) {
  dir.create(foldername, recursive = TRUE)
}


# Save the plot
name_file <- paste0("bp_optim_", M, "simu_", m, "rep_", ngrid^2,
                "s_", length(temp), "t_", param_str, "_",
                adv_str, ".png")

ggsave(plot = bplot, filename = paste0(foldername, name_file),
       width = 5, height = 5)
