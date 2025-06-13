
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
u <- 1
tau_vect <- 0:10

# Apply formatting
param_str <- format_value(true_param)
adv_str <- format_value(adv)
s0_str <- format_value(s0)
t0_str <- format_value(t0)

# Save the data
foldername <- paste0("./data/simulations_rpar/rpar_test_", param_str, "_",
                    adv_str, "/sim_", ngrid^2, "s_", length(temp), "t_s0_",
                    s0_str, "_t0_", t0_str, "/")

if (random_s0) {
  foldername <- paste0(foldername, "random_s0/")
} else {
  foldername <- paste0(foldername, "fixed_s0/")
}

if (is_anisotropic) {
  foldername <- paste0(foldername, "anisotropic/")
} else {
  foldername <- paste0(foldername, "isotropic/")
}

if (use_wind_data) {
  wind_df <- data.frame(vx = adv[1], vy = adv[2])
  foldername <- paste0(foldername, "wind/")
} else {
  wind_df <- NA
  foldername <- paste0(foldername, "no_wind/")
}

if (!dir.exists(foldername)) {
  print(paste0("Folder created: ", foldername))
  dir.create(foldername, recursive = TRUE)
}

sites_coords <- as.data.frame(generate_grid_coords(ngrid))

# Parallel configuration
num_cores <- detectCores() - 1
cl <- makeCluster(num_cores)

# Export variables and functions to the cluster
clusterExport(cl, varlist = c(
  "sim_rpareto", "beta1", "beta2", "alpha1", "alpha2",
  "spa", "temp", "adv", "t0", "m", "random_s0", "s0",
  "s0_radius", "is_anisotropic",
  "foldername", "generate_grid_coords", "save_simulations",
  "ngrid", "compute_st_gaussian_process",
  "compute_st_variogram", "convert_simulations_to_list",
  "empirical_excesses", "get_conditional_lag_vectors",
  "process_simulation_mem", "u", "t0", "tau_vect",
  "sites_coords", "theoretical_chi", "neg_ll_composite_mem",
  "neg_ll", "true_param", "empirical_excesses_rpar",
  "fixed_eta1", "fixed_eta2", "wind_df"
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
    s0 = s0,
    s0_radius = s0_radius,
    anisotropic = is_anisotropic
  )

  list_rpar <- convert_simulations_to_list(simu$Z, ngrid)
  s0_list <- simu$s0_used

  result <- process_simulation_mem(
    list_episodes = list_rpar,
    u = u,
    s0_list = s0_list,
    sites_coords = sites_coords,
    init_params = true_param,
    directional = is_anisotropic,
    hmax = 7,
    wind_df = wind_df,
    fixed_eta1 = fixed_eta1,
    fixed_eta2 = fixed_eta2
  )

  return(result)
})

stopCluster(cl) # Stop the cluster


# Combine results int a data frame
df_result_all <- do.call(rbind, result_list)
colnames(df_result_all) <- c("beta1", "beta2", "alpha1",
                             "alpha2", "adv1", "adv2")

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


df_bplot <- stack(df_bplot)

# Plot results
labels_latex <- c(
  beta1 = TeX("$\\beta_1$"),
  beta2 = TeX("$\\beta_2$"),
  alpha1 = TeX("$\\alpha_1$"),
  alpha2 = TeX("$\\alpha_2$"),
  eta1 = TeX("$\\eta_1$"),
  eta2 = TeX("$\\eta_2$"),
  adv1 = TeX("$v_x$"),
  adv2 = TeX("$v_y$")
)

initial_param <- c(
  beta1 = beta1,
  beta2 = beta2,
  alpha1 = alpha1,
  alpha2 = alpha2,
  eta1 = NA,
  eta2 = NA,
  adv1 = adv[1],
  adv2 = adv[2]
)

params_group1 <- c("beta1", "beta2", "alpha1", "alpha2")
params_group2 <- c("adv1", "adv2")

bplot1 <- ggplot(subset(df_bplot, ind %in% params_group1), aes(x = ind, y = values)) +
  geom_boxplot() +
  geom_point(aes(y = initial_param[match(ind, names(initial_param))]), color = "red", pch = 4) +
  scale_x_discrete(labels = labels_latex[params_group1]) +
  labs(x = "Parameters", y = "Estimated values") +
  theme_minimal()

bplot2 <- ggplot(subset(df_bplot, ind %in% params_group2), aes(x = ind, y = values)) +
  geom_boxplot() +
  geom_point(aes(y = initial_param[match(ind, names(initial_param))]), color = "red", pch = 4) +
  scale_x_discrete(labels = labels_latex[params_group2]) +
  labs(x = "Parameters", y = "Estimated values") +
  theme_minimal()


# folder to save the plot
foldername <- "./images/optim/rpar/"
if (random_s0) {
  foldername <- paste0(foldername, "random_s0/")
} else {
  foldername <- paste0(foldername, "fixed_s0/")
}

if(is_anisotropic) {
  foldername <- paste0(foldername, "anisotropic/")
} else {
  foldername <- paste0(foldername, "isotropic/")
}

if (!dir.exists(foldername)) {
  dir.create(foldername, recursive = TRUE)
  print(paste0("Folder created: ", foldername))
}

# Save the plot
name_file <- paste0("bp_optim_", M, "simu_", m, "rep_", ngrid^2,
                "s_", length(temp), "t_", param_str, "_beta_alpha.png")

ggsave(plot = bplot1, filename = paste0(foldername, name_file),
       width = 5, height = 5)

# Save the plot
name_file <- paste0("bp_optim_", M, "simu_", m, "rep_", ngrid^2,
              "s_", length(temp), "t_", param_str, "_adv.png")

ggsave(plot = bplot2, filename = paste0(foldername, name_file),
       width = 5, height = 5)
