
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
  im_folder <- "./images"
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
result_folder <- paste0(data_folder, "simulations_rpar/coords/")
if (!dir.exists(result_folder)) {
  dir.create(result_folder, recursive = TRUE)
}

# Configuration
wind_data <- adv
adv_real <- c(eta1 * abs(adv[1])^eta2 * sign(adv[1]),
              eta1 * abs(adv[2])^eta2 * sign(adv[2]))
true_param <- c(params, adv_real)
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
# adv_str <- format_value(adv)
s0_str <- format_value(s0)
t0_str <- format_value(t0)

# Subfolder names
s0_type <- if (random_s0) "random_s0" else "fixed_s0"

eta_type <- if (!is.na(fixed_eta1) && !is.na(fixed_eta2)) {
  "fixed_eta"
} else if (!is.na(fixed_eta1)) {
  "fixed_eta1"
} else if (!is.na(fixed_eta2)) {
  "fixed_eta2"
} else {
  "free_eta"
}

wind_type <- if (use_wind_data) {
  "wind"
} else {
  "no_wind"
}

if (use_wind_data) {
  wind_df <- data.frame(vx = adv[1], vy = adv[2])
} else {
  wind_df <- NA
  fixed_eta1 <- NA
  fixed_eta2 <- NA
  eta_type <- ""
}

# Folder name to save the data
foldername <- file.path(
  paste0(data_folder, "simulations_rpar/coords"),
  paste0("rpar_", param_str),
  paste0("sim_", ngrid^2, "s_", length(temp), "t_s0_", s0_str, "_t0_", t0_str),
  s0_type,
  wind_type,
  distance_type,
  eta_type
)

if (!dir.exists(foldername)) {
  message("Folder created: ", foldername)
  dir.create(foldername, recursive = TRUE)
}

tau_vect <- 0:10
u <- 1
sites_coords <- generate_grid_coords(ngrid)
sites_coords <- as.data.frame(sites_coords)

# Parallel configuration
num_cores <- detectCores() - 1
cl <- makeCluster(num_cores)

# Export variables and functions to the cluster
# Export variables and functions to the cluster
clusterExport(cl, varlist = c(
  "sim_rpareto_coords", "beta1", "beta2", "alpha1", "alpha2",
  "spa", "temp", "adv_real", "t0", "m", "random_s0", "s0",
  "s0_radius", "distance_type",
  "foldername", "generate_grid_coords", "save_simu",
  "ngrid", "compute_st_gaussian_process_coords",
  "compute_st_variogram_coords", "convert_simulations_to_list",
  "empirical_excesses_rpar", "get_conditional_lag_vectors",
  "t0", "u", "tau_vect", "sites_coords"
))

# Load libraries in the cluster
clusterEvalQ(cl, {
  library(RandomFields)
  library(dplyr)
  library(tidyr)
  library(data.table)
  RandomFields::RFoptions(
    spConform = FALSE,
    allow_duplicated_locations = TRUE,
    storing = FALSE,
    printlevel = 0
  )
})

# Simulate the process in parallel
result_list <- parLapply(cl, 1:M, function(i) {

  simu <- sim_rpareto_coords(
    beta1 = beta1,
    beta2 = beta2,
    alpha1 = alpha1,
    alpha2 = alpha2,
    coords = sites_coords,
    t = temp,
    adv = adv_real,
    t0 = t0,
    nres = m,
    random_s0 = FALSE,
    s0 = s0,
    s0_radius = s0_radius,
    distance = distance_type
  )
  
  list_s0 <- simu$s0_used
  list_rpar <- save_simu(simu$Z, sites_coords, folder = foldername)

  list_lags_excesses <- lapply(1:m, function(j) {
    s0_df <- list_s0[[j]]
    s0_x <- s0_df$x
    s0_y <- s0_df$y
    s0_coords <- sites_coords[sites_coords$Longitude == s0_x &
                              sites_coords$Latitude == s0_y, ]
    lags <- get_conditional_lag_vectors(sites_coords, s0_coords, t0,
                                        tau_vect, latlon = FALSE)
    excesses <- empirical_excesses_rpar(list_rpar[[j]], quantile = u,
                                        df_lags = lags,
                                        threshold = TRUE, t0 = t0)
    list(lags = lags, excesses = excesses)
  })

  # Extract lags and excesses from the list of results
  list_lags <- lapply(list_lags_excesses, `[[`, "lags")
  list_excesses <- lapply(list_lags_excesses, `[[`, "excesses")

  return(list(
    list_rpar = list_rpar,
    s0_used = simu$s0_used,
    list_excesses = list_excesses,
    list_lags = list_lags
  ))
})

stopCluster(cl) # Stop the cluster

# Get the list of data.frames
list_rpar <- do.call(c, lapply(result_list, function(res) res$list_rpar))
list_excesses <- do.call(c, lapply(result_list,
                                   function(res) res$list_excesses))
list_lags <- do.call(c, lapply(result_list, function(res) res$list_lags))
list_s0 <- do.call(c, lapply(result_list, function(res) res$s0_used))

### Optimization ###############################################################

# Initialization
if (use_wind_data) { # with eta1 and eta2
  true_param <- c(beta1, beta2, alpha1, alpha2, eta1, eta2)
} else { # with adv1 and adv2
  true_param <- c(beta1, beta2, alpha1, alpha2, adv[1], adv[2])
}

init_diff <- FALSE

if (init_diff) {
  # Get the initial values for the parameters with a small difference
  true_param <- true_param + c(0.1, 0.1, 0.1, 0.1, 0.1, 0.1)
  init_type <- "init_diff"
} else {
  init_type <- ""
}


result_list <- mclapply(1:M, process_simulation, m = m,
                        list_simu = list_rpar, u = u,
                        list_lags = list_lags,
                        list_excesses = list_excesses,
                        init_params = true_param,
                        distance = distance_type,
                        hmax = 7, wind_df = wind_df,
                        mc.cores = num_cores)


# Combine results int a data frame
df_result_all <- do.call(rbind, result_list)
df_result_all <- as.data.frame(df_result_all)
# Redefine the column names
if(use_wind_data) {
  colnames(df_result_all) <- c("beta1", "beta2", "alpha1",
                             "alpha2", "eta1", "eta2")
  # Add adv1 and adv2 with etas estimated
  df_result_all$adv1 <- df_result_all$eta1 * abs(wind_df$vx)^df_result_all$eta2 * sign(wind_df$vx)
  df_result_all$adv2 <- df_result_all$eta1 * abs(wind_df$vy)^df_result_all$eta2 * sign(wind_df$vy)
} else {
  colnames(df_result_all) <- c("beta1", "beta2", "alpha1",
                             "alpha2", "adv1", "adv2")
}

# Eta type (if it has changed in optimization)
eta_type <- if (!is.na(fixed_eta1) && !is.na(fixed_eta2)) {
  "fixed_eta"
} else if (!is.na(fixed_eta1)) {
  "fixed_eta1"
} else if (!is.na(fixed_eta2)) {
  "fixed_eta2"
} else {
  "free_eta"
}

# save data in csv
foldername_result <- file.path(
  data_folder,
  "optim_results/coords/rpar",
  paste0("rpar_", param_str),
  paste0("sim_", ngrid^2, "s_", length(temp), "t_s0_", s0_str, "_t0_", t0_str),
  s0_type,
  wind_type,
  distance_type,
  eta_type,
  init_type
)
if (!dir.exists(foldername_result)) {
  print(paste0("Folder created: ", foldername_result))
  dir.create(foldername_result, recursive = TRUE)
}

name_file <- paste0("/optim_rpar_", M, "simu_", m, "rep_", ngrid^2,
                "s_", length(temp), "t_", param_str, ".csv")

# Save the results as a CSV file
write.csv(df_bplot, paste0(foldername_result, name_file), row.names = FALSE)

# Get the number of simulations without convergence
number_no_convergence <- sum(is.na(df_result_all$beta1))
# Save the number into a text file
if (number_no_convergence > 0) {
  filename <- paste0(
    "/no_convergence_", M, "simu_", m, "rep_", ngrid^2,
    "s_", length(temp), "t_", param_str, ".txt"
  )
  text <- paste0(
    "Number of simulations without convergence: ", number_no_convergence,
    "\n",
    "Total number of simulations: ", M
  )
  writeLines(text, paste0(foldername_result, filename))
}

# Get data frame with the results for boxplot
# Convert the data frame to long format
df_bplot <- stack(df_result_all)

# Latex labels for the parameters
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

# Get the initial parameters for true values
initial_param <- c(
  beta1 = beta1,
  beta2 = beta2,
  alpha1 = alpha1,
  alpha2 = alpha2,
  eta1 = 2,
  eta2 = 1,
  adv1 = adv_real[1],
  adv2 = adv_real[2]
)

# First plot
params_group1 <- c("beta1", "beta2", "alpha1", "alpha2")

# Second plot
params_group2 <- if (use_wind_data) {
  c("eta1", "eta2", "adv1", "adv2")
} else {
  c("adv1", "adv2")
}

library(ggplot2)
# Create the boxplots
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

# folder to save the plots
foldername_image <- file.path(
  im_folder,
  "/optim_results/coords/rpar",
  paste0("rpar_", param_str),
  paste0("sim_", ngrid^2, "s_", length(temp), "t_s0_", s0_str, "_t0_", t0_str),
  s0_type,
  wind_type,
  distance_type,
  eta_type,
  init_type
)

if (!dir.exists(foldername_image)) {
  print(paste0("Folder created: ", foldername_image))
  dir.create(foldername_image, recursive = TRUE)
}

# Save the plots
name_file <- paste0("/bp_optim_", M, "simu_", m, "rep_", ngrid^2,
                "s_", length(temp), "t_", param_str, "_beta_alpha.png")

ggsave(plot = bplot1, filename = paste0(foldername_image, name_file),
       width = 5, height = 5)

suffix <- if (all(!is.na(wind_df))) "etas" else "adv"

name_file <- paste0( "/bp_optim_", M, "simu_", m, "rep_", ngrid^2,
              "s_", length(temp), "t_", param_str, "_", suffix, ".png")

ggsave(plot = bplot2, filename = paste0(foldername_image, name_file),
       width = 5, height = 5)
