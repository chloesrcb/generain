
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
  random_s0 <- TRUE
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
tau_vect <- 0:10
u <- 1

nsites <- ngrid^2 # if the grid is squared

# Number of realizations
nres <- M * m

# Apply formatting
param_str <- format_value(true_param)
adv_str <- format_value(adv)
s0_str <- format_value(s0)
t0_str <- format_value(t0)

# Save the data
foldername <- paste0("./data/simulations_rpar/rpar_test_", param_str,
                   "/sim_", ngrid^2, "s_", length(temp), "t_s0_",
                    s0_str, "_t0_", t0_str, "/")

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
  print(paste0("Folder created: ", foldername))
  dir.create(foldername, recursive = TRUE)
}

# CONFIG
tmax <- max(tau_vect)
folder_save <- "./data/processed_batches/"
dir.create(folder_save, recursive = TRUE, showWarnings = FALSE)

# Préparation du cluster
num_cores <- detectCores() - 1
cl <- makeCluster(num_cores)

clusterExport(cl, varlist = c(
  "sim_rpareto", "beta1", "beta2", "alpha1", "alpha2",
  "spa", "temp", "adv", "t0", "m", "random_s0", "s0",
  "s0_radius", "is_anisotropic", "folder_save",
  "foldername", "ngrid", "u", "tau_vect"
))

clusterEvalQ(cl, {
  # Charger tous les fichiers de ./R dans chaque worker
  functions_folder <- "./R"
  files <- list.files(functions_folder, full.names = TRUE)
  invisible(lapply(files, function(f) source(f, echo = FALSE)))

  # Charger les packages aussi ici si nécessaire
  library(RandomFields)
  library(dplyr)
  library(tidyr)
  RandomFields::RFoptions(spConform = FALSE)
})


# Fonction pour traiter un batch
simulate_and_save_batch <- function(i) {
  cat("Batch", i, "\n")
  simu <- sim_rpareto(
    beta1, beta2, alpha1, alpha2,
    x = 1:ngrid, y = 1:ngrid, t = temp,
    adv = adv, t0 = t0, nres = m,
    random_s0 = random_s0, s0 = s0,
    s0_radius = s0_radius,
    anisotropic = is_anisotropic
  )

  list_rpar <- convert_simulations_to_list(simu$Z, ngrid)
  s0_list <- simu$s0_used
  sites_coords <- as.data.frame(generate_grid_coords(ngrid))

  list_results <- lapply(1:m, function(j) {
    s0_coords <- sites_coords[
      sites_coords$Longitude == s0_list[[j]]$x &
      sites_coords$Latitude == s0_list[[j]]$y, ]
    
    episode <- list_rpar[[j]]
    lags <- get_conditional_lag_vectors(sites_coords, s0_coords, t0,
                                        tau_vect, latlon = FALSE)
    excesses <- empirical_excesses(episode, u, lags, type = "rpareto",
                                   t0 = t0, threshold = TRUE)
    list(lags = lags, excesses = excesses, episode = episode)
  })

  saveRDS(list_results, file = paste0(folder_save, "batch_", i, ".rds"))
  gc()
  return(NULL)
}

# Lancement parallèle
parLapply(cl, 1:M, simulate_and_save_batch)
stopCluster(cl)


# Load all batches
cat("Loading all batches...\n")
all_batches <- list()
for (i in 1:M) {
  batch <- readRDS(paste0(folder_save, "batch_", i, ".rds"))
  all_batches <- c(all_batches, batch)
}

cat("Extracting components...\n")
list_lags      <- lapply(all_batches, `[[`, "lags")
list_excesses  <- lapply(all_batches, `[[`, "excesses")
list_episodes  <- lapply(all_batches, `[[`, "episode")

# Optional: reconstruct full data structure for your optimizer
# This depends on how process_simulation() expects data

if (is_anisotropic) {
   directional <- TRUE
} else {
   directional <- FALSE
}

# Full global optimization without wind data
cat("Running optimization...\n")
result_list <- mclapply(1:M, function(i) {
  process_simulation(
    i = i,
    m = m,
    list_simu = list_episodes,
    u = u,
    list_lags = list_lags,
    list_excesses = list_excesses,
    init_params = true_param,
    directional = directional,
    hmax = 7,
    wind_df = NA
  )
}, mc.cores = detectCores() - 1)

# Combine results
df_result_all <- do.call(rbind, result_list)
colnames(df_result_all) <- c("beta1", "beta2", "alpha1",
                             "alpha2", "adv1", "adv2")

# Combine results int a data frame
df_result_all <- do.call(rbind, result_list)
colnames(df_result_all) <- c("beta1", "beta2", "alpha1",
                             "alpha2", "adv1", "adv2")

df_bplot <- as.data.frame(df_result_all)

# save data in csv
foldername <- "./data/optim/rpar/"
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


# With wind data == adv to estimate eta1 and eta2 to 1
cat("Running optimization with wind...\n")
wind_df <- data.frame(vx = adv[1], vy = adv[2])
result_list <- mclapply(1:M, function(i) {
  process_simulation(
    i = i,
    m = m,
    list_simu = list_episodes,
    u = u,
    list_lags = list_lags,
    list_excesses = list_excesses,
    init_params = true_param,
    directional = directional,
    hmax = 7,
    wind_df = wind_df
  )
}, mc.cores = detectCores() - 1)

# Combine results int a data frame
df_result_all <- as.data.frame(do.call(rbind, result_list))
colnames(df_result_all) <- c("beta1", "beta2", "alpha1",
                             "alpha2", "eta1", "eta2")


df_result_all$adv1 <- df_result_all$eta1 * abs(adv[1])^df_result_all$eta2 * sign(adv[1])
df_result_all$adv2 <- df_result_all$eta1 * abs(adv[2])^df_result_all$eta2 * sign(adv[2])


# Convert in latex
df_bplot <- as.data.frame(df_result_all)

# save data in csv
foldername <- "./data/optim/rpar/eta/"
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
}

name_file <- paste0("optim_rpar_", M, "simu_", m, "rep_", ngrid^2,
                "s_", length(temp), "t_", param_str, ".csv")
write.csv(df_bplot, paste0(foldername, name_file), row.names = FALSE)

# Plot results
df_bplot <- stack(df_bplot)

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

params_group1 <- c("beta1", "beta2", "alpha1", "alpha2")
params_group2 <- c("eta1", "eta2", "adv1", "adv2")

true_param <- c(
  beta1 = beta1,
  beta2 = beta2,
  alpha1 = alpha1,
  alpha2 = alpha2,
  eta1 = 1,
  eta2 = 1,
  adv1 = adv[1],
  adv2 = adv[2]
)

bplot1 <- ggplot(subset(df_bplot, ind %in% params_group1), aes(x = ind, y = values)) +
  geom_boxplot() +
  geom_point(aes(y = true_param[match(ind, names(true_param))]), color = "red", pch = 4) +
  scale_x_discrete(labels = labels_latex[params_group1]) +
  labs(x = "Parameters", y = "Estimated values") +
  theme_minimal()

bplot2 <- ggplot(subset(df_bplot, ind %in% params_group2), aes(x = ind, y = values)) +
  geom_boxplot() +
  geom_point(aes(y = true_param[match(ind, names(true_param))]), color = "red", pch = 4) +
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
name_file <- paste0("bp_optim_eta_", M, "simu_", m, "rep_", ngrid^2,
                "s_", length(temp), "t_", param_str, "_beta_alpha.png")

ggsave(plot = bplot1, filename = paste0(foldername, name_file),
       width = 5, height = 5)

# Save the plot
name_file <- paste0("bp_optim_eta_", M, "simu_", m, "rep_", ngrid^2,
              "s_", length(temp), "t_", param_str, "_eta_adv.png")

ggsave(plot = bplot2, filename = paste0(foldername, name_file),
       width = 5, height = 5)