# remove objects from the environment
rm(list = ls())

# clear console
cat("\014")

# Load libraries and set theme
source("./script/load_libraries.R")

# library(generain)
library(ggplot2)
library(reshape2)
library(animation)
library(RandomFields)
library(RandomFieldsUtils)

# Get all files in the folder "R"
functions_folder <- "./R"
files <- list.files(functions_folder, full.names = TRUE)
# load all functions in files
invisible(lapply(files, function(f) source(f, echo = FALSE)))
library(latex2exp)


# get grid coordinates
filename <- paste0(data_folder, "comephore/grid/grid_coords_5km.csv")
grid_comephore <- read.csv(filename, stringsAsFactors = FALSE)
colnames(grid_comephore) <- c("Site", "Longitude", "Latitude")

################################################################################
# Simulate data using the rpareto model and estimated parameters from 
# WLSE on the comephore dataset
################################################################################

ngrid <- 10
spa <- 1:ngrid
temp <- 0:29
m <- 800 # number of episodes
M <- 1 # number of simulations
n.res <- m * M
param <- c(0.0261, 0.5815, 1.1321, 0.6741) # WLSE estimated parameters
threshold <- 12 # threshold for the r-Pareto process
beta1 <- param[1]
beta2 <- param[2]
alpha1 <- param[3]
alpha2 <- param[4]
adv <- c(0, 0) # separable model with WLSE parameters
s0 <- c(1, 1)
s0_radius <- 7
t0 <- 0
random_s0 <- TRUE

# Simulate spatio-temporal r-Pareto process
simu_rpar <- sim_rpareto(param[1], param[2], param[3], param[4], spa, spa, temp,
                         adv, t0, n.res, random_s0, s0)

simu_rpar <- sim_rpareto_coords(alpha1 = alpha1, alpha2 = alpha2,
                                beta1 = beta1, beta2 = beta2,
                                coords = grid_comephore,
                                t = temp, adv = adv, t0 = t0, nres = n.res,
                                random_s0 = random_s0, s0 = s0,
                                s0_radius = Inf, distance = "euclidean",
                                threshold = threshold)

params <- c(param, adv)
# Formatage des paramètres pour la sauvegarde
param_str <- format_value(params)
adv_str <- format_value(adv)
s0_str <- format_value(s0)
t0_str <- format_value(t0)

if (random_s0) {
  s0_str <- "random_s0"
} else {
  s0_str <- paste0("s0_", s0_str)
}


Z_rpar <- simu_rpar$Z
s0_list <- simu_rpar$s0_used

foldername <- paste0(data_folder, "comephore/simulations/rpar_",
                                length(temp), "t/")
# Create folder if it does not exist
if (!dir.exists(foldername)) {
  dir.create(foldername, recursive = TRUE)
}

save_simu(Z_rpar, grid_comephore, folder = foldername)


# verify for all episode that Z_s0,t0 > threshold
check_Z_gtu <- sapply(seq_along(s0_list), function(i) {
  s0_used <- s0_list[[i]]
  s0_index <- which(grid_comephore$Site == s0_used$Site)
  all(Z_rpar[s0_index, t0 + 1, i] > threshold)
})

all(check_Z_gtu)


# Get data from file
files <- list.files(foldername, full.names = TRUE)
length(files)
list_rpar <- list()
for (j in 1:m) {
    file_name <- paste0(foldername, "episode_", j, ".csv")
    list_rpar[[j]] <- read.csv(file_name)
}

simu_df <- list_rpar[[1]] # first simulation
s0_used <- s0_list[[1]]
nsites <- ncol(simu_df) # number of sites
par(mfrow = c(1, 1))
plot(simu_df[, 1], main = "rpareto simulation")
# get grid coordinates
gif_folder <- "/user/cserreco/home/Documents/These/phd_extremes/thesis/resources/images/swg/comephore/rpar"
if (!dir.exists(gif_folder)) {
  dir.create(gif_folder, recursive = TRUE)
}
gif_folder <- paste0(gif_folder, "/", s0_str, "/")
create_simu_gif(simu_df, sites_coords, c(param, adv), type = "rpar",
                foldername = gif_folder, s0 = s0_used)

length(list_rpar)
# concat all simulations together
simu_all <- do.call(rbind, list_rpar)
create_simu_gif(simu_all, grid_comephore, c(param, adv), type = "rpar",
                foldername = gif_folder, forcedtemp = 200, threshold = threshold)


