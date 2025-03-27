
# remove objects from the environment
rm(list = ls())

# clear console
cat("\014")

# LOAD LIBRARIES ###############################################################
muse <- FALSE

if (muse) {
  # Get the muse folder
  folder_muse <- "/home/serrec/work_rainstsimu/downscaling"
  setwd(folder_muse)
  # Load libraries and set theme
  source("load_libraries.R")
  source("pinnEV.R")
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
M <- 100 # number of simulations
m <- 1000 # number of extreme episodes
nres <- M * m

# Simulate the process
set.seed(123)
simu <- sim_rpareto(beta1, beta2, alpha1, alpha2, spa, spa, temp, adv, s0,
                    t0, nres)

# Function to format values correctly, considering decimal places
format_value <- function(x) {
  # Check if x is an integer
  if (x == as.integer(x)) {
    return(sprintf("%d", x))
  } else {
    # Count number of decimal places
    num_decimals <- nchar(sub("^[^.]*\\.", "", as.character(x)))

    if (x >= 1) {
      # If the number is greater than or equal to 1
      if (num_decimals == 1) {
        return(sprintf("%02d", round(x * 10)))
      } else {
        return(sprintf("%03d", round(x * 100)))
      }
    } else {
      # If the number is less than 1
      if (num_decimals == 1) {
        return(sprintf("%02d", round(x * 10)))
      } else {
        return(sprintf("%03d", round(x * 100)))
      }
    }
  }
}
# Apply formatting
params_str <- paste(sapply(params, format_value), collapse = "_")
adv_str <- paste(sapply(adv, format_value), collapse = "_")
s0_str <- paste(sapply(s0, format_value), collapse = "_")
t0_str <- format_value(t0)

# Save the data
foldername <- paste0("./data/simulations_rpar/rpar_", params_str, "_", adv_str,
                   "/sim_", ngrid^2, "s_", length(temp), "t_s0_",
                    s0_str, "_t0_", t0_str, "/")

if (!dir.exists(foldername)) {
  dir.create(foldername, recursive = TRUE)
}

save_simulations(simu, ngrid, nres, folder = foldername,
        file = paste0("rpar_", ngrid^2, "s_", length(temp), "t"))

list_simuM <- list()
for (i in 1:nres) {
  file_name <- paste0(foldername, "rpar_", ngrid^2, "s_",
                        length(temp), "t_", i, ".csv")
  list_simuM[[i]] <- read.csv(file_name)
}


### Optimization ###############################################################
library(parallel)
num_cores <- detectCores() - 1  # Reserve 1 core for the OS

# Parallel execution
sites_coords <- generate_grid_coords(ngrid)
df_lags <- get_conditional_lag_vectors(sites_coords, s0, t0, tau_vect = 0:10)
u <- 1 # threshold corresponding to the r-pareto simulation
result_list <- mclapply(1:M, process_simulation, M = M, m = m,
                        list_simuM = list_simuM, u = u, df_lags = df_lags,
                        t0 = t0, true_param = true_param, hmax = sqrt(17),
                        mc.cores = num_cores)

# Combine results into a data frame
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

