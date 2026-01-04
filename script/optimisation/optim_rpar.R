
# remove objects from the environment
rm(list = ls())

# clear console
cat("\014")

# LOAD LIBRARIES ###############################################################
library(parallel)
library(abind)
muse <- FALSE
# PARAMETERS ###################################################################


if (muse) {
  # Get the muse folder
  folder_muse <- "/home/serrec/work_rainstsimu/rpareto/"
  setwd(folder_muse)
  # Load libraries and set theme
  source("load_libraries.R")
  im_folder <- "./images"
  source("config.R")
  data_folder <- "./data/"
  num_cores <- 27
} else {
  # Load libraries and set theme
  num_cores <- detectCores() - 1
  source("./script/load_libraries.R")
  source("./script/optimisation/config.R")
  M <- num_cores
}

# Get all files in the folder "R"
functions_folder <- "./R"
files <- list.files(functions_folder, full.names = TRUE)
# load all functions in files
invisible(lapply(files, function(f) source(f, echo = FALSE)))
library(latex2exp)
print("All functions loaded")

# SIMULATION ###################################################################
result_folder <- paste0(data_folder, "optim_results/rpar/")
if (!dir.exists(result_folder)) {
  dir.create(result_folder, recursive = TRUE)
}


# get random advection for each simulation
set.seed(123)
adv_df <- data.frame(
  adv1 = rnorm(m, mean = 1, sd = 0.5),
  adv2 = rnorm(m, mean = 1, sd = 0.2)
)

set.seed(123)

wind_list <- lapply(1:M, function(i) {
  data.frame(
    vx = rnorm(m, mean = 1, sd = 0.5),
    vy = rnorm(m, mean = 1, sd = 0.2)
  )
})


# params <- c(0.4, 0.7, 0.3, 0.8) # true parameters for the variogram
# eta1 <- 0.5
# eta2 <- 1.5
true_param <- c(params, eta1, eta2)
beta1 <- params[1]
beta2 <- params[2]
alpha1 <- params[3]
alpha2 <- params[4]
spa <- 1:ngrid

nsites <- ngrid^2 # if the grid is squared

# Number of realizations
nres <- M * m

# Apply formatting
format_value <- function(x) {
  formatted_values <- sapply(x, function(val) {
    # Check if val is negative
    is_negative <- val < 0
    val <- abs(val)  # Work with the absolute value for formatting

    # Check if val is an integer
    if (val == as.integer(val)) {
      formatted_value <- sprintf("%d", val)
    } else {
      # Count the number of decimals
      num_decimals <- nchar(sub("^[^.]*\\.", "", as.character(val)))

      if (val >= 1) {
        # If the number is greater than or equal to 1
        if (num_decimals == 1) {
          formatted_value <- sprintf("%02d", round(val * 10))
        } else {
          formatted_value <- sprintf("%03d", round(val * 100))
        }
      } else {
        # If the number is less than 1
        if (num_decimals == 1) {
          formatted_value <- sprintf("%02d", round(val * 10))
        } else {
          formatted_value <- sprintf("%03d", round(val * 100))
        }
      }
    }

    # Add "neg" prefix if the number was negative
    if (is_negative) {
      formatted_value <- paste0("neg", formatted_value)
    }

    return(formatted_value)
  })

  # Concatenate values
  return(paste(formatted_values, collapse = "_"))
}
param_str <- format_value(true_param)
# adv_str <- format_value(adv)
s0_str <- format_value(s0)
t0_str <- format_value(t0)

# Subfolder names
s0_type <- if (random_s0) "random_s0" else "fixed_s0"
fixed_eta1 <- eta1
fixed_eta2 <-eta2
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
use_wind_data <- TRUE
if (use_wind_data) {
  wind_df <- adv_df
  fixed_eta1 <- eta1
  fixed_eta2 <- eta2
} else {
  wind_df <- NA
  fixed_eta1 <- NA
  fixed_eta2 <- NA
  eta_type <- ""
}

# Folder name to save the data
foldername <- file.path(
  result_folder,
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

radial_adv_from_wind <- function(wind_df, eta1, eta2, eps = 1e-12) {
  vx <- wind_df$vx
  vy <- wind_df$vy
  r <- sqrt(vx^2 + vy^2)
  scale <- eta1 * (pmax(r, eps)^(eta2 - 1))
  cbind(scale * vx, scale * vy)
}


# Parallel configuration
print("Starting simulations...")

cl <- makeCluster(num_cores)

# Export variables and functions to the cluster
clusterExport(cl, varlist = c(
  "sim_rpareto", "beta1", "beta2", "alpha1", "alpha2",
  "spa", "temp", "adv_df", "t0", "m", "random_s0", "s0",
  "s0_radius", "distance_type",
  "foldername", "generate_grid_coords", "save_simulations",
  "ngrid", "compute_st_gaussian_process",
  "compute_st_variogram", "convert_simulations_to_list",
  "empirical_excesses_rpar", "get_conditional_lag_vectors",
  "t0", "u", "tau_vect", "sites_coords", "eta1", "eta2",
  "radial_adv_from_wind", "wind_list"
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

  wind_df_i <- wind_list[[i]]
  adv_mat_i <- radial_adv_from_wind(wind_df_i, eta1, eta2)  # m x 2

  simu <- sim_rpareto(
    beta1 = beta1, beta2 = beta2,
    alpha1 = alpha1, alpha2 = alpha2,
    x = spa, y = spa, t = temp,
    adv = adv_mat_i,
    t0 = t0,
    nres = m,
    random_s0 = random_s0,
    s0 = s0,
    s0_radius = s0_radius,
    distance = distance_type
  )


  list_rpar <- convert_simulations_to_list(simu$Z, sites_coords)

  list_lags_excesses <- lapply(1:m, function(j) {
    s0_x <- simu$s0_used[[j]][[1]]$x
    s0_y <- simu$s0_used[[j]][[1]]$y
    s0_coords <- sites_coords[sites_coords$Longitude == s0_x &
                              sites_coords$Latitude == s0_y, ]
    lags <- get_conditional_lag_vectors(sites_coords, s0_coords, t0,
                                        tau_vect, latlon = FALSE)
    excesses <- empirical_excesses_rpar(list_rpar[[j]],
                                        df_lags = lags,
                                        t0 = t0, threshold = u)
    lags$tau <- lags$tau
    list(lags = lags, excesses = excesses)
  })

  # Extract lags and excesses from the list of results
  list_lags <- lapply(list_lags_excesses, `[[`, "lags")
  list_excesses <- lapply(list_lags_excesses, `[[`, "excesses")

  return(list(
    list_rpar = list_rpar,
    s0_used = simu$s0_used,
    list_excesses = list_excesses,
    list_lags = list_lags,
    wind_df = wind_df_i
  ))
})
stopCluster(cl) # Stop the cluster

print("Simulations done")
# Get the list of data.frames
list_rpar <- do.call(c, lapply(result_list, function(res) res$list_rpar))
list_excesses <- do.call(c, lapply(result_list,
                                   function(res) res$list_excesses))
list_lags <- do.call(c, lapply(result_list, function(res) res$list_lags))

list_wind <- lapply(result_list, function(res) res$wind_df)


# check excesses
sums_kij <- sapply(list_excesses, function(ex) sum(ex$kij, na.rm = TRUE))
sums_kij

### Optimization ###############################################################

# Initialization
if (use_wind_data) { # with eta1 and eta2
  true_param <- c(beta1, beta2, alpha1, alpha2, eta1, eta2)
  wind_df <- adv_df
} else { # with adv1 and adv2
  true_param <- c(beta1, beta2, alpha1, alpha2, adv[1], adv[2])
  wind_df <- NA
}

init_diff <- FALSE

if (init_diff) {
  # Get the initial values for the parameters with a small difference
  true_param <- true_param + c(0.1, 0.1, 0.1, 0.1, 0.1, 0.1)
  init_type <- "init_diff"
} else {
  init_type <- ""
}

init_params <- true_param

result_list_opt <- mclapply(
  1:M,
  function(i) {
    process_simulation(
      i = i, m = m,
      list_simu = list_rpar,
      u = u,
      list_lags = list_lags,
      list_excesses = list_excesses,
      init_params = init_params,
      distance = distance_type,
      hmax = 7,
      wind_df = list_wind[[i]],
      fixed_eta1 = fixed_eta1,
      fixed_eta2 = fixed_eta2,
      normalize = FALSE
    )
  },
  mc.cores = num_cores
)


print("Optimization done")
# Combine results_list into a data.frame
df_result_all <- do.call(rbind, result_list_opt)

# Redefine the column names
if (use_wind_data) {
  colnames(df_result_all) <- c("beta1", "beta2", "alpha1",
                             "alpha2", "eta1", "eta2")
  
} else {
  colnames(df_result_all) <- c("beta1", "beta2", "alpha1",
                             "alpha2", "adv1", "adv2")
}

print(head(df_result_all))

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
  "./data/optim_results/rpar",
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

write.csv(df_result_all, paste0(foldername_result, name_file), row.names = FALSE)
print(paste0("Results saved in: ", foldername_result, name_file))


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
} else {
  print("All simulations converged")
}

# PLOTTING RESULTS #############################################################
df_bplot <- stack(df_result_all)
# write.csv(df_bplot, paste0(foldername_result, name_file), row.names = FALSE)

# Latex labels for the parameters
labels_latex <- c(
  beta1 = TeX("$\\beta_1$"),
  beta2 = TeX("$\\beta_2$"),
  alpha1 = TeX("$\\alpha_1$"),
  alpha2 = TeX("$\\alpha_2$"),
  eta1 = TeX("$\\eta_1$"),
  eta2 = TeX("$\\eta_2$")
)

# Get the initial parameters for true values
initial_param <- c(
  beta1 = beta1,
  beta2 = beta2,
  alpha1 = alpha1,
  alpha2 = alpha2,
  eta1 = 1,
  eta2 = 1
)

# First plot
params_group1 <- c("beta1", "beta2", "alpha1", "alpha2")

# Second plot
params_group2 <- if (use_wind_data) {
  c("eta1", "eta2")
}

library(ggplot2)
# Create the boxplots
bplot1 <- ggplot(subset(df_bplot, ind %in% params_group1), aes(x = ind, y = values)) +
  geom_boxplot() +
  geom_point(aes(y = initial_param[match(ind, names(initial_param))]), color = "red", pch = 4) +
  scale_x_discrete(labels = labels_latex[params_group1]) +
  labs(x = "Parameters", y = "Estimated values") +
  theme_minimal()


# folder to save the plots
foldername_image <- file.path(
  im_folder,
  "/optim_results/rpar",
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
                "s_", length(temp), "t_", param_str, "_beta_alpha.pdf")

ggsave(plot = bplot1, filename = paste0(foldername_image, name_file),
       width = 5, height = 5)