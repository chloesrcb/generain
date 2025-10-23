
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
  num_cores <- 27
} else {
  # Load libraries and set theme
  num_cores <- detectCores() - 1
  source("./script/load_libraries.R")
  source("./script/optimisation/config.R")
}

# Get all files in the folder "R"
functions_folder <- "./R"
files <- list.files(functions_folder, full.names = TRUE)
# load all functions in files
invisible(lapply(files, function(f) source(f, echo = FALSE)))
library(latex2exp)
print("All functions loaded")

# SIMULATION ###################################################################
result_folder <- paste0(data_folder, "simulations_rpar/")
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
  paste0(data_folder, "simulations_rpar"),
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
print("Starting simulations...")

cl <- makeCluster(num_cores)

# Export variables and functions to the cluster
clusterExport(cl, varlist = c(
  "sim_rpareto", "beta1", "beta2", "alpha1", "alpha2",
  "spa", "temp", "adv_real", "t0", "m", "random_s0", "s0",
  "s0_radius", "distance_type",
  "foldername", "generate_grid_coords", "save_simulations",
  "ngrid", "compute_st_gaussian_process",
  "compute_st_variogram", "convert_simulations_to_list",
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
  simu <- sim_rpareto(
    beta1 = beta1,
    beta2 = beta2,
    alpha1 = alpha1,
    alpha2 = alpha2,
    x = spa,
    y = spa,
    t = temp,
    adv = adv_real,
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

print("Simulations done")
# Get the list of data.frames
list_rpar <- do.call(c, lapply(result_list, function(res) res$list_rpar))
list_excesses <- do.call(c, lapply(result_list,
                                   function(res) res$list_excesses))
list_lags <- do.call(c, lapply(result_list, function(res) res$list_lags))

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

# Optimization
init_params <- c(params[1:4], 0, 0)
result_list <- mclapply(1:M, process_simulation, m = m,
                        list_simu = list_rpar, u = u,
                        list_lags = list_lags,
                        list_excesses = list_excesses,
                        init_params = init_params,
                        distance = distance_type,
                        hmax = 7, wind_df = wind_df,
                        mc.cores = num_cores)



print("Optimization done")
# Combine results int a data frame
df_result_all <- do.call(rbind, result_list)
df_result_all <- as.data.frame(df_result_all)
# Redefine the column names
if (use_wind_data) {
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


# Theoretical variogram
variogram <- function(lags_df, beta1, beta2, alpha1, alpha2, adv,
                      distance_type = "euclidean") {
  df <- as.data.frame(lags_df)
  tau <- df$tau

  if (distance_type == "lalpha") {
    df$gamma <- beta1 * abs(df$hxV)^alpha1 +
                beta1 * abs(df$hyV)^alpha1 +
                beta2 * abs(tau)^alpha2
  } else {
    df$gamma <- beta1 * sqrt(df$hxV^2 + df$hyV^2)^alpha1 +
                beta2 * abs(tau)^alpha2
  }

  return(df)
}


lag_all <- do.call(rbind, list_lags)
df_lags_all <- as.data.frame(lag_all)

# keep only unique rows
df_lags_all <- unique(df_lags_all)

# Compute hnorm (spatial lag norm) if not present
df_lags_all$hxV <- df_lags_all$hx - adv_real[1] * df_lags_all$tau
df_lags_all$hyV <- df_lags_all$hy - adv_real[2] * df_lags_all$tau
df_lags_all$hnormV <- sqrt(df_lags_all$hxV^2 + df_lags_all$hyV^2)



# Plot all empirical variogram curves for each tau, with one theoretical curve
for (tau in tau_vect) {
  # Filter lags for current tau
  df_lags_tau <- df_lags_all[df_lags_all$tau == tau, ]
  vario_th_tau <- variogram(df_lags_tau, beta1, beta2, alpha1, alpha2, adv_real, 
                            distance_type)
  vario_th_tau <- vario_th_tau %>%
    group_by(hnormV) %>%
    summarise(gamma = mean(gamma), .groups = "drop")
  
  # Prepare empirical curves for all simulations
  vario_emp_all <- lapply(1:nrow(df_result_all), function(j) {
    beta1_hat <- df_result_all$beta1[j]
    beta2_hat <- df_result_all$beta2[j]
    alpha1_hat <- df_result_all$alpha1[j]
    alpha2_hat <- df_result_all$alpha2[j]
    if (use_wind_data) {
      eta1_hat <- df_result_all$eta1[j]
      eta2_hat <- df_result_all$eta2[j]
      adv1_hat <- eta1_hat * abs(wind_df$vx)^eta2_hat * sign(wind_df$vx)
      adv2_hat <- eta1_hat * abs(wind_df$vy)^eta2_hat * sign(wind_df$vy)
      adv_hat <- c(adv1_hat, adv2_hat)
    } else {
      adv_hat <- c(df_result_all$adv1[j], df_result_all$adv2[j])
    }
    vario_emp_j <- variogram(df_lags_tau, beta1_hat, beta2_hat, alpha1_hat,
                              alpha2_hat, adv_hat, distance_type)
    vario_emp_j$sim <- j
    vario_emp_j
  })
  vario_emp_all_df <- bind_rows(vario_emp_all)

  # get mean empirical value for each hnorm
  vario_emp_all_df <- vario_emp_all_df %>%
    group_by(sim, hnormV) %>%
    summarise(gamma = mean(gamma), .groups = "drop")
  
  vario_emp_all_df$type <- "Empirical"
  vario_th_tau$type <- "Theoretical"
  vario_plot_df <- bind_rows(vario_emp_all_df, vario_th_tau)

  smooth_loess <- function(x, y, span = 0.3) {
    fit <- loess(y ~ x, span = span)
    tibble(hnormV = x, gamma = predict(fit, x))
  }

  vario_emp_smooth <- vario_plot_df %>%
    filter(type == "Empirical") %>%
    group_by(sim) %>%
    arrange(hnormV) %>%
    do(smooth_loess(.$hnormV, .$gamma, span = 0.4)) %>%
    ungroup()
  
  vario_th_smooth <- smooth_loess(vario_th_tau$hnormV, vario_th_tau$gamma,
                                  span = 0.4)

  p <- ggplot() +
  geom_line(data = subset(vario_plot_df, type == "Empirical"),
            aes(x = hnormV, y = gamma, group = sim),
            color = btfgreen, alpha = 0.2) +
  geom_line(data = vario_emp_smooth,
            aes(x = hnormV, y = gamma, group = sim),
            color = btfgreen, alpha = 0.5) +
  geom_line(data = subset(vario_plot_df, type == "Theoretical"),
            aes(x = hnormV, y = gamma),
            color = "#e05e5e", linewidth = 0.5, alpha=0.4) +
  geom_line(data = vario_th_smooth,
            aes(x = hnormV, y = gamma),
            color = "#e05e5e", linewidth = 1.2) +
  labs(x = "Spatial lag shifted by advection", y = "Variogram") +
  theme_minimal()

  # save the plot
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
  name_file <- paste0("/vario_optim_tau", tau, "_", M, "simu_", m, "rep_",
                      ngrid^2, "s_", length(temp), "t_", param_str, ".png")
  ggsave(plot = p, filename = paste0(foldername_image, name_file),
         width = 6, height = 4)
}




# Get data frame with the results for boxplot
# Convert the data frame to long format
df_bplot <- stack(df_result_all)
write.csv(df_bplot, paste0(foldername_result, name_file), row.names = FALSE)

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
  eta1 = 1,
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
                "s_", length(temp), "t_", param_str, "_beta_alpha.png")

ggsave(plot = bplot1, filename = paste0(foldername_image, name_file),
       width = 5, height = 5)

suffix <- if (all(!is.na(wind_df))) "etas" else "adv"

name_file <- paste0( "/bp_optim_", M, "simu_", m, "rep_", ngrid^2,
              "s_", length(temp), "t_", param_str, "_", suffix, ".png")

ggsave(plot = bplot2, filename = paste0(foldername_image, name_file),
       width = 5, height = 5)
