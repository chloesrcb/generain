# Get all files in the folder "R"
functions_folder <- "./R"
files <- list.files(functions_folder, full.names = TRUE)
# # load all functions in files
invisible(lapply(files, function(f) source(f, echo = FALSE)))
library(latex2exp)
library(sf)
library(dplyr)
library(ggplot2)
library(gridExtra)


# First model
foldername <- paste0(data_folder,
                    "/comephore/optim_results/lalpha/free_eta/")

filename <- paste0(foldername, "combined_optim_results.csv")

res_mod1 <- read.csv(filename)

# remove rows with NA
res_mod1 <- res_mod1[complete.cases(res_mod1), ]

# Second model with omega
foldername <- paste0(data_folder,
                    "/comephore/optim_results/omega_etas/")

filename <- paste0(foldername, "combined_optim_results_omega.csv")

res_mod2 <- read.csv(filename)

# remove rows with NA
res_mod2 <- res_mod2[complete.cases(res_mod2), ]
# remove nll column
res_mod2 <- res_mod2[, !names(res_mod2) %in% c("nll")]

################################################################################
# Fictive advection vectors for representation
fictive_adv <- data.frame( # km/h, advection
  adv_x = c(  8, 12,   5,   2,  15, -5,  3),
  adv_y = c(  2,  6, -10,   1,   8, -3, -2)
)
# Estimated parameters from optimization
theta_hat <- c(res_mod1[2, "beta1"], res_mod1[2, "beta2"],
               res_mod1[2, "alpha1"], res_mod1[2, "alpha2"])
eta1_hat <- res_mod1[2, "eta1"]
eta2_hat <- res_mod1[2, "eta2"]

# Lag distances (km)
h_vals <- seq(0, 100, by = 0.5)

# Time lags (hours)
tau_vals <- c(0, 1, 3, 5, 7, 10)

# Standard spatial directions (unit vectors)
directions_named <- list(
  EW = c(1, 0),
  NS = c(0, 1),
  Diagonal = c(1, 1) / sqrt(2)
)

# Output folder for saving plots (update path as needed)
foldername <- paste0(im_folder, "optim/comephore/variogram/q", q * 100,
                     "_delta", delta, "_dmin", min_spatial_dist, "_fictive_adv/")
if (!dir.exists(foldername)) dir.create(foldername, recursive = TRUE)
foldername_hnorm <- paste0(foldername, "hnorm/")
if (!dir.exists(foldername_hnorm)) dir.create(foldername_hnorm, recursive = TRUE)

# Loop over all fictive advection vectors
for (i in seq_len(nrow(fictive_adv))) {
  fictive_v <- c(fictive_adv$adv_x[i], fictive_adv$adv_y[i])
  
  # Normalize advection direction vector
  norm_v <- sqrt(sum(fictive_v^2))
  if (norm_v > 0) {
    dir_adv <- fictive_v / norm_v
  } else {
    dir_adv <- c(1, 0)  # default if zero vector
  }
  
  # Add the advection direction to standard directions
  directions_all <- c(directions_named, list(Advection = dir_adv))
  
  # Loop over all directions
  for (dname in names(directions_all)) {
    direction <- directions_all[[dname]]
    
    # Compute gamma grid with corrected distance for each direction and advection vector
    df_gamma <- compute_gamma_grid(h_vals, tau_vals, direction,
                                   theta_hat,
                                   eta1 = eta1_hat, eta2 = eta2_hat,
                                   fictive_v = fictive_v)
    
    # Plot gamma as function of corrected distance (dist_corr)
    subtitle_txt <- paste0("Advection: V = (",
                           round(fictive_v[1], 3), ", ", round(fictive_v[2], 3), "), direction: ", dname)
    
    p <- ggplot(df_gamma, aes(x = h_norm, y = gamma, color = factor(tau))) +
        geom_line(size = 1.2) +
        labs(
          subtitle = subtitle_txt,
          x = expression("||h|| (km)"),
          y = expression(gamma(h, tau)),
          color = expression(tau ~ "(h)")
        ) +
        theme_minimal()

    filename <- paste0(foldername_hnorm, "variogram_fictive", i, "_dir_", dname, ".pdf")
    ggsave(filename, plot = p, width = 20, height = 15, units = "cm", dpi = 600)
  }
}


foldername_omega <- paste0(foldername, "omega/")
if (!dir.exists(foldername_omega)) dir.create(foldername_omega, recursive = TRUE)

fictive_adv <- data.frame( # km/h, advection
  adv_x = c(  8, 12,   5,   2,  15, -5,  3),
  adv_y = c(  2,  6, -10,   1,   8, -3, -2)
)

fictive_wind <- data.frame( # km/h, vent "classique"
  wx = c( 15,   5,   8,   3,  12, -4,  2),
  wy = c(  3,   4,  -7,   2,   9, -2, -1)
)

theta_hat <- c(res_mod2[2, "beta1"], res_mod2[2, "beta2"],
               res_mod2[2, "alpha1"], res_mod2[2, "alpha2"])
omega_hat <- res_mod2[2, "omega"]
eta1_hat <- res_mod2[2, "eta1"]
eta2_hat <- res_mod2[2, "eta2"]
# Loop over all fictive advection vectors
for (i in seq_len(nrow(fictive_adv))) {
  fictive_v <- c(fictive_adv$adv_x[i], fictive_adv$adv_y[i])
  fictive_w <- c(fictive_wind$wx[i], fictive_wind$wy[i])
  # Normalize advection direction vector
  norm_v <- sqrt(sum(fictive_v^2))
  if (norm_v > 0) {
    dir_adv <- fictive_v / norm_v
  } else {
    dir_adv <- c(1, 0)  # default if zero vector
  }
  
  # Add the advection direction to standard directions
  directions_all <- c(directions_named, list(Advection = dir_adv))
  
  # Loop over all directions
  for (dname in names(directions_all)) {
    direction <- directions_all[[dname]]
    
    # Compute gamma grid with corrected distance for each direction and advection vector
    df_gamma <- compute_gamma_grid_omega(h_vals, tau_vals, direction,
                                   theta_hat, omega = omega_hat,
                                   eta1 = eta1_hat, eta2 = eta2_hat,
                                   fictive_v = fictive_v, fictive_w = fictive_w)
    
    # Plot gamma as function of corrected distance (dist_corr)
    subtitle_txt <- paste0("Advection: V = (",
                           round(fictive_v[1], 3), ", ", round(fictive_v[2], 3), "), Wind: W = (",
                            round(fictive_w[1], 3), ", ", round(fictive_w[2], 3), "), direction: ", dname)
    
    p <- ggplot(df_gamma, aes(x = h_norm, y = gamma, color = factor(tau))) +
        geom_line(size = 1.2) +
        labs(
          subtitle = subtitle_txt,
          x = expression("||h|| (km)"),
          y = expression(gamma(h, tau)),
          color = expression(tau ~ "(h)")
        ) +
        theme_minimal()

    filename <- paste0(foldername_omega, "variogram_fictive", i, "_dir_", dname, ".pdf")
    ggsave(filename, plot = p, width = 20, height = 15, units = "cm", dpi = 600)
  }
}
