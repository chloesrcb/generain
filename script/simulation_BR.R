library(generain)
library(reshape2)
library(ggplot2)
# spatial and temporal structures
ngrid <- 5
spa <- 1:ngrid
nsites <- ngrid^2 # if the grid is squared
temp <- 1:300

# beta1, beta2, alpha1, alpha2
param <- c(0.8, 0.4, 1.5, 1) # true parameters for the variogram
beta1 <- param[1] / 2
beta2 <- param[2] / 2
# beta3 <- param[2] / 2
alpha1 <- param[3]
alpha2 <- param[4]
# alpha3 <- param[4]
n.BR <- 10
adv <- c(2, 1)
true_param <- c(beta1, beta2, alpha1, alpha2, adv)
BR <- sim_BR(param[1], param[2], param[3], param[4], spa, spa, temp, n.BR, adv)
BR <- sim_BR_3D(param[1], param[2], param[3], param[4], spa, spa, temp, n.BR, adv)

plot(BR[5,4,,1], main = "BR simulation")






# save simulations to CSV files
save_simulations(BR, ngrid, n.BR,
                 folder = "../data/simulations_BR/",
                 file = paste0("br_3D_adv_", ngrid^2, "s_",
                                length(temp), "t"))

file_path <- paste0("../data/simulations_BR/br_3D_adv_", ngrid^2, "s_",
                                length(temp), "t_1.csv")

simulation_data <- read.csv(file_path)
ngrid <- sqrt(ncol(simulation_data))  # Number of grid points in each dimension

simulation_data$Time <- rownames(simulation_data) # Add a time column
simulation_data_long <- melt(simulation_data) # Convert to long format
simulation_data_long$Time <- as.numeric(simulation_data_long$Time)
# Create a dataframe to represent grid points
grid <- expand.grid(x = 1:ngrid, y = 1:ngrid)

plots <- list()
# for each time step
for (i in unique(simulation_data_long$Time)) {
  # Add the simulated values to the grid dataframe
  grid$value <- simulation_data_long$value[simulation_data_long$Time == i]

  # Plot
  p <-  ggplot(data = grid, aes(x = x, y = y, fill = value)) +
    geom_tile() +
    scale_fill_gradient(low = "#70a7ae", high = "#9d503d",
                        name = "Rainfall in mm",
                        limits = c(min(simulation_data_long$value),
                        max(simulation_data_long$value))
                       ) +
    labs(title = paste0("t =", i, " | Betas: ", beta1, ", ", beta2,
                        " | Alphas: ",
                        alpha1, ", ", alpha2, " | Advection: ", adv[1],
                        ", ", adv[2])) +
    theme_minimal() +
    theme(plot.background = element_rect(fill = "#F9F8F6", color = "#F9F8F6"),
          panel.border = element_blank(),
          panel.grid = element_blank(),
          axis.text.x = element_blank(),
          axis.text.y = element_blank(),
          axis.title.x = element_blank(),
          axis.title.y = element_blank())

  plots[[i]] <- p
}

library(animation)
# Save the plots as a gif
ani.options(interval = 0.5) # time between frames
saveGIF({
  for (i in temp) {
    print(plots[[i]])
  }
}, movie.name = paste0("/user/cserreco/home/Documents/These/generain/images",
                       "/simu_gif/sim_br/br_3D_adv_", ngrid^2, "s_",
                                length(temp), ".gif"),
    ani.width = 700, ani.height = 600, ani.units = "px", ani.type = "cairo")


# Simulation
list_BR <- list()
for (i in 1:n.BR) {
  file_path <- paste0("../data/simulations_BR/br_3D_adv_", ngrid^2, "s_",
                                length(temp), "t_", i, ".csv")
  df <- read.csv(file_path)
  list_BR[[i]] <- df
}

# dependece buhl
simu_df <- BR
simu_df <- list_BR[[1]] # first simulation
nsites <- ncol(simu_df) # number of sites
# get grid coordinates
sites_coords <- generate_grid_coords(sqrt(nsites))
dist_mat <- get_dist_mat(sites_coords, adv = adv, tau = 1:10,
                         latlon = FALSE) # distance matrix
df_dist <- reshape_distances(dist_mat) # reshape the distance matrix


params <- c(0.4, 0.2, 1.5, 1, 0.1, 0.1)
nsites <- ncol(simu_df) # number of sites
# get grid coordinates
sites_coords <- generate_grid_coords(sqrt(nsites))
# sites_coords$id <- 1:nsites
h_vect <- get_lag_vectors(sites_coords, params,
                          hmax = sqrt(17), tau_vect = 1:10)

# Evaluate the estimates with Buhl method
spa_estim_25 <- evaluate_vario_estimates(list_BR, 0.8,
                                  spatial = TRUE, df_dist = df_dist,
                                  hmax = sqrt(17))

temp_estim_25 <- evaluate_vario_estimates(list_BR, 0.7,
                                          spatial = FALSE, tmax = 10)

df_result <- cbind(spa_estim_25, temp_estim_25)
colnames(df_result) <- c("beta1", "alpha1", "beta2", "alpha2")

df_valid <- get_criterion(df_result, true_param)


# get the number of simulations
n_res <- length(list_BR)
# create a dataframe to store the results
df_result <- data.frame(beta = rep(NA, n_res), alpha = rep(NA, n_res))
# for all simulations
for (n in 1:n_res) {
  simu_df <- as.data.frame(list_BR[[n]])
  if (spatial) {
    chi <- spatial_chi_alldist(df_dist, simu_df, quantile = 0.9,
                                hmax = hmax)
    print(n)
    params <- get_estimate_variospa(chi, weights = "none", summary = TRUE)
    print(n)
  } else {
    chi <- temporal_chi(simu_df, tmax = tmax, quantile = 0.9)
    params <- get_estimate_variotemp(chi, tmax, npoints = ncol(simu_df),
                                    weights = "exp", summary = FALSE)
  }
  df_result$beta[n] <- params[1]
  df_result$alpha[n] <- params[2]
}

