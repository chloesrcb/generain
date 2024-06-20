
library(generain)
library(ggplot2)
library(reshape2)
library(animation)

# Simulate data using the rpareto model
ngrid <- 3
spa <- 1:ngrid
temp <- 1:50
n.res <- 2
param <- c(0.2, 0.2, 1.5, 0.8) # true parameters for the variogram
beta1 <- param[1] / 2
beta2 <- param[2] / 2
alpha1 <- param[3]
alpha2 <- param[4]
adv <- c(0.1, 0.1)

# Simulate spatio-temporal r-Pareto process
simu_rpar <- sim_rpareto(param[1], param[2], param[3], param[4], spa, spa, temp,
                          n.res, adv = adv)
plot(simu_rpar[1,,,], main = "BR simulation")
# Save the simulations
save_simulations(simu_rpar, ngrid, n.res,
                 folder = "../data/simulations_rpar/",
                 file = paste0("rpar_", ngrid^2, "s_",
                                length(temp), "t"))

file_path <- paste0("../data/simulations_rpar/rpar_", ngrid^2, "s_",
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
                        max(simulation_data_long$value))) +
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

# Save the plots as a gif
ani.options(interval = 0.5) # time between frames
saveGIF({
  for (i in temp) {
    print(plots[[i]])
  }
}, movie.name = paste0("/user/cserreco/home/Documents/These/generain/images",
                       "/simu_gif/simu_rpar/rpar_", ngrid^2, "s_",
                                length(temp), "t.gif"),
    ani.width = 700, ani.height = 600, ani.units = "px", ani.type = "cairo")



# validation 

# Simulation
list_rpar <- list()
for (i in 1:n.res) {
  file_path <- paste0("../data/simulations_rpar/rpar_", ngrid^2, "s_",
                                length(temp), "t_", i, ".csv")
  df <- read.csv(file_path)
  list_rpar[[i]] <- df
}

# dependece buhl
simu_df <- list_rpar[[1]] # first simulation
nsites <- ncol(simu_df) # number of sites

# get grid coordinates
sites_coords <- generate_grid_coords(sqrt(nsites))

# df_dist <- distances_regular_grid(nsites) # distance matrix
dist_mat <- get_dist_mat(sites_coords, adv = adv, tau = 1:10, latlon=FALSE) # distance matrix
df_dist <- reshape_distances(dist_mat) # reshape the distance matrix

# # Evaluate the estimates
# spa_estim_25 <- evaluate_vario_estimates(list_rpar, 0.9,
#                                   spatial = TRUE, df_dist = df_dist,
#                                   hmax = sqrt(17))

# temp_estim_25 <- evaluate_vario_estimates(list_rpar, 0.9,
#                                           spatial = FALSE, tmax = 10)

# df_result <- cbind(spa_estim_25, temp_estim_25)
# colnames(df_result) <- c("beta1", "alpha1", "beta2", "alpha2")

# df_valid <- get_criterion(df_result, true_param)

# optim

# df_dist <- distances_regular_grid(nsites) # distance matrix
