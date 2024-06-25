library(generain)
library(reshape2)
library(ggplot2)
setwd("./script")
# load libraries
source("load_libraries.R")

# spatial and temporal structures
ngrid <- 7
spa <- 1:ngrid
nsites <- ngrid^2 # if the grid is squared
temp <- 1:100

# beta1, beta2, alpha1, alpha2
param <- c(0.8, 0.4, 1.5, 1) # true parameters for the variogram
beta1 <- param[1] / 2
beta2 <- param[2] / 2
# beta3 <- param[2] / 2
alpha1 <- param[3]
alpha2 <- param[4]
# alpha3 <- param[4]

adv <- c(0, 0)
true_param <- c(beta1, beta2, alpha1, alpha2, adv)

library(foreach)
library(doParallel)

n.BR <- 1
num_iterations <- 100

# Create a parallel backend with the desired number of cores
cl <- makeCluster(detectCores())
registerDoParallel(cl)

# Run the simulation in parallel
results <- foreach(i = 1:num_iterations, .combine = "c") %dopar% {
  library(generain)
  BR <- sim_BR(param[1], param[2], param[3], param[4], spa, spa, temp, n.BR, 
               adv)
  save_simulations(BR, ngrid, n.BR,
                  folder = paste0("../data/simulations_BR/sim_", ngrid^2, "s_",
                                length(temp), "t/"),
                  file = paste0("br_", ngrid^2, "s_",
                                length(temp), "t"), forcedind = i)

}

# Stop the parallel backend
stopCluster(cl)


# plot(BR[5,4,,1], main = "BR simulation")


# save simulations to CSV files
save_simulations(BR, ngrid, n.BR,
                 folder = "../data/simulations_BR/",
                 file = paste0("br_", ngrid^2, "s_",
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
                       "/simu_gif/sim_br/br_", ngrid^2, "s_",
                                length(temp), ".gif"),
    ani.width = 700, ani.height = 600, ani.units = "px", ani.type = "cairo")


# Simulation
num_iterations <- 100
list_BR <- list()
for (i in 1:num_iterations) {
  file_path <- paste0("../data/simulations_BR/sim_25s_300t_adv/br_",
                      ngrid^2, "s_", length(temp), "t_", i, ".csv")
  df <- read.csv(file_path)
  list_BR[[i]] <- df
}

# dependece buhl
# simu_df <- BR
simu_df <- list_BR[[1]] # first simulation
nsites <- ncol(simu_df) # number of sites
plot(simu_df[, 1])
# get grid coordinates
sites_coords <- generate_grid_coords(sqrt(nsites))
dist_mat <- get_dist_mat(sites_coords,
                         latlon = FALSE) # distance matrix
df_dist <- reshape_distances(dist_mat) # reshape the distance matrix


params <- true_param
nsites <- ncol(simu_df) # number of sites
# get grid coordinates
sites_coords <- generate_grid_coords(sqrt(nsites))
# sites_coords$id <- 1:nsites
h_vect <- get_lag_vectors(sites_coords, params,
                          hmax = sqrt(17), tau_vect = 1:10)

# Evaluate the estimates with Buhl method
spa_estim <- evaluate_vario_estimates(list_BR, 0.6,
                                  spatial = TRUE, df_dist = df_dist,
                                  hmax = sqrt(17))

temp_estim <- evaluate_vario_estimates(list_BR, 0.9,
                                          spatial = FALSE, tmax = 10)

df_result <- cbind(spa_estim, temp_estim)
colnames(df_result) <- c("beta1", "alpha1", "beta2", "alpha2")

df_valid <- get_criterion(df_result, true_param)


# save the results in latex table
write.table(round(df_valid, 3), file = "valid_br.tex",
            row.names = c("Beta1", "Alpha1", "Beta2", "Alpha2"),
            col.names = c("Mean", "RMSE", "MAE"),
            quote = FALSE)
 