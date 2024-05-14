library(generain)

# spatial and temporal structures
ngrid <- 5
spa <- 1:ngrid
temp <- 1:30

# beta1, beta2, alpha1, alpha2
param <- c(0.6, 0.2, 1.5, 0.4) # true parameters for the variogram
beta1 <- param[1]
beta2 <- param[2]
alpha1 <- param[3]
alpha2 <- param[4]

n.BR <- 1
adv <- c(-0.5, -0.5)

BR <- sim_BR(param[1], param[2], param[3], param[4], spa, spa, temp, n.BR, adv)

# # save simulations to CSV files
# path <- paste0("../data/simulations_BR/sim_", ngrid, "s_", length(temp), "t/")
# save_simulations(BR, ngrid, n.BR, path)

save_simulations(BR, ngrid, n.BR,
                 folder = "./data/simulations_BR/",
                 file = "rain_br")

file_path <- paste0("./data/simulations_BR/rain_br_1.csv")
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
                       "/simu_gif/sim_br/br_25s_adv_4.gif"),
    ani.width = 700, ani.height = 600, ani.units = "px", ani.type = "cairo")
