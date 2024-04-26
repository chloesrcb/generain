
library(generain)
library(ggplot2)
library(reshape2)
library(animation)

# Simulate data using the rpareto model
ngrid <- 7
spa <- 1:ngrid
temp <- 1:100
n.res <- 2
beta1 <- 0.2
beta2 <- 0.1
alpha1 <- 1.5
alpha2 <- 0.5
adv <- c(-1, -1)

simu_rpar <- sim_rpareto(beta1, beta2, alpha1, alpha2, spa, spa, temp,
                          n.res, adv = adv)

save_simulations(simu_rpar, ngrid, n.res,
                 folder = "./data/simulations_rpar/test/",
                 file = "rain_rpar")


file_path <- paste0("./data/simulations_rpar/test/rain_rpar_1.csv")
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
ani.options(interval = 0.4) # time between frames
saveGIF({
  for (i in 1:100) {
    print(plots[[i]])
  }
}, movie.name = paste0("/user/cserreco/home/Documents/These/generain/images",
                       "/simu_gif/variations_adv/rpar_6.gif"),
    ani.width = 700, ani.height = 600, ani.units = "px", ani.type = "cairo")
