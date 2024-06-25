
library(generain)
library(ggplot2)
library(reshape2)
library(animation)

# Simulate data using the rpareto model
ngrid <- 5
spa <- 1:ngrid
temp <- 1:300
n.res <- 2
param <- c(0.4, 0.2, 1.5, 0.8) # true parameters for the variogram
beta1 <- param[1] / 2
beta2 <- param[2] / 2
alpha1 <- param[3]
alpha2 <- param[4]
adv <- c(0.1, 0.1)

# Simulate spatio-temporal r-Pareto process
simu_rpar <- sim_rpareto(param[1], param[2], param[3], param[4], spa, spa, temp,
                          n.res, adv = adv)
plot(simu_rpar[5,5,,1], main = "BR simulation")
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

params <- c(param, adv)

h_vect <- get_lag_vectors(sites_coords, params,
                          hmax = sqrt(17), tau_vect = 1:10) # lag vectors

tau <- 1:10 # temporal lags
quantile <- 0.9
nmin <- 5

# get the empirical excesses
excesses <- empirical_excesses(simu_df, quantile, tau, h_vect,
                              nmin)

parscale <- c(1, 1, 1, 1, 1, 1)  # scale parameters
lower.bound <- c(1e-6, 1e-6, 1e-6, 1e-6, -1e-6, -1e-6)
upper.bound <- c(Inf, Inf, 1.999, 1.999, Inf, Inf)

result <- optimr(par = params, method = "Rcgmin",
                  gr = "grfwd", fn = function(par) {
                  neg_ll(par, simu = simu_df, quantile = quantile,
                        h_vect = h_vect, tau = tau,
                        locations = sites_coords)
                  }, lower = lower.bound, upper = upper.bound,
                  control = list(parscale = parscale,
                                 maxit = 10000))


result <- optimr(par = params, method = "CG", fn = function(par) {
                  neg_ll(par, simu = simu_df, quantile = quantile,
                        h_vect = h_vect, tau = tau,
                        locations = sites_coords)
                  },
                  control = list(parscale = parscale,
                                 maxit = 10000, maximise = FALSE))