
library(generain)
library(ggplot2)
library(reshape2)
library(animation)


################################################################################
# Simulate data using the rpareto model
ngrid <- 10
spa <- 1:ngrid
temp <- 1:1000
n.res <- 10
param <- c(0.4, 0.2, 1.5, 1) # true parameters for the variogram
beta1 <- param[1]
beta2 <- param[2]
alpha1 <- param[3]
alpha2 <- param[4]
adv <- c(0.5, 0.3)
s0 <- c(1, 1)
t0 <- 1

# Simulate spatio-temporal r-Pareto process
simu_rpar <- sim_rpareto(param[1], param[2], param[3], param[4], spa, spa, temp,
                          adv, s0, t0, n.res)

# Save the simulations
save_simulations(simu_rpar, ngrid, n.res,
                 folder = "../data/simulations_rpar/",
                 file = paste0("rpar_", ngrid^2, "s_",
                                length(temp), "t"))

file_path <- paste0("../data/simulations_rpar/rpar_", ngrid^2, "s_",
                                length(temp), "t_1.csv")
simulation_data <- read.csv(file_path)


create_simu_gif <- function(simulation_data, params, type = "rpar", 
                            forcedtemp = NA) {
  # type = "rpar" or "br"
  ngrid <- sqrt(ncol(simulation_data)) # Number of grid points in each dimension
  Tmax <- nrow(simulation_data) # Number of time steps
  if (is.na(forcedtemp)) {
    temp <- 1:Tmax
  } else {
    temp <- 1:forcedtemp
  }

  beta1 <- params[1]
  beta2 <- params[2]
  alpha1 <- params[3]
  alpha2 <- params[4]
  if (length(params) == 6) {
    adv <- params[5:6]
  } else {
    adv <- c(0, 0)
  }

  if (length(params) == 6) {
    adv <- params[5:6]
  } else {
    adv <- c(0, 0)
  }

  simulation_data$Time <- rownames(simulation_data) # Add a time column
  simulation_data_long <- melt(simulation_data) # Convert to long format
  simulation_data_long$Time <- as.numeric(simulation_data_long$Time)
  # Create a dataframe to represent grid points
  grid <- expand.grid(x = 1:ngrid, y = 1:ngrid)

  plots <- list()
  cropped_data <- simulation_data_long[simulation_data_long$Time %in% temp, ]
  # for each time step
  for (i in unique(cropped_data$Time)) {
    # Add the simulated values to the grid dataframe
    grid$value <- cropped_data$value[cropped_data$Time == i]

    # Plot
    p <-  ggplot(data = grid, aes(x = x, y = y, fill = value)) +
      geom_tile() +
      scale_fill_gradient(low = "#70a7ae", high = "#9d503d",
                          name = "Rainfall in mm",
                          limits = c(min(cropped_data$value),
                                     max(cropped_data$value))) +
      labs(title = paste0("t =", i, " | Betas: ", beta1, ", ", beta2,
                          " | Alphas: ",
                          alpha1, ", ", alpha2, " | Advection: ", adv[1],
                          ", ", adv[2])) +
      theme_minimal() +
      theme(plot.background = element_rect(fill = "#F9F8F6",
                                         color = "#F9F8F6"),
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
                         "/simu_gif/simu_", type, "/", type, "_", ngrid^2, "s_",
                         Tmax, "t.gif"),
  ani.width = 700, ani.height = 600, ani.units = "px", ani.type = "cairo")

}

create_simu_gif(simulation_data, c(param, adv), type = "rpar", forcedtemp = 30)

################################################################################
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
par(mfrow = c(1, 1))
plot(simu_df[, 1], main = "rpareto simulation")
# get grid coordinates
sites_coords <- generate_grid_coords(sqrt(nsites))

params <- c(param, adv)

df_lags <- get_conditional_lag_vectors(sites_coords, params, s0, t0,
                          hmax = sqrt(17), tau_vect = 0:10)

chi_theorical <- theorical_chi(params, df_lags)
chi <- unique(chi_theorical$chi)
plot(chi)
tau <- 0:10 # temporal lags
quantile <- 0.6

# get the empirical excesses
excesses <- empirical_excesses(simu_df, quantile, df_lags)
# plot(density(excesses$kij), main = "Excesses")


result <- optim(par = c(params), fn = neg_ll,
                  data = simu_df,
                  quantile = quantile,
                  df_lags = df_lags,
                  excesses = excesses,
                  locations = sites_coords,
                  hmax = sqrt(17),
                  s0 = s0,
                  t0 = t0,
                  method = "BFGS",
                  control = list(parscale = c(1, 1, 1, 1, 1, 1),
                                 maxit = 10000))

rmse_optim <- sqrt((result$par - params)^2)
print(rmse_optim)

################################################################################
# Verification
################################################################################

library(extRemes)
library(ismev)
library(POT)

rpar <- simulation_data$S1
threshold <- quantile(rpar, probs = 0.98)
rpar_exc <- rpar[rpar > threshold]
fit_gpd <- gpd.fit(rpar, threshold)
sigma <- fit_gpd$mle[1]
xi <- fit_gpd$mle[2]

theorical_qgpd <- qgpd(ppoints(rpar_exc), loc=min(rpar_exc),
                       shape=xi, scale=sigma)

qqplot(rpar_exc, theorical_qgpd, main = "GPD Q-Q plot",
  xlab = "Empirical quantiles",
  ylab = "Theoretical quantiles")

