
library(geoR)

grid = as.matrix(expand.grid(1:20, 1:20))
s0 = 110 # Index of conditioning location
distmat = as.matrix(dist(grid)) # Matrix of distances
# Define the power variogram
vario = function(h,beta=1,alpha=1.5){(h/beta)^alpha}
gamma = vario(distmat) # Variogram matrix
# Simulation of a Gaussian random field:
set.seed(1234)
G = geoR::grf(grid = grid, cov.pars = c(1, 1.5), cov.model="power")$data
# Compute the corresponding spectral process:
Y = exp(G-G[s0]-gamma[,s0])
# Simulation of a standard Pareto variable:
R = evd::rgpd(n = 1, loc = 1, scale = 1, shape = 1)
# Construct the simple Pareto process:
Z = R*Y

par(mfrow=c(1,1))
plot(Z)

grid_df <- as.data.frame(grid)
colnames(grid_df) <- c("x", "y")
grid_df$value <- Z

ggplot(grid_df, aes(x = x, y = y, fill = value)) +
  geom_tile() +
  scale_fill_gradient(low = "#80a5c6", high = "#be2035") +
  theme_minimal() +
  labs(title = "Downscaling de la pluie avec EGPD")


# spatial and temporal
# spatial and temporal structures
ngrid <- 2
x <- 1:ngrid
y <- 1:ngrid
z <- 1:2

alpha1 <- 1.5
beta1 <- 0.8
alpha2 <- 1
beta2 <- 0.4

RandomFields::RFoptions(spConform = FALSE, install="no")
lx <- length(sx <- seq_along(x))  # spatial
ly <- length(sy <- seq_along(y))  # spatial
lz <- length(sz <- seq_along(z))  # temporal
## Model-Variogram BuhlCklu
modelBuhlCklu <- RandomFields::RMfbm(alpha = alpha1, var = beta1, proj = 1) +
                 RandomFields::RMfbm(alpha = alpha1, var = beta1, proj = 2) +
                 RandomFields::RMfbm(alpha = alpha2, var = beta2, proj = 3)

## Construct spatio-temporal grid
Nxy <- lx * ly # spatial
N <- Nxy * lz # spatio-temporal
grid <- matrix(0, nrow=N, ncol=3) # (N,3)-matrix

for (i in sx) {
  for (j in seq_len(ly*lz))
    grid[i+(j-1)*ly, 1] <- i
}

for (i in sy) {
  for (j in sx)
    for(k in sz)
      grid[j+lx*(i-1)+(k-1)*Nxy, 2] <- i
}

for (i in sz) {
  for (j in seq_len(Nxy))
    grid[j+Nxy*(i-1), 3] <- i
}

adv <- c(0,0)
## Construct shifted variogram
Varm1 <- vapply(seq_len(N), function(n)
            RandomFields::RFvariogram(modelBuhlCklu,
                        x = (sx - grid[n, 1]) - adv[1] * grid[n, 3],
                        y = (sy - grid[n, 2]) - adv[2] * grid[n, 3],
                        z = sz - grid[n, 3]),
            array(NA_real_, dim = c(lx, ly, lz)))

sx0 = 1 # Index of conditioning location
sy0 = 1 # Index of conditioning location
t0 = 1 # Index of conditioning time
# Simulation of a Gaussian random field:
set.seed(1234)
W <- RandomFields::RFsimulate(modelBuhlCklu, x, y, z, n=1) 
# Compute the corresponding spectral process:
Y = exp(W - W[1] - Varm1[,,, 1])
# Simulation of a standard Pareto variable:
R = evd::rgpd(n = 1, loc = 1, scale = 1, shape = 1)
# Construct the simple Pareto process:
Z = R*Y


grid_df <- as.data.frame(grid)
colnames(grid_df) <- c("x", "y", "t")
grid_df$value <- Z

ggplot(grid_df, aes(x = x, y = y, fill = value)) +
  geom_tile() +
  scale_fill_gradient(low = "#80a5c6", high = "#be2035") +
  theme_minimal() +
  labs(title = "Downscaling de la pluie avec EGPD")


Z1 = rep(0, nrow(grid))
Z1 <- array(0, dim = c(lx, ly, lz))
while(max(Z1) < nrow(Z1)){
    W <- RandomFields::RFsimulate(modelBuhlCklu, x, y, z, n=1) 
    # Compute the corresponding spectral process:
    Y = exp(W - W[1] - Varm1[,,, 1])
    R = evd::rgpd(n = 1, loc = 1, scale = 1, shape = 1)
    Z1 = R*Y/mean(Y)
}
Z2 = Z1/nrow(Z1)

ngrid <- 9
spa <- 1:ngrid
temp <- 1:5
n.BR <- 1
simu_rpar <- sim_rpareto(0.8, 0.4, 1.5, 1, spa, spa, temp, n.BR, adv = c(0, 0))
save_simulations(simu_rpar, ngrid, n.BR,
                 folder = "../data/simulations_rpar/test/",
                 file = "rain_rpar")

library(animation)

file_path <- paste0("../data/simulations_rpar/test/rain_rpar_1.csv")
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
                        name = "Rainfall in mm") +
    labs(title = paste0("t =", i)) +
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
}, movie.name = "simulation_rpar.gif",
    ani.width = 700, ani.height = 600, ani.units = "px", ani.type = "cairo")

