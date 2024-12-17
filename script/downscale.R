library(generain)
setwd("./script")
# load libraries
source("load_libraries.R")
library(MASS)
library(ggplot2)


generate_egpd <- function(n, kappa, scale, shape) {
  u <- runif(n)
  z <- dgpdExt1(u, kappa, scale, shape) # EGPD model 1
  return(z)
}

egpd_par_com <- c(0.6, 1, 0.5) # kappa, sigma, xi
egpd_par_ohsm <- c(0.5, 0.25, 0.5)


downscale_egpd <- function(x_grid, y_grid, kappa_grid, scale_grid, shape_grid,
                           x_points, y_points, kappa_points, scale_points, shape_points) {
  
  # Générer des valeurs EGPD pour la grille régulière
  z_simulated <- numeric(length(x_grid))
  for (i in 1:length(x_grid)) {
    z_simulated[i] <- generate_egpd(1, kappa_grid[i], scale_grid[i], shape_grid[i])
  }
  
  # Ajustement des valeurs simulées à l'échelle des points
  kappa_diff <- kappa_points - mean(kappa_grid)
  scale_ratio <- scale_points / mean(scale_grid)
  shape_diff <- shape_points - mean(shape_grid)
  
  z_downscaled <- (z_simulated + mean(kappa_diff) - kappa_diff) * (scale_ratio) + mean(shape_diff) - shape_diff
  
  return(z_downscaled)
}

# load data: regular grid
df_comephore <- read.csv("../data/comephore/inside_zoom.csv", sep = ",")
comephore <- df_comephore[-1] # remove dates column
# get location of each pixel
loc_px <- read.csv("../data/comephore/loc_pixels_zoom.csv", sep = ",")

# Get the longitude and latitude coordinates of the pixels
x_lon <- loc_px$Latitude
y_lat <- loc_px$Longitude

# Create a grid from longitude and latitude coordinates
grid <- expand.grid(x_lon = x_lon, y_lat = y_lat)

# EGPD param for the grid
kappa_grid <- runif(21, 0.5, 0.7)
scale_grid <- runif(21, 0.8, 1.2)
shape_grid <- runif(21, 0.3, 0.7)


# load data: non regular grid
load("../data/PluvioMontpellier_1min/rain_mtp_5min_2019_2022.RData")
ohsm <- rain.all5[c(1, 6:ncol(rain.all5))]
# get location of each rain gauge
location_gauges <- read.csv("../data/PluvioMontpellier_1min/pluvio_mtp_loc.csv")
location_gauges$codestation <- c("iem", "mse", "poly", "um", "cefe", "cnrs",
                                 "crbm", "archiw", "archie", "um35", "chu1",
                                 "chu2", "chu3", "chu4", "chu5", "chu6", "chu7")

x_points <- location_gauges$Latitude
y_points <- location_gauges$Longitude


kappa_points <- mean(kappa_grid)
scale_points <- mean(scale_grid)
shape_points <- mean(shape_grid)

# Downscaling
z_downscaled <- downscale_egpd(x_lon, y_lat, kappa_grid, scale_grid, shape_grid,
                               x_points, y_points, kappa_points, scale_points, shape_points)

# Add the downscaled values to the grid
grid$z <- z_downscaled

# Visualisation
ggplot(grid, aes(x = x, y = y, fill = z)) +
  geom_tile() +
  scale_fill_gradient(low = "white", high = "blue") +
  theme_minimal() +
  labs(title = "Downscaling de la pluie avec EGPD")


# revoir gamlss
# John Richards : 