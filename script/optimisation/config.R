adv <- c(2, 1) # advection
params <- c(0.01, 0.2, 1.5, 1)
ngrid <- 5
temp <- 0:29
s0 <- c(1, 1)
t0 <- 0
random_s0 <- TRUE
M <- 7 # number of simulations
m <- 500 # number of extreme episodes
is_anisotropic <- FALSE
s0_radius <- Inf
fixed_eta1 <- TRUE
fixed_eta2 <- TRUE
use_wind_data <- TRUE
eta1 <- 2
eta2 <- 1