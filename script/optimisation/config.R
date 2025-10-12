adv <- c(0.5, 2) # advection
params <- c(0.4, 0.2, 1.5, 1)
ngrid <- 5
temp <- 0:29
s0 <- c(1, 1)
t0 <- 0
random_s0 <- TRUE
M <- 10 # number of simulations
m <- 500 # number of extreme episodes
distance_type <- "lalpha" # "euclidean" or "lalpha"
s0_radius <- Inf
eta1 <- 1
eta2 <- 1
fixed_eta1 <- NA
fixed_eta2 <- NA
use_wind_data <- TRUE