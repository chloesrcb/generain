adv <- c(0.5, 0.2) # advection
params <- c(0.3, 0.7, 0.1, 0.7)
ngrid <- 5
temp <- 0:11
s0 <- c(1, 1)
t0 <- 0
random_s0 <- TRUE
M <- 50 # number of simulations
m <- 300 # number of extreme episodes
distance_type <- "euclidean" # "euclidean" or "lalpha"
s0_radius <- Inf
eta1 <- 4
eta2 <- 2
fixed_eta1 <- NA
fixed_eta2 <- NA
use_wind_data <- TRUE

