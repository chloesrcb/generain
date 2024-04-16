library(generain)

# spatial and temporal structures
ngrid <- 7
spa <- 1:ngrid
temp <- 1:100

# beta1, beta2, alpha1, alpha2
param <- c(0.8, 0.4, 1.5, 1) # true parameters for the 2*semivariogram
n.BR <- 100
BR <- sim_BR(param[1], param[2], param[3], param[4], spa, spa, temp, n.BR)

# save simulations to CSV files
path <- paste0("../data/simulations_BR/sim_", ngrid, "s_", length(temp), "t/")
save_simulations(BR, ngrid, n.BR, path)