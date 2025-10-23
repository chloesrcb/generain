

convert_params <- function(beta1, beta2, alpha1, alpha2, c_x = 1, c_t = 1) {
  beta1_new <- beta1 / (c_x^alpha1)
  beta2_new <- beta2 / (c_t^alpha2)
  list(beta1 = beta1_new, beta2 = beta2_new)
}

beta1_hat <- result$par[1]
beta2_hat <- result$par[2]
alpha1_hat <- result$par[3]
alpha2_hat <- result$par[4]
# eta1_hat <- result$par[5]
# eta2_hat <- result$par[6]

# Conversion factors
c_x_km <- 1      # for km/h
c_t_h <- 1
c_x_m <- 1000    # for m/5min
c_t_5min <- 12   # 1 hour = 12 * 5min
# with etas fixed at (1, 1)
# res : [1] 0.51154558 5.59374130 0.09803423 0.80481171 1.00000000 1.00000000
# with etas fixed at (5, 2)
# 0.91978615 4.93468493 0.06751808 0.72966117 8.26000000 2.06000000
beta1_hat <- 0.9
beta2_hat <- 4.47
alpha1_hat <- 0.448
alpha2_hat <- 0.678
# Parameters in km/h
# params_kmh <- convert_params(beta1_hat, beta2_hat, alpha1_hat, alpha2_hat, c_x = c_x_km, c_t = c_t_h)
# Parameters in m/5min
params_m5min <- convert_params(beta1_hat, beta2_hat, alpha1_hat, alpha2_hat, c_x = c_x_m, c_t = c_t_5min)
# 2 m/5min 0.2598853 0.7571199 0.09803423 0.8048117    1    1
# 2 m/5min 0.2478567 0.7491045 0.1056372 0.8067743    1    1

# Table
param_table <- data.frame(
  Unit = c("km/h", "m/5min"),
  beta1 = c(params_kmh$beta1, params_m5min$beta1),
  beta2 = c(params_kmh$beta2, params_m5min$beta2),
  alpha1 = rep(alpha1_hat, 2),
  alpha2 = rep(alpha2_hat, 2),
  eta1 = rep(eta1_hat, 2),
  eta2 = rep(eta2_hat, 2)
)
print(param_table)

# save results to csv
results_filename <- paste0(data_folder, "omsev/optim_results/optim_results_q",
                           q * 100, "_delta", delta, "_dmin", min_spatial_dist, ".csv")
write.csv(param_table, results_filename, row.names = FALSE)