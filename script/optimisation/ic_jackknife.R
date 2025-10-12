foldername_results <- paste0(data_folder,
                    "/comephore/optim_results/lalpha/free_eta/")

filenames <- list.files(foldername_results, full.names = TRUE, pattern = "\\.csv$")
# keep only those that begin with "results_jk" not "all_results_jk"
file <- filenames[5]
# q95_delta30_dmin5.csv
result <- read.csv(file)
head(result)

foldername <- paste0(data_folder,
                    "/comephore/optim_results/jackknife_estimates/")

filenames <- list.files(foldername, full.names = TRUE, pattern = "\\.csv$")
# keep only those that begin with "results_jk" not "all_results_jk"
filenames <- filenames[grepl("all_results_jk", filenames)]

file <- filenames[1]

jack_estimates <- read.csv(file)
head(jack_estimates)

jack_estimates <- na.omit(jack_estimates)
n_eff <- nrow(jack_estimates)

# Save the jackknife estimates
filename <- paste0(data_folder, "comephore/optim_results/jackknife_estimates/all_results_jk_by_month_n", 
                   n_eff, "_q", q*100, "_delta", delta, "_dmin", min_spatial_dist, ".csv")
# write.csv(jack_estimates, filename, row.names = FALSE)

theta_hat <- as.numeric(result[1, ]) 
names(theta_hat) <- colnames(jack_estimates) %||% names(theta_hat)
p <- length(theta_hat)
n_eff <- nrow(jack_estimates)
pseudo_values <- matrix(NA, nrow = n_eff, ncol = p)
colnames(pseudo_values) <- names(theta_hat)

for (i in seq_len(n_eff)) {
  est_minus_i <- as.numeric(jack_estimates[i, ])
  # theta_i* = n * theta_hat - (n - 1) * theta_{(-i)}
  pseudo_values[i, ] <- n_eff * theta_hat - (n_eff - 1) * est_minus_i
}

pv_mean <- colMeans(pseudo_values)
pv_se   <- apply(pseudo_values, 2, sd) / sqrt(n_eff)

if (n_eff <= 30) {
  z <- qt(0.975, df = n_eff - 1)
} else {
  z <- qnorm(0.975)
}

ci_lower <- pv_mean - z * pv_se
ci_upper <- pv_mean + z * pv_se

results_jk <- data.frame(
  Parameter = names(theta_hat),
  Estimate  = theta_hat,
  PJ_mean   = pv_mean,
  StdError  = pv_se,
  CI_lower  = ci_lower,
  CI_upper  = ci_upper,
  row.names = NULL
)

print(results_jk)
