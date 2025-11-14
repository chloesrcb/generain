
foldername <- file.path("./data/comephore/jackknife")
list_files <- list.files(foldername)
file <- list_files[2]
jack_estimates <- read.csv(file.path(foldername, file))
colnames(jack_estimates) <- c("beta1", "beta2", "alpha1", "alpha2", "eta1", "eta2")
n_eff <- nrow(jack_estimates)
jack_mean <- colMeans(jack_estimates)
jack_estimates <- as.matrix(jack_estimates)
theta_full <- round(c(0.308153284869901, 0.734418815078333, 
                      0.388850070714932,	0.694475390953746,	1.45355613609783,	
                      5.42017737398839,	270800.229103946), 4)[1:6]

pseudo_values <- matrix(NA, nrow = n_eff, ncol = 6)
for (i in 1:n_eff) {
  pseudo_values[i, ] <- n_eff * jack_mean - (n_eff - 1) * jack_estimates[i, ]
}

jack_mean_pseudo <- colMeans(pseudo_values)
jack_se <- apply(pseudo_values, 2, sd) / sqrt(n_eff)

z <- qnorm(0.975)
lower_ci <- jack_mean_pseudo - z * jack_se
upper_ci <- jack_mean_pseudo + z * jack_se

jackknife_seasonyear_results <- data.frame(
  Parameter = c("beta1", "beta2", "alpha1", "alpha2", "eta1", "eta2"),
  Estimate_full = theta_full,
  Estimate_jk   = jack_mean_pseudo,
  StdError = jack_se,
  CI_lower = lower_ci,
  CI_upper = upper_ci
)
