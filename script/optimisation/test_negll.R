beta1_grid <- seq(0.1, 1, length.out = 10)
alpha1_grid <- seq(0.05, 1, length.out = 10)
vals <- matrix(NA, length(beta1_grid), length(alpha1_grid))
wind_df <- V_episodes_noadv
for (i in seq_along(beta1_grid)) {
  for (j in seq_along(alpha1_grid)) {
    pars <- c(beta1_grid[i], 4, alpha1_grid[j], 0.7, 1, 1)
    vals[i, j] <- neg_ll_composite(
      pars, list_episodes, list_excesses, list_lags, wind_df = wind_df
    )
  }
}
contour(beta1_grid, alpha1_grid, vals,
    xlab = expression(beta[1]),
    ylab = expression(alpha[1]),
    main = "Surface du Neg LL")
dev.copy(png, filename = "negll_surface_eta1_0.png")
dev.off()


start <- log(c(beta1 = 0.5, beta2 = 4, alpha1 = 0.1, alpha2 = 0.7, eta1 = 1, eta2 = 1)) # valeurs initiales plausibles
res <- optim(start, objfun, method = "BFGS")
best_pars <- exp(res$par)  # back-transform


neg_ll_composite(test_par, list_episodes, list_excesses, list_lags, wind_df)

for (eps in c(-0.1, -0.05, 0.05, 0.1)) {
  par2 <- test_par
  par2[1] <- par2[1] + eps
  cat("beta1 =", par2[1], " -> Neg LL =", 
      neg_ll_composite(par2, list_episodes, list_excesses, list_lags, wind_df), "\n")
}


objfun <- function(p) {
  neg_ll_composite_fixed_eta(p, list_episodes, list_excesses, list_lags, wind_df = wind_df,
  latlon = FALSE, distance = "lalpha", fixed_eta1 = init_params_com[5], fixed_eta2 = init_params_com[6])
}



result <- optim(
  par = init_params_com[1:4],
  fn = neg_ll_composite_fixed_eta,
  list_lags = list_lags,
  list_episodes = list_episodes,
  list_excesses = list_excesses,
  hmax = hmax,
  wind_df = V_episodes,
  latlon = FALSE,
  distance = "lalpha",
  fixed_eta1 = init_params_com[5],
  fixed_eta2 = init_params_com[6],
  method = "L-BFGS-B",
  lower = c(1e-08, 1e-08, 1e-08, 1e-08),
  upper = c(10, 10, 1.999, 1.999),
  control = list(maxit = 20000, trace = 1)
)

result

grad(objfun, result$par)

objfun(result$par)
objfun(result$par + c(1e-4, 0, 0, 0))
objfun(result$par - c(1e-4, 0, 0, 0))


# etas = c(1,1)
# > grad(objfun, result$par)
# [1]   25.22033   20.46536 -126.02358  -93.68068
# > result
# $par
# [1] 0.5743809 4.7635971 0.1575893 0.7107742
# $value
# [1] 19237.07
# > objfun(result$par)
# [1] 19182.59

# etas = init_params_com[5:6] = c(0.8, 3.8)
# > grad(objfun, result$par)
# [1]  0.040076029  0.003505651  0.036382788 -0.013105235
# > result
# $par
# [1] 0.5950533 4.6880691 0.1853842 0.7144501
# $value
# [1] 19180.16
# objfun(result$par)
# # [1] 19180.16


# eta 1 fixed at 0, init = init_params_com[1:4]
# $par
# [1] 0.50886955 4.87047484 0.02650187 0.70188455
# $value
# [1] 19168.69 (pas exactement la negll)
# objfun(result$par)
# # [1] 19248.59 # negll minimum

# > result
# $par
# [1] 0.5696485 4.6626892 0.1437467 0.7249424
# $value
# [1] 19172.37
# > grad(objfun, result$par)
# [1] -197.34782  -38.86865  -98.66691  183.47429
# > objfun(result$par)
# [1] 19186.13
# > objfun(result$par + c(1e-4, 0, 0, 0))
# [1] 19186.11
# > objfun(result$par - c(1e-4, 0, 0, 0))
# [1] 19186.15

# > result
# $par
# [1] 0.5950533 4.6880691 0.1853842 0.7144501

# $value
# [1] 19180.16

# $counts
# function gradient 
#       25       25 

# $convergence
# [1] 0

# $message
# [1] "CONVERGENCE: REL_REDUCTION_OF_F <= FACTR*EPSMCH"

# > grad(objfun, result$par)
# [1]  0.040076029  0.003505651  0.036382788 -0.013105235
# > objfun(result$par)
# [1] 19180.16
# > objfun(result$par + c(1e-4, 0, 0, 0))
# [1] 19180.16
# > objfun(result$par - c(1e-4, 0, 0, 0))
# [1] 19180.16

res_fix <- optim(
  par = init_param[1:4],
  fn = objfun_fixed_eta,
  method = "L-BFGS-B",
  lower = c(1e-8, 1e-8, 1e-8, 1e-8),
  upper = c(10, 10, 1.999, 1.999),
  control = list(maxit = 10000)
)



# verif par episodes
for (i in 1:m) {
  # extract episode and excesses from i-th r-pareto process from data
    excesses <- list_excesses[[i]]
    lags <- list_lags[[i]]
    adv <- as.vector(adv_df[i,])
    params_adv <- c(params[1:4], adv) # Add advection parameters

    nll_i <- neg_ll(params = params_adv,
                df_lags = lags,
                hmax = hmax, excesses = excesses,
                latlon = FALSE, distance = "lalpha",
                threshold = threshold, rpar = TRUE)
    cat("Episode", i, "Neg LL =", nll_i, "\n")
}

  nll_composite <- 0
  for (i in 1:m) {
    # extract episode and excesses from i-th r-pareto process from data
    excesses <- list_excesses[[i]]
    lags <- list_lags[[i]]

    if (!all(is.na(wind_df)) && length(adv_df) != 2) {
      adv <- as.vector(adv_df[i,])
    }

    params_adv <- c(params[1:4], adv) # Add advection parameters
    if (!rpar) {
      data <- list_episodes[[i]]
      quantile <- quantile
    }

    nll_i <- neg_ll(params = params_adv,
                    df_lags = lags,
                    hmax = hmax, excesses = excesses,
                    latlon = latlon, distance = distance,
                    threshold = threshold, rpar = rpar,
                    data = data, quantile = quantile)

    nll_composite <- nll_composite + nll_i
  }
