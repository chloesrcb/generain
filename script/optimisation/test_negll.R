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

wind_df = V_episodes
fixed_eta1 <- 0
objfun <- function(p) {
  neg_ll_composite_fixed_eta(p, list_episodes, list_excesses, list_lags, wind_df = wind_df,
  latlon = FALSE, distance = "lalpha", fixed_eta1 = fixed_eta1, fixed_eta2 = init_params_com[6])
}

# 12189.63 eta1 =etacom
#  eta1 =0
# > grad(objfun, result$par)
# [1] 0.10359255 0.01939044 0.05390137 0.03406592
# > grad(objfun, result$par)
# [1] -0.272812247  0.001803819  0.088287371 -0.027072175
# > objfun(result$par)
# [1] 12177.71
# > objfun(result$par + c(1e-4, 0, 0, 0))
# [1] 12177.71
# > objfun(result$par - c(1e-4, 0, 0, 0))
# [1] 12177.71

init_adv0 <- c(0.6788816, 4.6916841, 0.0115919, 0.6726300)
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
  fixed_eta1 = fixed_eta1,
  fixed_eta2 = init_params_com[6],
  method = "L-BFGS-B",
  lower = c(1e-08, 1e-08, 1e-08, 1e-08),
  upper = c(10, 10, 1.999, 1.999),
  control = list(maxit = 20000, trace = 1)
)

result

#$par
# [1] 0.7869492 4.5210071 0.1478401 0.6985827
# $value
# [1] 12204.73
# > grad(objfun, result$par)
# objfun(result$par) == result$value
# [1]  0.022664280  0.003596251  0.020450349 -0.015970467
# > objfun(result$par) == result$value
# [1] TRUE
# > 

result <- optim(
  par = c(init_adv0, 0, 1),
  fn = neg_ll_composite,
  list_lags = list_lags,
  list_episodes = list_episodes,
  list_excesses = list_excesses,
  hmax = hmax,
  wind_df = V_episodes,
  latlon = FALSE,
  distance = "lalpha",
  method = "L-BFGS-B",
  lower = c(1e-08, 1e-08, 1e-08, 1e-08),
  upper = c(10, 10, 1.999, 1.999),
  control = list(maxit = 20000, trace = 1)
)

result

# prendre comme init ce que j'obtiens avec eta1 = 0
# > grad(objfun, result$par)
# [1] -0.11902888 -0.03091853  0.17253076  0.17182187
# > result$par
# [1] 0.7959731 4.5932411 0.1535142 0.6982904

grad(objfun, result$par)

objfun(result$par) == result$value

help(optim)


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
