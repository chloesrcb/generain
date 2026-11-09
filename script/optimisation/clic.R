

theta_full <- c(
  0.1370892, # beta1
  1.1812283, # beta2
  0.6566290, # alpha1
  0.6639425 # alpha2
)

# compute the Hessian matrix at the full model estimates
H_full <- numDeriv::hessian(
  func = function(theta) {
    neg_ll_composite(params = theta, list_episodes = list_episodes,
      list_excesses = list_excesses, list_lags = list_lags, wind_df = NA,
      hmax = hmax,  latlon = FALSE, distance = "euclidean"
    )}, x = theta_full
)


m <- length(list_episodes)

ncl_one_episode <- function(params, i, list_episodes, list_excesses,
                            list_lags, wind_df, hmax) {

  neg_ll_composite(
    params = params,
    list_episodes = list(list_episodes[[i]]),
    list_excesses = list(list_excesses[[i]]),
    list_lags = list(list_lags[[i]]),
    wind_df = NA,
    hmax = hmax,
    latlon = FALSE,
    distance = "euclidean",
  )
}

scores_full <- parallel::mclapply(seq_len(m),
                  function(i) {
                      numDeriv::grad(
                        func = function(theta) {
                          ncl_one_episode( params = theta, i = i,
                            list_episodes = list_episodes, list_excesses = list_excesses,
                            list_lags = list_lags, wind_df = NA, hmax = hmax
                          )
                        }, x = theta_full) 
                        }, mc.cores = ncores
)

scores_full <- do.call(rbind, scores_full)

# compute J
J_full <- crossprod(scores_full)


p_eff_full <- sum(diag(solve(H_full, J_full)))

# CLIC = 2 * ncl + 2 * p_eff
ncl_full <- neg_ll_composite(params = theta_full, list_episodes = list_episodes,
                             list_excesses = list_excesses, list_lags = list_lags,
                             wind_df = NA, hmax = hmax, latlon = FALSE,
                             distance = "euclidean")
CLIC_full <- 2 * ncl_full + 2 * p_eff_full

p_eff_full
CLIC_full

# without adv
# 1264368


# > p_eff_full
# [1] 783.9725
# > CLIC_full
# [1] 1242562

# > p_eff_eta2_1
# [1] 611.6029
# CLIC_eta2_1
# > CLIC_eta2_1
# [1] 1244470

# > p_eff_full
# [1] 783.9725
# > CLIC_full
# [1] 1242562
