adv <- c(0., 0.) # advection
params <- c(0.4, 0.2, 1.5, 1)
ngrid <- 5
temp <- 0:29
s0 <- c(1, 1)
t0 <- 0
random_s0 <- TRUE
M <- 2 # number of simulations
m <- 500 # number of extreme episodes
distance_type <- "lalpha" # "euclidean" or "lalpha"
s0_radius <- Inf
eta1 <- 1
eta2 <- 1
fixed_eta1 <- NA
fixed_eta2 <- NA
use_wind_data <- TRUE

wind_data <- adv
adv_real <- c(eta1 * abs(adv[1])^eta2 * sign(adv[1]),
              eta1 * abs(adv[2])^eta2 * sign(adv[2]))
true_param <- c(params, adv_real)
beta1 <- params[1]
beta2 <- params[2]
alpha1 <- params[3]
alpha2 <- params[4]
spa <- 1:ngrid

simu <- sim_rpareto(
    beta1 = beta1,
    beta2 = beta2,
    alpha1 = alpha1,
    alpha2 = alpha2,
    x = spa,
    y = spa,
    t = temp,
    adv = adv_real,
    t0 = t0,
    nres = m,
    random_s0 = random_s0,
    s0 = s0,
    s0_radius = s0_radius,
    distance = distance_type
)

list_rpar <- convert_simulations_to_list(simu$Z, ngrid)
head(list_rpar[[1]])
S5 <- list_rpar[[1]]$S5
head(S5)
library(ggplot2)

df <- data.frame(Index = 1:length(S5), S5 = S5)

ggplot(df, aes(x = Index, y = S5)) +
  geom_line(color = btfgreen, linewidth = 1) +
  theme_minimal() 



list_lags_excesses <- lapply(1:m, function(j) {
    s0_x <- simu$s0_used[[j]][[1]]$x
    s0_y <- simu$s0_used[[j]][[1]]$y
    s0_coords <- sites_coords[sites_coords$Longitude == s0_x &
                                sites_coords$Latitude == s0_y, ]
    lags <- get_conditional_lag_vectors(sites_coords, s0_coords, t0,
                                        tau_vect, latlon = FALSE)
    excesses <- empirical_excesses_rpar(list_rpar[[j]], quantile = u,
                                        df_lags = lags,
                                        t0 = t0, threshold = TRUE)
  list(lags = lags, excesses = excesses)
})


list_lags <- lapply(list_lags_excesses, `[[`, "lags")
list_excesses <- lapply(list_lags_excesses, `[[`, "excesses")


# Optimization
init_params <- c(params[1:4], 1, 1)
result_list <- mclapply(1:M, process_simulation, m = m,
                        list_simu = list_rpar, u = u,
                        list_lags = list_lags, t0 = t0,
                        list_excesses = list_excesses,
                        init_params = init_params,
                        hmax = 7, wind_df = wind_df,
                        mc.cores = num_cores)


result <- optim(
    par = params[1:4],
    fn = neg_ll_composite,
    list_episodes = list_rpar,
    list_lags = list_lags,
    list_excesses = list_excesses,
    hmax = 7,
    threshold = TRUE,
    latlon = FALSE,
    directional = TRUE,
    method = "L-BFGS-B",
    lower = c(1e-8, 1e-8, 1e-8, 1e-8, -Inf, -Inf),
    upper = c(Inf, Inf, 1.999, 1.999, Inf, Inf),
    control = list(maxit = 10000)
)
