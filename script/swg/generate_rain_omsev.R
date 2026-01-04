rm(list = ls())
cat("\014")  # clear console

source("./script/load_libraries.R")

library(ggplot2)
library(reshape2)
library(animation)
library(RandomFields)
library(RandomFieldsUtils)
library(dplyr)
library(tidyr)
library(lubridate)
library(purrr)
library(viridis)
library(magick)
library(latex2exp)
library(sf)

functions_folder <- "./R"
files <- list.files(functions_folder, full.names = TRUE)
invisible(lapply(files, function(f) source(f, echo = FALSE)))


################################################################################
## Global parameters
################################################################################
q     <- 0.95
delta <- 12
dmin  <- 1200
Nsim  <- 1000
calm_threshold <- 0

################################################################################
## Load data
################################################################################

filename_rain <- paste0(
  data_folder,
  "omsev/omsev_5min/rain_mtp_5min_2019_2024_cleaned.csv"
)

filename_egpd <- paste0(
  data_folder,
  "../thesis/resources/images/EGPD/OMSEV/2019_2024/egpd_results.csv"
)

filename_loc <- paste0(data_folder, "omsev/loc_rain_gauges.csv")

filename_adv_episode <- paste0(
  data_folder,
  "omsev/adv_estim/combined_comephore_omsev/episode_advection_q",
  q * 100, "_delta", delta, "_dmin", dmin, ".csv"
)

filename_adv_classes <- paste0(
  data_folder,
  "omsev/adv_estim/combined_comephore_omsev/omsev_results_adv_classes.csv"
)

rain_omsev <- read.csv(filename_rain)
rownames(rain_omsev) <- rain_omsev$dates
rain <- rain_omsev[, -1]

egpd_params <- read.csv(filename_egpd)

location_gauges <- read.csv(filename_loc)
location_gauges$Station <- c(
  "iem","mse","poly","um","cefe","cnrs","crbm","archiw","archie","um35",
  "chu1","chu2","chu3","chu4","chu5","chu6","chu7","cines","brives","hydro"
)

stations_to_remove <- c("cines", "hydro", "brives")
rain <- rain[, !(colnames(rain) %in% stations_to_remove)]
location_gauges <- location_gauges[!(
  location_gauges$Station %in% stations_to_remove
), ]

dist_mat <- get_dist_mat(location_gauges)
df_dist  <- reshape_distances(dist_mat)

sites_coords <- location_gauges[, c("Longitude", "Latitude")]
rownames(sites_coords) <- location_gauges$Station

sites_coords_sf <- st_as_sf(
  sites_coords,
  coords = c("Longitude", "Latitude"),
  crs = 4326
)
sites_coords_sf <- st_transform(sites_coords_sf, crs = 2154)
coords_m <- st_coordinates(sites_coords_sf)

grid_omsev <- as.data.frame(coords_m / 1000)
colnames(grid_omsev) <- c("Longitude", "Latitude")
rownames(grid_omsev) <- rownames(sites_coords)

################################################################################
## Spatio-temporal extremes and episodes
################################################################################

set_st_excess <- get_spatiotemp_excess(
  rain,
  quantile     = q,
  remove_zeros = TRUE
)

s0t0_set <- get_s0t0_pairs(
  grid_omsev,
  rain,
  min_spatial_dist = dmin,
  episode_size     = delta,
  set_st_excess    = set_st_excess,
  n_max_episodes   = 10000,
  latlon           = FALSE
)

selected_points <- s0t0_set %>%
  mutate(
    t0_date = as.POSIXct(
      t0_date,
      format = "%Y-%m-%d %H:%M:%S",
      tz = "UTC"
    )
  )

t0_list <- selected_points$t0
s0_list <- selected_points$s0
u_list  <- selected_points$u_s0

episodes_list <- get_extreme_episodes(
  selected_points,
  rain,
  episode_size = delta,
  unif = FALSE
)

list_episodes <- episodes_list$episodes

################################################################################
## Marginal parameters
################################################################################

p0_values <- sapply(colnames(rain), function(s) {
  mean(rain[, s] == 0, na.rm = TRUE)
})

kappa_vect <- egpd_params[match(colnames(rain), egpd_params$Site), "kappa"]
xi_vect    <- egpd_params[match(colnames(rain), egpd_params$Site), "xi"]
sigma_vect <- egpd_params[match(colnames(rain), egpd_params$Site), "sigma"]

params_margins <- list(
  xi    = xi_vect,
  sigma = sigma_vect,
  kappa = kappa_vect,
  p0    = p0_values
)

times <- 0:(delta - 1)

################################################################################
## Advection classification
################################################################################

adv_df_raw <- read.csv(filename_adv_episode)

adv_df_transfo <- adv_df_raw
adv_df_transfo$vnorm <- sqrt(adv_df_raw$vx^2 + adv_df_raw$vy^2)
adv_df_transfo$vnorm_t <- etas_estimates[1] * adv_df_transfo$vnorm^etas_estimates[2]

adv_df_transfo$vx <- ifelse(
  adv_df_transfo$vnorm > 0,
  adv_df_transfo$vx / adv_df_transfo$vnorm * adv_df_transfo$vnorm_t,
  0
)

adv_df_transfo$vy <- ifelse(
  adv_df_transfo$vnorm > 0,
  adv_df_transfo$vy / adv_df_transfo$vnorm * adv_df_transfo$vnorm_t,
  0
)

plot(adv_df_raw$vx, adv_df_raw$vy, pch = 19,
     xlab = "vx", ylab = "vy")

plot(adv_df_transfo$vx, adv_df_transfo$vy, pch = 19,
     xlab = "vx", ylab = "vy", xlim = c(-10,10), ylim = c(-10,10))
speed <- sqrt(adv_df_transfo$vx^2 + adv_df_transfo$vy^2)

speed_class <- 
  as.character(cut(
    speed,
    breaks = c(0, 0.2, 1, 3, Inf),
    labels = c("still", "weak", "moderate", "strong"),
    include.lowest = TRUE
  ))
table(speed_class)

adv_class <- factor(
  speed_class,
  levels = c("still", "weak", "moderate", "strong")
)

selected_points$adv_group <- adv_class

################################################################################
## Initial intensity classes
################################################################################

X_s0_t0 <- sapply(seq_along(list_episodes), function(i) {
  ep  <- list_episodes[[i]]
  ep[1, s0_list[i]]
})

selected_points$X_s0_t0 <- X_s0_t0

q_init <- quantile(X_s0_t0, probs = c(0.25, 0.75))

selected_points$init_class <- cut(
  X_s0_t0,
  breaks = c(-Inf, q_init[1], q_init[2], Inf),
  labels = c("init_low", "init_mid", "init_high"),
  include.lowest = TRUE
)

################################################################################
## Utility functions
################################################################################

get_adv_matrix_for_group <- function(adv_df_raw, indices, etas_estimates) {
  adv_sub <- adv_df_raw[indices, ]

  adv_sub$vnorm <- sqrt(adv_sub$vx^2 + adv_sub$vy^2)
  adv_sub$vx_trans <- ifelse(
    adv_sub$vnorm > 0,
    adv_sub$vx / adv_sub$vnorm * adv_sub$vnorm, 0
  )
  adv_sub$vy_trans <- ifelse(
    adv_sub$vnorm > 0,
    adv_sub$vy / adv_sub$vnorm * adv_sub$vnorm, 0
  )

  # adv_sub$vx_final <- adv_sub$vx_trans * (1000 / 60) * 5
  # adv_sub$vy_final <- adv_sub$vy_trans * (1000 / 60) * 5

  as.matrix(adv_sub[, c("vx_final", "vy_final")])
}

summarise_simulations <- function(df_sims) {
  df_sims %>%
    group_by(time) %>%
    summarise(
      sim_mean   = mean(rainfall),
      sim_median = median(rainfall),
      sim_sd     = sd(rainfall),
      sim_low    = quantile(rainfall, 0.025),
      sim_high   = quantile(rainfall, 0.975),
      .groups    = "drop"
    )
}

summarise_real_group <- function(df_real_group) {
  df_real_group %>%
    group_by(time) %>%
    summarise(
      real_mean   = mean(rainfall),
      real_median = median(rainfall),
      real_low    = quantile(rainfall, 0.025),
      real_high   = quantile(rainfall, 0.975),
      .groups     = "drop"
    )
}

build_real_group_df <- function(list_episodes_group, s0_list_group) {
  do.call(rbind, lapply(seq_along(list_episodes_group), function(i) {
    ep  <- list_episodes_group[[i]]
    s0i <- s0_list_group[i]
    data.frame(
      time       = seq_len(nrow(ep)),
      rainfall   = ep[, s0i],
      episode_id = i
    )
  }))
}

simulate_group_episodes <- function(
  list_episodes_group, s0_list_group, u_list_group,
  adv_group_matrix, params_vario, params_margins,
  grid_omsev, times, Nsim
) {
  sims_group <- vector("list", Nsim)
  s0_sim     <- integer(Nsim)

  for (i in seq_len(Nsim)) {
    idx <- sample(seq_along(list_episodes_group), 1)

    sims_group[[i]] <- simulate_many_episodes(
      N              = 1,
      u              = 1000,
      u_emp          = u_list_group[idx],
      params_vario   = params_vario,
      params_margins = params_margins,
      coords         = grid_omsev,
      times          = times,
      adv            = adv_group_matrix,
      t0             = 0,
      s0             = s0_list_group[idx]
    )[[1]]
    s0_sim[i] <- s0_list_group[idx]
  }

  df_sims <- do.call(
    rbind,
    lapply(seq_len(Nsim), function(i) {
      data.frame(
        time     = seq_len(ncol(sims_group[[i]])),
        rainfall = sims_group[[i]][s0_sim[i], ],
        sim_id   = i
      )
    })
  )

  list(df_sims = df_sims)
}

simulate_group_from_mean <- function(
  list_episodes_group, s0_list_group,
  params_vario, params_margins, grid_omsev,
  times, adv_mean, Nsim
) {
  X_s0_t0 <- sapply(seq_along(s0_list_group), function(i) {
    list_episodes_group[[i]][1, s0_list_group[i]]
  })

  X_s0_t0_mean <- mean(X_s0_t0, na.rm = TRUE)

  sims <- vector("list", Nsim)

  for (i in seq_len(Nsim)) {
    s0i <- sample(s0_list_group, 1)

    sims[[i]] <- simulate_many_episodes(
      N              = 1,
      u              = 1000,
      u_emp          = X_s0_t0_mean,
      params_vario   = params_vario,
      params_margins = params_margins,
      coords         = grid_omsev,
      times          = times,
      adv            = adv_mean,
      t0             = 0,
      s0             = s0i
    )[[1]]
  }

  df_sims <- do.call(
    rbind,
    lapply(seq_len(Nsim), function(i) {
      s0i <- sample(s0_list_group, 1)
      data.frame(
        time     = seq_len(ncol(sims[[i]])),
        rainfall = sims[[i]][s0i, ],
        sim_id   = i
      )
    })
  )

  list(df_sims = df_sims)
}

################################################################################
## Load class parameters
################################################################################

adv_classes_params <- read.csv(filename_adv_classes)
# params <- c(1.766, 7.063, 0.420, 0.902, 0.868, 1.688)  # from omsev_jk.R
params <- c(1.5129632, 4.4833122, 0.2138693, 0.7012997, 0.6732099, 5.9844737)

convert_params <- function(beta1, beta2, alpha1, alpha2, c_x = 1, c_t = 1) {
  beta1_new <- beta1 / (c_x^alpha1)
  beta2_new <- beta2 / (c_t^alpha2)
  list(beta1 = beta1_new, beta2 = beta2_new)
}

# convert params from km/h to m/5min
c_x_m <- 1000    # for m
c_t_5min <- 12   # 1 hour = 12 * 5min

# convert params and ci to m/5min
params_m5min <- convert_params(result$par[1], result$par[2],
                               result$par[3], result$par[4],
                               c_x = c_x_m, c_t = c_t_5min)
params_vario <- list(
  beta1 = params_m5min$beta1,
  beta2 = params_m5min$beta2,
  alpha1 = result$par[3],
  alpha2 = result$par[4]
)

################################################################################
## Loop over all advection classes
################################################################################

params_vario <- list(
    beta1  = params[1],
    beta2  = params[2],
    alpha1 = params[3],
    alpha2 = params[4]
)

etas_estimates <- c(params[5], params[6])
group_adv <- "strong"
indices_group <- which(selected_points$adv_group == group_adv)


list_episodes_group <- list_episodes[indices_group]
s0_list_group       <- s0_list[indices_group]
u_list_group        <- u_list[indices_group]

adv_group_matrix <- get_adv_matrix_for_group(
    adv_df_raw, indices_group, etas_estimates
)

speed <- sqrt(adv_group_matrix[,1]^2 + adv_group_matrix[,2]^2)
print(summary(speed))


sim_res <- simulate_group_episodes(
    list_episodes_group,
    s0_list_group,
    u_list_group,
    adv_group_matrix,
    params_vario,
    params_margins,
    grid_omsev,
    times,
    Nsim
)

df_sims         <- sim_res$df_sims
df_sim_summary  <- summarise_simulations(df_sims)
df_real_group   <- build_real_group_df(list_episodes_group, s0_list_group)
df_real_summary <- summarise_real_group(df_real_group)

adv_mean <- colMeans(adv_group_matrix)

sim_mean_res <- simulate_group_from_mean(
  list_episodes_group,
  s0_list_group,
  params_vario,
  params_margins,
  grid_omsev,
  times,
  adv_mean,
  Nsim
)

df_sims_mean        <- sim_mean_res$df_sims
df_sim_summary_mean <- summarise_simulations(df_sims_mean)

# plot simulations vs real
p <- ggplot() +
  geom_ribbon(
    data = df_sim_summary_mean,
    aes(x = time, ymin = sim_low, ymax = sim_high),
    fill = "lightblue",
    alpha = 0.5
  ) +
  geom_line(
    data = df_sim_summary_mean,
    aes(x = time, y = sim_mean),
    color = "blue",
    size = 1
  ) +
  geom_line(
    data = df_real_summary,
    aes(x = time, y = real_mean),
    color = "red",
    size = 1
  ) +
  labs(
    title = paste0("Simulated vs Real Rainfall - Advection class: ", group_adv,
                   " (mean adv)"),
    x = "Time (5 min intervals)",
    y = "Rainfall"
  ) +
  theme_minimal()
p 
