rm(list = ls())

# clear console
cat("\014")

# Load libraries and set theme
source("./script/load_libraries.R")

# library(generain)
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


# Get all files in the folder "R"
functions_folder <- "./R"
files <- list.files(functions_folder, full.names = TRUE)
# load all functions in files
invisible(lapply(files, function(f) source(f, echo = FALSE)))
library(latex2exp)

################################################################################
# FOCUS ON A SINGLE EPISODE
################################################################################


# get rain data
filename_rain <- paste0(data_folder, "omsev/omsev_5min/rain_mtp_5min_2019_2024_cleaned.csv")
rain_omsev <- read.csv(filename_rain)
head(rain_omsev)

# egpd fit
filename_egpd <- paste0(data_folder, "../thesis/resources/images/EGPD/OMSEV/2019_2024/egpd_results.csv")
egpd_params <- read.csv(filename_egpd)

# put dates as rownames
rownames(rain_omsev) <- rain_omsev$dates
rain <- rain_omsev[-1] # remove dates column

# get location of each rain gauge
filename_loc <- paste0(data_folder, "omsev/loc_rain_gauges.csv")
location_gauges <- read.csv(filename_loc)
location_gauges$Station <- c("iem", "mse", "poly", "um", "cefe", "cnrs",
                             "crbm", "archiw", "archie", "um35", "chu1",
                             "chu2", "chu3", "chu4", "chu5", "chu6", "chu7",
                             "cines", "brives", "hydro")


rain <- rain[, !(colnames(rain) %in% c("cines", "hydro", "brives"))]
location_gauges <- location_gauges[location_gauges$Station != "cines" &
                                   location_gauges$Station != "hydro" &
                                   location_gauges$Station != "brives", ]
dist_mat <- get_dist_mat(location_gauges)
df_dist <- reshape_distances(dist_mat)

sites_names <- colnames(rain)

sites_coords <- location_gauges[, c("Longitude", "Latitude")]

rownames(sites_coords) <- location_gauges$Station
sites_coords_sf <- st_as_sf(sites_coords, coords = c("Longitude", "Latitude"),
                            crs = 4326)
sites_coords_sf <- st_transform(sites_coords_sf, crs = 2154)
coords_m <- st_coordinates(sites_coords_sf)
grid_omsev <- as.data.frame(coords_m / 1000)
colnames(grid_omsev) <- c("Longitude", "Latitude")
rownames(grid_omsev) <- rownames(sites_coords)

grid_omsev_km <- grid_omsev

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

#######################################################################################################
# GET EPISODES
#######################################################################################################

q <- 0.95
delta <- 12
dmin <- 1200  # m
times <- 0:(delta - 1)

# in rain remove when all data are NA
set_st_excess <- get_spatiotemp_excess(rain, quantile = q, remove_zeros = TRUE)


# Spatio-temporal neighborhood parameters
s0t0_set <- get_s0t0_pairs(grid_omsev, rain,
                            min_spatial_dist = dmin,
                            episode_size = delta,
                            set_st_excess = set_st_excess,
                            n_max_episodes = 10000,
                            latlon = FALSE)

selected_points <- s0t0_set
selected_points <- selected_points %>%
  mutate(t0_date = as.POSIXct(t0_date, format = "%Y-%m-%d %H:%M:%S", tz = "UTC"))

n_episodes <- length(selected_points$s0)
t0_list <- selected_points$t0
s0_list <- selected_points$s0
u_list <- selected_points$u_s0
list_episodes_points <- get_extreme_episodes(selected_points, rain,
                                     episode_size = delta, unif = FALSE)
list_episodes <- list_episodes_points$episodes
tau_vect <- 0:10
tmax <- max(tau_vect)
df_coords <- as.data.frame(grid_omsev)
# Compute the lags and excesses for each conditional point
list_results <- mclapply(1:length(s0_list), function(i) {
  s0 <- s0_list[i]
  row_s0 <- which(rownames(df_coords) == s0)
  s0_coords <- df_coords[row_s0, ]
  episode <- list_episodes[[i]]
  ind_t0_ep <- 0 # index of t0 in the episode
  lags <- get_conditional_lag_vectors(df_coords, s0_coords, ind_t0_ep,
                                tau_vect, latlon = FALSE)
  u <- u_list[i]
  excesses <- empirical_excesses_rpar(episode, threshold = u,
                                      df_lags = lags, t0 = ind_t0_ep)
  # tau is in 5 minutes
  lags$tau <- lags$tau * 5 / 60 # convert to hours
  # lags$hnorm <- lags$hnorm / 1000 # convert to km
  list(lags = lags, excesses = excesses)
}, mc.cores = detectCores() - 1)

list_lags <- lapply(list_results, `[[`, "lags")
list_excesses <- lapply(list_results, `[[`, "excesses")
df_lags <- list_lags[[10]] # km/h
df_excesses <- list_excesses[[13]]
sum(df_excesses$kij)


#######################################################################################################
# VARIOS PARAMETERS FROM KM/H TO M/5MIN
#######################################################################################################
params_est <- c(0.4, 0.2, 1.5, 1, 1, 6)  # km/h
etas_estimates <- params_est[5:6]

params_est <- c(1.507, 4.466 , 0.206, 0.703, 1.092, 5.666)  # km/h
etas_estimates <- params_est[5:6]

params_vario_kmh <- list(
  beta1 = params_est[1],
  beta2 = params_est[2],
  alpha1 = params_est[3],
  alpha2 = params_est[4]
)
params_kmh <- c(
  params_vario_kmh$beta1,
  params_vario_kmh$beta2,
  params_vario_kmh$alpha1,
  params_vario_kmh$alpha2
)



##################################################################################
# TRANSFORM ADVECTION SPEEDS
##################################################################################
adv_filename <- paste(data_folder, "/omsev/adv_estim/combined_comephore_omsev/episode_advection_q",
                          q * 100, "_delta", delta, "_dmin", dmin,
                          ".csv", sep = "")
adv_df <- read.csv(adv_filename, sep = ",")
head(adv_df)
# count number of advection speed between 0 and 0.1 km/h and with mean_dy_comephore == 0
n_between_0_and_01 <- sum((sqrt(adv_df$vx_final^2 + adv_df$vy_final^2) > 0) &
                           (sqrt(adv_df$vx_final^2 + adv_df$vy_final^2) < 0.1))
cat("Number of episodes with advection between 0 and 0.1 km/h:", n_between_0_and_01, "\n")

index_01_com0 <- which((sqrt(adv_df$vx_final^2 + adv_df$vy_final^2) > 0) &
                      (sqrt(adv_df$vx_final^2 + adv_df$vy_final^2) < 0.1) &
                      (adv_df$mean_dy_kmh_comephore == 0))
adv_df[index_01_com0, c("vx_final", "vy_final")] <- 0

adv_df_transfo <- adv_df
adv_df_transfo$vnorm <- sqrt(adv_df$vx_final^2 + adv_df$vy_final^2)
adv_df_transfo$vnorm_t <- etas_estimates[1] * adv_df_transfo$vnorm^etas_estimates[2]

adv_df_transfo$vx_t <- ifelse(
  adv_df_transfo$vnorm > 0,
  adv_df_transfo$vx_final / adv_df_transfo$vnorm * adv_df_transfo$vnorm_t,
  0
)

adv_df_transfo$vy_t <- ifelse(
  adv_df_transfo$vnorm > 0,
  adv_df_transfo$vy_final / adv_df_transfo$vnorm * adv_df_transfo$vnorm_t,
  0
)

plot(adv_df$vx_final, adv_df$vy_final, pch = 19,
     xlab = "vx", ylab = "vy")

plot(adv_df_transfo$vx_t, adv_df_transfo$vy_t, pch = 19,
     xlab = "vx", ylab = "vy")

# Compute speed
speed_raw <- sqrt(adv_df$vx_final^2 + adv_df$vy_final^2)
summary(speed_raw)
speed_transfo <- sqrt(adv_df_transfo$vx_t^2 + adv_df_transfo$vy_t^2)
summary(speed_transfo)

direction <- atan2(adv_df$vy_final, adv_df$vx_final) * (180 / pi)
# make sure direction is in [0, 360]
angle_deg <- (direction + 360) %% 360
is_calm <- speed_raw == 0
speed_class <-as.character(cut(
    speed_raw,
    breaks = c(0, 0.5, 1, Inf),
    labels = c("still","weak","significant"),
    include.lowest = TRUE
  ))

speed_class <- factor(speed_class, levels = c("still","weak","significant"))
direction_class <- ifelse(
  is_calm,
  "none",
  as.character(cut(
    angle_deg,
    breaks = c(0, 90, 180, 270, 360),
    labels = c("E", "S", "W", "N"),
    include.lowest = TRUE
  ))
)
direction_class <- factor(direction_class, levels = c("none","E","S","W","N"))

adv_class <- data.frame(
  speed = speed_raw,
  speed_transfo = speed_transfo,
  direction = direction,
  direction_class = direction_class,
  group = paste0(speed_class, "_", direction_class)
)

table(adv_class$group)
adv_df <- adv_df_transfo # km/h

selected_points$adv_group <- adv_class$group


X_s0_t0 <- numeric(length(list_episodes))
for (i in seq_along(list_episodes)) {
  ep <- list_episodes[[i]]
  s0 <- s0_list[i]
  X_s0_t0[i] <- ep[1, s0]
}

selected_points$X_s0_t0 <- X_s0_t0

q_init <- quantile(selected_points$X_s0_t0, probs = c(0.25, 0.75))


selected_points$init_class <- cut(
  selected_points$X_s0_t0,
  breaks = c(-Inf, q_init[1], q_init[2], Inf),
  labels = c("init_low", "init_mid", "init_high"),
  include.lowest = TRUE
)

head(selected_points)




################################################################################
# Simulation for similar episodes advection
################################################################################
adv_matrix <- as.matrix(adv_df_transfo[, c("vx_t", "vy_t")])
all_group_names <- unique(selected_points$adv_group)
group_adv <- all_group_names[12]  # choose one group to simulate
#get list_episodes and s0_list for this group
indices_group <- which(selected_points$adv_group == group_adv)
list_episodes_group <- list_episodes[indices_group]
length(list_episodes_group)
s0_list_group <- s0_list[indices_group]
u_list_group <- u_list[indices_group]
adv_group <- adv_matrix[indices_group, ]
Nsim <- 1000
s0_sim <- integer(Nsim)
sims_group <- vector("list", Nsim)
adv_sim <- matrix(0, nrow = Nsim, ncol = 2)
for (i in seq_len(Nsim)) {
  idx <- sample(seq_along(list_episodes_group), 1)
  s0_i <- s0_list_group[idx]
  s0_sim[i] <- s0_i
  u_i <- u_list_group[idx]
  adv_i <- adv_group[idx, ]
  adv_sim[i, ] <- adv_i
  sims_group[[i]] <- simulate_many_episodes(
    N = 1,
    u = 1000,
    u_emp = u_i,
    params_vario = params_vario_kmh,
    params_margins = params_margins,
    coords = grid_omsev,
    times = times * 5 / 60,  # in hours
    adv = adv_i,
    t0 = 0,
    s0 = s0_i
  )[[1]]
}

# check that we have exceedances at s0, t0 for each episode
u_sim <- sapply(sims_group, function(sim) sim$u_latent)
all(sapply(seq_len(Nsim), function(i) {
  sims_group[[i]]$Z[s0_sim[i], 1] >= u_sim[i]
}))


df_sims <- do.call(rbind, lapply(seq_len(Nsim), function(i) {
  s0_i <- s0_sim[i]
  data.frame(
    time = seq_len(ncol(sims_group[[i]]$X)),
    rainfall = sims_group[[i]]$X[s0_i, ],
    sim_id = i
  )
}))


df_sim_summary <- df_sims %>%
  group_by(time) %>%
  summarise(
    sim_mean = mean(rainfall),
    sim_median = median(rainfall),
    sim_sd  = sd(rainfall),
    sim_low  = quantile(rainfall, 0.025),
    sim_high = quantile(rainfall, 0.975)
  )


df_real_group <- do.call(rbind, lapply(seq_along(list_episodes_group), function(idx) {
  ep <- list_episodes_group[[idx]]
  s0_i <- s0_list_group[idx]
  data.frame(
    time = seq_len(nrow(ep)),
    rainfall = ep[, s0_i],
    episode_id = idx
  )
}))


df_real_summary <- df_real_group %>%
  group_by(time) %>%
  summarise(
    real_mean = mean(rainfall),
    real_median = median(rainfall),
    real_low  = quantile(rainfall, 0.025),
    real_high = quantile(rainfall, 0.975)
  )

ggplot() +
  geom_ribbon(
    data = df_sim_summary,
    aes(x = time, ymin = sim_low, ymax = sim_high, fill = "Simulations"),
    alpha = 0.25
  ) +
  geom_line(
    data = df_sim_summary,
    aes(x = time, y = sim_median, color = "Simulations"),
    linewidth = 1
  ) +
  geom_line(
    data = df_real_summary,
    aes(x = time, y = real_median, color = "Observed episodes"),
    linewidth = 1
  ) +  geom_ribbon(
    data = df_real_summary,
    aes(x = time, ymin = real_low, ymax = real_high, fill = "Observed episodes"),
    alpha = 0.25
  ) +
  ylim(0, 11)+
  scale_color_manual(values = c("Simulations"="#6b6b6b","Observed episodes"="#c73535")) +
  scale_fill_manual(values = c("Simulations"="#bdbdbd","Observed episodes"="#e28b8b")) +
  labs(x="Time step", y="Rainfall (mm/5min)", color="", fill="") +
  theme_minimal() +
  btf_theme


ggplot() +
  # --- Simulations: IC 95 %
  geom_ribbon(
    data = df_sim_summary,
    aes(x = time, ymin = sim_low, ymax = sim_high, fill = "Simulations"),
    alpha = 0.30
  ) +

  # --- Simulations: médiane
  geom_line(
    data = df_sim_summary,
    aes(x = time, y = sim_median, color = "Simulations"),
    linewidth = 1.1
  ) +

  # --- Épisodes observés: séries individuelles (fin, transparent)
  geom_line(
    data = df_real_group,
    aes(x = time, y = rainfall, group = episode_id, color = "Observed episodes"),
    linewidth = 0.5,
    alpha = 0.25
  ) +

  ylim(0, 11) +

  scale_color_manual(
    values = c(
      "Simulations" = "#5c5c5c",
      "Observed episodes" = "#c73535"
    )
  ) +
  scale_fill_manual(
    values = c(
      "Simulations" = "#bdbdbd",
      "Observed episodes" = "#e28b8b"
    )
  ) +

  labs(
    x = "Time step",
    y = "Rainfall (mm / 5 min)",
    color = "",
    fill  = ""
  ) +

  theme_minimal() +
  btf_theme


# save plot

foldername_plot <- paste0(
  im_folder,
  "swg/omsev/all/"
)
if (!dir.exists(foldername_plot)) {
  dir.create(foldername_plot, recursive = TRUE)
}
filename <- paste0(
  foldername_plot,
  "simulated_vs_observed_rainfall_advgroup_",
  group_adv,
  "_n", Nsim,
  ".png"
)

ggsave(
  filename = filename,
  plot = last_plot(),
  width = 7,
  height = 5,
  dpi = 300
)

# plot just simulated rainfall at s0
ggplot(df_sims, aes(x = time, y = rainfall, group = sim_id)) +
  geom_line(alpha = 0.7, color = "gray") +
  labs(x="Time step (5min)", y="Rainfall (mm/5min)") +
  theme_minimal() +
  btf_theme

# save plot
filename_simonly <- paste0(
  foldername_plot,
  "simulated_rainfall_only_advgroup_",
  group_adv,
  "_n", Nsim,
  ".png"
)
ggsave(
  filename = filename_simonly,
  plot = last_plot(),
  width = 8,
  height = 6,
  dpi = 300
)

dist_mat <- get_dist_mat(location_gauges) / 1000
df_dist <- reshape_distances(dist_mat)
n_hbins <- 10
h_all <- df_dist$value

h_breaks <- quantile(
  h_all,
  probs = seq(0, 1, length.out = n_hbins + 1),
  na.rm = TRUE
)

h_breaks <- unique(as.numeric(h_breaks))
h_breaks[length(h_breaks)] <- 1.6


# on simulated episodes
list_results_sim <- lapply(seq_along(sims_group), function(i) {

  Z_ep <- sims_group[[i]]$Z
  s0 <- s0_sim[i]
  row_s0 <- which(rownames(df_coords) == s0)
  s0_coords <- df_coords[row_s0, ]
  ind_t0_ep <- 0

  u_0 <- sims_group[[i]]$u_latent

  lags <- get_conditional_lag_vectors(
    df_coords, s0_coords, ind_t0_ep,
    tau_vect, latlon = FALSE
  )

  excesses <- empirical_excesses_rpar(
    data_rain   = t(Z_ep),
    threshold = u_0,
    df_lags   = lags,
    t0        = ind_t0_ep
  )

  lags$tau <- lags$tau * 5/60 # convert to hours
  if(sum(excesses$kij) == 0) {
    cat("No exceedances for simulated episode", i, "\n")
  }
  list(lags = lags, excesses = excesses)
})



list_lags_sim <- lapply(list_results_sim, `[[`, "lags")
list_excesses_sim <- lapply(list_results_sim, `[[`, "excesses")
lags_i <- list_lags_sim[[10]]
summary(lags_i$hnorm)
summary(lags_i$tau)
adv_i <- adv_sim[10, ]
adv_x <- adv_i[1]
adv_y <- adv_i[2]
sqrt(adv_x^2 + adv_y^2) * max(lags_i$tau)


list_sumkij <- sapply(list_excesses_sim, function(df) sum(df$kij))
summary(list_sumkij)
table(list_sumkij > 0)


adv_df_sim <- as.data.frame(adv_sim)
colnames(adv_df_sim) <- c("vx", "vy")

params_estimates <- as.numeric(c(params_vario_kmh, etas_estimates))
p_sim <- plot_th_emp_chi(
  list_lags_sim,
  list_excesses_sim,
  wind_df_filtered = adv_df_sim,
  params_estimates = params_estimates,
  tau_min = 0,
  h_breaks = h_breaks,
  latlon = FALSE,
  adv_transform = FALSE
)


# save plot
foldername_plot <- paste0(
  im_folder,
  "swg/omsev/simu/"
)

if (!dir.exists(foldername_plot)) {
  dir.create(foldername_plot, recursive = TRUE)
}

filename_plot <- paste0(
  foldername_plot,
  "chi_vs_h_bytau_simulated_advgroup_",
  group_adv,
  "_n", Nsim,
  ".png"
)


# plot by tau
res_cmp_sim <- p_sim$res_cmp
res <- p_sim$res
# chi theoretical vs chi empirical with y=x line
ggplot(res_cmp_sim, aes(x = chi_theo_bar, y = chi_emp_bar)) +
  geom_point(alpha = 0.7, color=btfgreen) +
  geom_abline(slope = 1, intercept = 0, color = "red", linetype = "dashed") +
  labs(
    x = "Chi theoretical",
    y = "Chi empirical"
  ) +
  theme_minimal()

ggsave(
  filename = filename_plot,
  width = 7,
  height = 7,
  dpi = 300
)


res_boot <- res %>%
  dplyr::filter(.data$tau >= 0) %>%
  dplyr::group_by(.data$tau, .data$hbin) %>%
  dplyr::summarise(
    chi_emp_bar  = mean(.data$chi_emp, na.rm = TRUE),
    chi_theo_bar = mean(.data$chi_theo, na.rm = TRUE),
    ci_low  = bootstrap_ci(.data$chi_emp)[1],
    ci_high = bootstrap_ci(.data$chi_emp)[2],
    n_pairs = dplyr::n(),
    .groups = "drop"
  )


plot_all <- ggplot(res_boot, aes(
  x = chi_theo_bar,
  y = chi_emp_bar
)) +
  geom_errorbar(
    aes(ymin = ci_low, ymax = ci_high),
    width = 0,
    alpha = 0.4,
    color = btfgreen
  ) +
  geom_point(alpha = 0.7, color = btfgreen) +
  geom_abline(slope = 1, intercept = 0,
              color = "red", linetype = "dashed") +
  labs(
    title = "",
    x = "Theoretical Chi",
    y = "Empirical Chi"
  ) +
  theme_minimal() +
  btf_theme

plot_all

filename_plot <- paste0(
  foldername_plot,
  "chi_vs_h_bytau_simulated_advgroup_CI_",
  group_adv,
  "_n", Nsim,
  ".png"
)

ggsave(
  filename = filename_plot,
  plot = plot_all,
  width = 7,
  height = 7,
  dpi = 300
)

################################################################################

# look at marginals
sims_all_X <- do.call(cbind, lapply(sims_group, function(sim) sim$X))

# crps
library(scoringRules)

# sims_group : list of simulations
# each sim$X : matrix (n_sites x lt) with rownames = sites

sites <- rownames(sims_group[[1]]$X)
lt <- ncol(sims_group[[1]]$X)
Nsim <- length(sims_group)

# Stack all values for each site across sims and times
# Result: matrix n_sites x (Nsim*lt)
X_by_site <- do.call(cbind, lapply(sims_group, function(sim) sim$X))
# X_by_site rows = sites, cols = concatenated times across sims
dim(X_by_site)


Hgpd <- function(x, sigma, xi) {
  x <- pmax(x, 0)
  if (abs(xi) < 1e-10) return(1 - exp(-x / sigma))
  1 - (1 + xi * x / sigma)^(-1/xi)
}

Fegpd_pos <- function(x, kappa, sigma, xi) {
  stopifnot(all(x > 0))
  Hgpd(x, sigma, xi)^kappa
}

Qgpd <- function(p, sigma, xi) {
  p <- pmin(pmax(p, 1e-12), 1 - 1e-12)
  if (abs(xi) < 1e-10) return(-sigma * log(1 - p))
  sigma/xi * ((1 - p)^(-xi) - 1)
}

r_egpd_pos <- function(n, kappa, sigma, xi) {
  u <- runif(n)
  Qgpd(u^(1 / kappa), sigma, xi)
}

crps_from_samples <- function(y, xs) {
  mean(abs(xs - y)) - 0.5 * mean(abs(outer(xs, xs, "-")))
}

crps_egpd_mc_pos <- function(y, kappa, sigma, xi, m = 500) {
  if (y <= 0) return(NA_real_)
  xs <- r_egpd_pos(m, kappa, sigma, xi)
  crps_from_samples(y, xs)
}


X_by_site <- do.call(cbind, lapply(sims_group, function(sim) sim$X))
sites <- rownames(X_by_site)


crps_site <- sapply(seq_along(sites), function(i) {

  x <- as.numeric(X_by_site[i, ])
  x <- x[is.finite(x) & x > 0]   # ⬅️ IMPORTANT

  if (length(x) < 50) return(NA_real_)

  mean(vapply(
    x,
    function(y) crps_egpd_mc_pos(
      y,
      kappa = params_margins$kappa[i],
      sigma = params_margins$sigma[i],
      xi    = params_margins$xi[i],
      m = 400
    ),
    numeric(1)
  ))
})

data.frame(site = sites, crps_egpd = crps_site)


library(scoringRules)

crps_site_emp <- sapply(seq_along(sites), function(i) {
  x <- as.numeric(X_by_site[i, ])
  x <- x[is.finite(x) & x > 0]
  if (length(x) < 50) return(NA_real_)
  mean(vapply(x, function(y) crps_sample(y = y, dat = x), numeric(1)))
})

res <- data.frame(
  site = sites,
  crps_egpd = crps_site,
  crps_emp  = crps_site_emp,
  ratio     = crps_site / crps_site_emp
)

res[order(res$ratio), ]

# save in csv
foldername_crps <- paste0(
  data_folder,
  "swg/omsev/crps/"
)
if (!dir.exists(foldername_crps)) {
  dir.create(foldername_crps, recursive = TRUE)
}
filename_crps <- paste0(
  foldername_crps,
  "crps_ratio_advgroup_",
  group_adv,
  "_n", Nsim,
  ".csv"
)

write.csv(res, filename_crps, row.names = FALSE)

pit_site_pos <- function(x, kappa, sigma, xi) {
  x <- x[x > 0]
  Fegpd_pos(x, kappa, sigma, xi)
}

# Exemple site i
i <- 1
x <- as.numeric(X_by_site[i, ])
u <- pit_site_pos(
  x,
  params_margins$kappa[i],
  params_margins$sigma[i],
  params_margins$xi[i]
)


res$group <- "All sites"

ggplot(res, aes(x = group, y = ratio)) +
  geom_boxplot(
    fill = "#69b3a2",
    alpha = 0.6,
    width = 0.3
  ) +
  labs(
    y = "CRPS(EGPD) / CRPS(empirical)",
    x = ""
  ) +
  ylim(0.9, 1.2) +
  theme_minimal()


# save plot
foldername_plot_crps <- paste0(
  im_folder,
  "swg/omsev/simu/crps/"
)

if (!dir.exists(foldername_plot_crps)) {
  dir.create(foldername_plot_crps, recursive = TRUE)
}

filename_plot_crps <- paste0(
  foldername_plot_crps,
  "crps_ratio_advgroup_",
  group_adv,
  "_n", Nsim,
  ".png"
)






compute_crps_marginal_group <- function(sims_group, params_margins) {

  # sims_group : list of simulations (un seul groupe d’advection)

  sites <- rownames(sims_group[[1]]$X)

  # Stack all values for each site across sims and times
  X_by_site <- do.call(cbind, lapply(sims_group, function(sim) sim$X))

  crps_site <- sapply(seq_along(sites), function(i) {

    x <- as.numeric(X_by_site[i, ])
    x <- x[is.finite(x) & x > 0]

    if (length(x) < 50) return(NA_real_)

    mean(vapply(
      x,
      function(y) crps_egpd_mc_pos(
        y,
        kappa = params_margins$kappa[i],
        sigma = params_margins$sigma[i],
        xi    = params_margins$xi[i],
        m = 400
      ),
      numeric(1)
    ))
  })

  crps_site_emp <- sapply(seq_along(sites), function(i) {

    x <- as.numeric(X_by_site[i, ])
    x <- x[is.finite(x) & x > 0]

    if (length(x) < 50) return(NA_real_)

    mean(vapply(
      x,
      function(y) crps_sample(y = y, dat = x),
      numeric(1)
    ))
  })

  data.frame(
    site = sites,
    crps_egpd = crps_site,
    crps_emp  = crps_site_emp,
    ratio     = crps_site / crps_site_emp
  )
}


results_all_groups <- list()

for (g in all_group_names) {

  cat("Processing group:", g, "\n")

  indices_group <- which(selected_points$adv_group == g)
  if (length(indices_group) < 10) next

  list_episodes_group <- list_episodes[indices_group]
  s0_list_group <- s0_list[indices_group]
  u_list_group  <- u_list[indices_group]
  adv_group     <- adv_matrix[indices_group, ]

  Nsim <- 100
  s0_sim <- integer(Nsim)
  sims_group <- vector("list", Nsim)

  for (i in seq_len(Nsim)) {
    idx <- sample(seq_along(list_episodes_group), 1)
    s0_i <- s0_list_group[idx]
    s0_sim[i] <- s0_i

    sims_group[[i]] <- simulate_many_episodes(
      N = 1,
      u = 1000,
      u_emp = u_list_group[idx],
      params_vario = params_vario_kmh,
      params_margins = params_margins,
      coords = grid_omsev,
      times = times * 5 / 60,
      adv = adv_group[idx, ],
      t0 = 0,
      s0 = s0_i
    )[[1]]
  }

  # 🔥 CRPS marginal EXACTEMENT comme avant
  res_group <- compute_crps_marginal_group(
    sims_group,
    params_margins
  )

  res_group$group <- g
  results_all_groups[[g]] <- res_group
}


crps_all <- do.call(rbind, results_all_groups)
crps_all <- crps_all[is.finite(crps_all$ratio), ]

library(ggplot2)

ggplot(crps_all, aes(x = group, y = ratio)) +
  geom_boxplot(fill = "#69b3a2", alpha = 0.6) +
  coord_flip() +
  labs(
    x = "Advection group",
    y = "CRPS(EGPD) / CRPS(empirical)",
    title = ""
  ) +
  ylim(0.9, 1.3) +
  theme_minimal()

# save plot
foldername_plot_crps_groups <- paste0(
  im_folder,
  "swg/omsev/simu/crps/"
)

if (!dir.exists(foldername_plot_crps_groups)) {
  dir.create(foldername_plot_crps_groups, recursive = TRUE)
}

filename_plot_crps_groups <- paste0(
  foldername_plot_crps_groups,
  "crps_ratio_all_advgroups_n", Nsim,
  ".png"
)

ggsave(
  filename = filename_plot_crps_groups,
  width = 8,
  height = 6,
  dpi = 300
)
