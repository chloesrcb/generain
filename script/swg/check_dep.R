rm(list = ls())
cat("\014")

source("./script/load_libraries.R")

library(data.table)
library(dplyr)
library(tidyr)
library(purrr)
library(ggplot2)

# Get all files in the folder "R"
functions_folder <- "./R"
files <- list.files(functions_folder, full.names = TRUE)
# load all functions in files
invisible(lapply(files, function(f) source(f, echo = FALSE)))
library(latex2exp)


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
grid_coords_km <- sites_coords
grid_coords_m <- sites_coords
grid_coords_m$x_m <- (coords_m[, "X"] - min(coords_m[, "X"]))
grid_coords_m$y_m <- (coords_m[, "Y"] - min(coords_m[, "Y"]))
grid_coords_km$x_km <- (coords_m[, "X"] - min(coords_m[, "X"])) / 1000
grid_coords_km$y_km <- (coords_m[, "Y"] - min(coords_m[, "Y"]))  / 1000

# get distance matrix
grid_coords_m <- grid_coords_m[, c("x_m", "y_m")]
grid_coords_km <- grid_coords_km[, c("x_km", "y_km")]
colnames(grid_coords_m) <- c("Longitude", "Latitude")
colnames(grid_coords_km) <- c("Longitude", "Latitude")
dist_mat <- get_dist_mat(grid_coords_m, latlon = FALSE)

# Spatial chi
df_dist <- reshape_distances(dist_mat)
df_dist_km <- df_dist
df_dist_km$value <- df_dist$value / 1000

stopifnot(nrow(grid_coords_m) == ncol(rain))
rownames(grid_coords_m) <- colnames(rain)


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
s0t0_set <- get_s0t0_pairs(grid_coords_m, rain,
                            min_spatial_dist = dmin,
                            episode_size = delta,
                            set_st_excess = set_st_excess,
                            n_max_episodes = 10000,
                            latlon = FALSE)

selected_points <- s0t0_set
selected_points <- selected_points %>%
  mutate(t0_date = as.POSIXct(t0_date, format = "%Y-%m-%d %H:%M:%S", tz = "UTC"))

n_episodes <- length(selected_points$s0)
table(selected_points$s0) # number of episodes per s0
# number of episodes par t0 date
n_s0_per_t0 <- table(selected_points$t0_date)
# see the number of episodes per date in a figure

df_mult <- as.data.frame(n_s0_per_t0) %>%
    dplyr::filter(Freq > 0) %>%
    dplyr::rename(n_s0 = Freq)

ggplot(df_mult, aes(x = factor(n_s0))) +
  geom_bar(fill = btfgreen, alpha = 0.7) +
  labs(
    x = TeX("Number of episodes for the same $t_0$"),
    y = TeX("Number of $t_0$"),
    title = ""
  ) +
  theme_minimal() +
  btf_theme

# save plot
foldername_plot <- paste0(
  im_folder,
  "swg/omsev/"
)

filename_plot <- paste0(
  foldername_plot,
  "n_episodes_per_t0_",
  q * 100, "_dmin", dmin,
    "_delta", delta, ".png"
)

ggsave(
  filename = filename_plot,
  plot = last_plot(),
  width = 7,
  height = 5,
  dpi = 300
)

t0_list <- selected_points$t0
s0_list <- selected_points$s0
u_list <- selected_points$u_s0
list_episodes_points <- get_extreme_episodes(selected_points, rain,
                                     episode_size = delta, unif = FALSE)
list_episodes <- list_episodes_points$episodes
tau_vect <- 0:10
tmax <- max(tau_vect)
df_coords <- as.data.frame(grid_coords_m)
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
params_est <- c(1.090, 4.628, 0.225, 0.713, 1.621, 5.219)  # km/h
# params_est <- c(1.5, 4.4, 0.2, 0.7, 1.621, 5.219)  # km/h
# [1] 1.0871823 4.6266826 0.2238579 0.7127469 1.6210000 5.2190000
etas_estimates <- params_est[5:6]

params_kmh <- list(
  beta1 = params_est[1],
  beta2 = params_est[2],
  alpha1 = params_est[3],
  alpha2 = params_est[4]
)

beta1 <- params_kmh$beta1
beta2 <- params_kmh$beta2
alpha1 <- params_kmh$alpha1
alpha2 <- params_kmh$alpha2



##################################################################################
# TRANSFORM ADVECTION SPEEDS
##################################################################################
adv_filename <- paste(data_folder, "/omsev/adv_estim/combined_comephore_omsev/episode_advection_q",
                          q * 100, "_delta", delta, "_dmin", dmin,
                          ".csv", sep = "")
adv_df_raw <- read.csv(adv_filename, sep = ",")
head(adv_df_raw)
setDT(selected_points)
setDT(adv_df_raw)

adv_df_raw[, t0_omsev := as.POSIXct(t0_omsev, format="%Y-%m-%d %H:%M:%S", tz="UTC")]
selected_points[, t0_date := as.POSIXct(t0_date, format="%Y-%m-%d %H:%M:%S", tz="UTC")]

adv_df_t0 <- adv_df_raw[, .(
  vx_final = vx_final[1],
  vy_final = vy_final[1]
), by = t0_omsev]

setkey(adv_df_t0, t0_omsev)

# join to get the advections for selected points and 
selected_episodes <- adv_df_t0[selected_points,
  on   = .(t0_omsev = t0_date),
]

setnames(selected_episodes, c("vx_final","vy_final"), c("adv_x","adv_y"))

V_episodes <- data.frame(
  v_x = selected_episodes$adv_x,
  v_y = selected_episodes$adv_y
)


adv_df <- V_episodes

colnames(adv_df) <- c("vx_final", "vy_final")
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
direction <- atan2(adv_df$vy_final, adv_df$vx_final) * (180 / pi)
# make sure direction is in [0, 360]
angle_deg <- (direction + 360) %% 360
is_calm <- speed_raw == 0
speed_class <- as.character(cut(
  speed_raw,
  breaks = c(0, 0.5, 1, Inf),
  labels = c("still","weak","significant"),
  include.lowest = TRUE
  ))

speed_class <- factor(speed_class, levels = c("still","weak","significant"))
table(speed_class)

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
  angle_deg = angle_deg,
  direction_class = direction_class,
  speed_class = speed_class,
  group = paste0(speed_class, "_", direction_class)
)


# number of null advections
n_null_adv <- sum(adv_class$direction_class == "none")
# Wind rose plot
wind_rose_data <- adv_class %>%
  filter(direction_class != "none") %>%
  group_by(direction_class, speed_class) %>%
  summarise(count = n(), .groups = "drop")

ggplot(wind_rose_data, aes(x = direction_class, y = count, fill = speed_class)) +
  geom_col(position = "stack", width = 0.7) +
  scale_fill_manual(values = c("still" = "#fee5d9", "weak" = "#fcae91", "significant" = "#e6550d")) +
  labs(
  title = paste0("Number of null advections without direction: ", n_null_adv),
  x = "Direction",
  y = "Episodes count (advective only)",
  fill = "Speed class"
  ) +
  theme_minimal() +
  coord_polar(theta = "x", start = pi/4, direction = 1)

# save plot
foldername_plot <- paste0(
  im_folder,
  "swg/omsev/"
)

filename_plot <- paste0(
  foldername_plot,
  "adv_wind_rose_95q1200dmin12delta_NSEW.png"
)

ggsave(
  filename = filename_plot,
  plot = last_plot(),
  width = 7,
  height = 7,
  dpi = 300
)

table(adv_class$group)
selected_points$speed_class <- adv_class$speed_class
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
# Simulations
################################################################################
grid_omsev <- grid_coords_km
adv_matrix <- as.matrix(adv_df_transfo[, c("vx_t", "vy_t")])
# group_adv_names <- unique(selected_points$speed_class)
# group <- group_adv_names[3]
# idx_episodes_group <- which(selected_points$speed_class == group)
# s0_list_group <- s0_list[idx_episodes_group]
# u_list_group <- u_list[idx_episodes_group]
# adv_matrix_group <- adv_matrix[idx_episodes_group, , drop = FALSE]
M <- 10
Nsim <- length( list_episodes)  # number of episodes to simulate

# list_episodes_group <- list_episodes[idx_episodes_group]

sims_by_ep <- vector("list", Nsim)

for (j in seq_len(Nsim)) {
  ep_idx <- j
  s0_j  <- s0_list[ep_idx]
  u_j   <- u_list[ep_idx]
  adv_j <- adv_matrix[ep_idx, ]

  sims_by_ep[[j]] <- vector("list", M)

  for (m in seq_len(M)) {
    sims_by_ep[[j]][[m]] <- simulate_many_episodes(
      N = 1,
      u = 1000,
      u_emp = u_j,
      params_vario = params_kmh,
      params_margins = params_margins,
      coords = grid_omsev,
      times = times * 5 / 60,
      adv = adv_j,
      t0 = 0,
      s0 = s0_j
    )[[1]]
  }
}



sites_names <- colnames(rain)

sum_over_time_by_site <- function(ep_mat, sites_ref) { 
  # ep_mat can be [time x site] or [site x time] 
  x <- as.matrix(ep_mat) 
  if (!is.null(colnames(x)) && all(colnames(x) %in% sites_ref)) { 
    # time x site 
    s <- colSums(x, na.rm = TRUE) 
    return(s) 
  } 
  if (!is.null(rownames(x)) && all(rownames(x) %in% sites_ref)) { 
    # site x time 
    s <- rowSums(x, na.rm = TRUE) 
    return(s) 
  } 
  if (ncol(x) == length(sites_ref)) { 
      s <- colSums(x, na.rm = TRUE) 
      names(s) <- sites_ref
   } else if (nrow(x) == length(sites_ref)) { 
    s <- rowSums(x, na.rm = TRUE) 
    names(s) <- sites_ref 
    } 
  return(s) 
} 
   
   
make_qq_df <- function(x_obs, x_sim, min_n = 50) { 
  x_obs <- x_obs[is.finite(x_obs)] 
  x_sim <- x_sim[is.finite(x_sim)] 
  n <- min(length(x_obs), length(x_sim)) 
  if (n < min_n) return(NULL) 
  p <- (seq_len(n) - 0.5) / n 
  data.frame( 
    q_obs = as.numeric(quantile(x_obs, probs = p, names = FALSE)), 
    q_sim = as.numeric(quantile(x_sim, probs = p, names = FALSE)) 
  ) 
}

cum_obs_mat <- t(vapply(
  seq_along(list_episodes),
  function(i) sum_over_time_by_site(list_episodes[[i]], sites_names),
  FUN.VALUE = setNames(numeric(length(sites_names)), sites_names)
))

cum_obs_long <- as.data.frame(cum_obs_mat) |>
  mutate(ep_id = seq_len(nrow(cum_obs_mat))) |>
  pivot_longer(-ep_id, names_to = "site", values_to = "cum_obs")

cum_sim_long <- purrr::map_dfr(seq_along(sims_by_ep), function(j) {
  purrr::map_dfr(seq_along(sims_by_ep[[j]]), function(m) {
    v <- sum_over_time_by_site(sims_by_ep[[j]][[m]]$X, sites_names)
    tibble(
      ep_id = j,
      rep   = m,
      site  = names(v),
      cum_sim = as.numeric(v)
    )
  })
})

meta_ep <- as.data.frame(selected_points) |>
  mutate(ep_id = seq_len(nrow(selected_points))) |>
  select(ep_id, speed_class)
head(meta_ep)

df_cum <- cum_obs_long |>
  left_join(meta_ep, by = "ep_id") |>
  left_join(cum_sim_long, by = c("ep_id", "site"))


qq_adv <- df_cum |>
  filter(
    is.finite(cum_obs),
    is.finite(cum_sim),
    speed_class %in% c("still", "weak", "significant")
  ) |>
  group_by(speed_class) |>
  summarise(
    qq = list(make_qq_df(cum_obs, cum_sim, min_n = 300)),
    .groups = "drop"
  ) |>
  unnest(qq)


ggplot(qq_adv, aes(q_obs, q_sim)) +
  geom_point(alpha = 0.25, size = 0.9) +
  geom_abline(slope = 1, intercept = 0, linetype = 2) +
  facet_wrap(~ speed_class, scales = "free") +
  labs(
    x = "Quantiles des cumuls observés",
    y = "Quantiles des cumuls simulés",
    title = "QQ-plots des cumuls de pluie par épisode\nselon le régime d’advection"
  ) +
  theme_minimal()


df_cum_site <- df_cum |>
  group_by(ep_id, site) |>
  summarise(
    cum_obs = first(cum_obs),
    cum_sim = median(cum_sim, na.rm = TRUE),
    .groups = "drop"
  )

qq_by_site_rep <- df_cum |>
  filter(is.finite(cum_obs), is.finite(cum_sim)) |>
  group_by(site, rep) |>
  summarise(qq = list(make_qq_df(cum_obs, cum_sim, min_n = 100)), .groups="drop") |>
  unnest(qq)

ggplot(qq_by_site_rep, aes(q_obs, q_sim)) +
  geom_point(alpha = 0.12, size = 0.6, color=btfgreen) +
  geom_abline(slope = 1, intercept = 0, linetype = 2, color="red", alpha=0.5) +
  facet_wrap(~ site, scales = "free") +
  labs(
    x = "Quantiles observed cumuls",
    y = "Quantiles simulated cumuls",
  ) +
  theme_minimal()

# save plot
foldername_plot <- paste0(
  im_folder,
  "swg/omsev/qq_plots/"
)

filename_plot <- paste0(
  foldername_plot,
  "qq_plots_cumuls_by_site_allreps.png"
)

ggsave(
  filename = filename_plot,
  plot = last_plot(),
  width = 10,
  height = 8,
  dpi = 300
)

sum_over_time_by_site <- function(ep_mat, sites_ref) {

  x <- as.matrix(ep_mat)

  # ensure orientation: time x site
  if (ncol(x) != length(sites_ref) && nrow(x) == length(sites_ref)) {
    x <- t(x)
  }

  stopifnot(ncol(x) == length(sites_ref))

  apply(x, 2, function(v) {
    if (all(is.na(v))) {
      NA_real_
    } else {
      sum(v, na.rm = TRUE)
    }
  }) |> setNames(sites_ref)
}


# cum_obs_mat: [episode x site]
cum_obs_mat <- t(vapply(
  seq_along(list_episodes),
  function(i) sum_over_time_by_site(list_episodes[[i]], sites_names),
  FUN.VALUE = setNames(numeric(length(sites_names)), sites_names)
))

cum_obs_long <- as.data.frame(cum_obs_mat) %>%
  mutate(ep_id = seq_len(nrow(cum_obs_mat))) %>%
  pivot_longer(-ep_id, names_to = "site", values_to = "cum_obs")

# cum_sim_long_rep: one value per (episode, replicate, site)
cum_sim_long_rep <- purrr::map_dfr(seq_along(sims_by_ep), function(j) {
  purrr::map_dfr(seq_along(sims_by_ep[[j]]), function(m) {

    # Each replicate is assumed to contain the simulated episode matrix in $X
    epX <- sims_by_ep[[j]][[m]]$X
    epX[epX < 0.22] <- 0  # thresholding at 0.22 mm (measurement precision)
    v <- sum_over_time_by_site(epX, sites_names)

    tibble(
      ep_id   = j,
      rep     = m,
      site    = names(v),
      cum_sim = as.numeric(v)
    )
  })
})

# Aggregate over replicates: median cum_sim per (episode, site)
cum_sim_long <- cum_sim_long_rep %>%
  group_by(ep_id, site) %>%
  summarise(
    cum_sim = median(cum_sim, na.rm = TRUE),
    .groups = "drop"
  )

meta_ep <- as.data.frame(selected_points) %>%
  mutate(ep_id = seq_len(nrow(selected_points))) %>%
  select(ep_id, s0, speed_class)

# Join OBS + SIM at the (episode, site) level
df_cum <- cum_obs_long %>%
  left_join(meta_ep, by = "ep_id") %>%
  left_join(cum_sim_long, by = c("ep_id", "site"))

################################################################################
# QQ of cumuls at sites != s0, grouped by distance to s0,
# and stratified by advection regime.
################################################################################
cum_obs_mat <- t(vapply(
  seq_along(list_episodes),
  function(i) sum_over_time_by_site(list_episodes[[i]], sites_names),
  FUN.VALUE = setNames(numeric(length(sites_names)), sites_names)
))

cum_obs_long <- as.data.frame(cum_obs_mat) |>
  mutate(ep_id = seq_len(nrow(cum_obs_mat))) |>
  pivot_longer(-ep_id, names_to = "site", values_to = "cum_obs")

# 2) Simulated cumuls per episode x site x replicate
cum_sim_long_rep <- purrr::map_dfr(seq_along(sims_by_ep), function(j) {
  purrr::map_dfr(seq_along(sims_by_ep[[j]]), function(m) {
    v <- sum_over_time_by_site(sims_by_ep[[j]][[m]]$X, sites_names)
    tibble(ep_id = j, rep = m, site = names(v), cum_sim = as.numeric(v))
  })
})

cum_sim_long <- cum_sim_long_rep |>
  group_by(ep_id, site) |>
  summarise(cum_sim = median(cum_sim, na.rm = TRUE), .groups = "drop")

meta_ep <- as.data.frame(selected_points) |>
  mutate(ep_id = seq_len(nrow(selected_points))) |>
  select(ep_id, speed_class)

df_cum <- cum_obs_long |>
  left_join(meta_ep, by = "ep_id") |>
  left_join(cum_sim_long, by = c("ep_id", "site"))


df_cum_rep <- cum_obs_long %>%
  left_join(meta_ep, by = "ep_id") %>%
  left_join(cum_sim_long_rep, by = c("ep_id", "site"))

qq_by_site_rep <- df_cum_rep %>%
  filter(is.finite(cum_obs), is.finite(cum_sim)) %>%
  group_by(site, rep) %>%
  summarise(
    qq = list(make_qq_df(cum_obs, cum_sim, min_n = 100)),
    .groups = "drop"
  ) %>%
  unnest(qq)

p_site_rep <- ggplot(qq_by_site_rep, aes(q_obs, q_sim)) +
  geom_point(alpha = 0.5, size = 0.6, color = btfgreen) +
  geom_abline(slope = 1, intercept = 0, linetype = 2, color = "red", alpha = 0.5) +
  facet_wrap(~ site, scales = "free") +
  labs(
    x = "Observed cumuls (quantiles)",
    y = "Simulated cumuls (quantiles)",
  ) +
  theme_minimal() +
  btf_theme

print(p_site_rep)


foldername_plot <- paste0(im_folder, "swg/omsev/qq_plots/")
if (!dir.exists(foldername_plot)) dir.create(foldername_plot, recursive = TRUE)

ggsave(
  filename = paste0(foldername_plot, "qq_dependence_distbin_by_adv.png"),
  plot = p_dep,
  width = 10, height = 7, dpi = 300
)

ggsave(
  filename = paste0(foldername_plot, "qq_cumuls_by_adv.png"),
  plot = p_adv,
  width = 7, height = 5, dpi = 300
)

ggsave(
  filename = paste0(foldername_plot, "qq_cumuls_by_site_allreps.pdf"),
  plot = p_site_rep,
  width = 10, height = 8, dpi = 300
)



# Join obs with per-replicate sims + episode metadata (speed_class)
df_cum_rep <- cum_obs_long %>%
  left_join(meta_ep, by = "ep_id") %>%                 # adds speed_class (and s0 if present)
  left_join(cum_sim_long_rep, by = c("ep_id", "site")) %>%
  filter(is.finite(cum_obs), is.finite(cum_sim)) %>%
  filter(speed_class %in% c("still","weak", "significant"))    # optional: drop "still"

# QQ per site x rep x advection class
qq_by_site_rep_adv <- df_cum_rep %>%
  group_by(speed_class, site, rep) %>%
  summarise(
    qq = list(make_qq_df(cum_obs, cum_sim, min_n = 10)),
    .groups = "drop"
  ) %>%
  unnest(qq)

# Plot: same style as yours, but faceted by advection class
p_site_rep_adv <- ggplot(qq_by_site_rep_adv, aes(q_obs, q_sim)) +
  geom_point(alpha = 0.20, size = 0.55, color = btfgreen) +
  geom_abline(slope = 1, intercept = 0, linetype = 2, color = "red", alpha = 0.5) +
  facet_grid(speed_class ~ site, scales = "free") +
  labs(
    x = "Observed cumuls (quantiles)",
    y = "Simulated cumuls (quantiles)",
    title = "QQ-plots of episode cumuls by site (all MC replicates), stratified by advection regime"
  ) +
  theme_minimal() + btf_theme

print(p_site_rep_adv)




meta_ep2 <- as.data.frame(selected_points) %>%
  mutate(ep_id = seq_len(nrow(selected_points))) %>%
  select(ep_id, adv_group)

df_cum_rep2 <- cum_obs_long %>%
  left_join(meta_ep2, by = "ep_id") %>%
  left_join(cum_sim_long_rep, by = c("ep_id", "site")) %>%
  filter(is.finite(cum_obs), is.finite(cum_sim)) %>%
  filter(!is.na(adv_group))

qq_by_site_rep_advgroup <- df_cum_rep2 %>%
  group_by(adv_group, site, rep) %>%
  summarise(
    qq = list(make_qq_df(cum_obs, cum_sim, min_n = 10)),
    .groups = "drop"
  ) %>%
  unnest(qq)

p_site_rep_advgroup <- ggplot(qq_by_site_rep_advgroup, aes(q_obs, q_sim)) +
  geom_point(alpha = 0.18, size = 0.55, color = btfgreen) +
  geom_abline(slope = 1, intercept = 0, linetype = 2, color = "red", alpha = 0.5) +
  facet_grid(adv_group ~ site, scales = "free") +
  labs(
    x = "Observed cumuls (quantiles)",
    y = "Simulated cumuls (quantiles)",
    title = "QQ-plots by site (all MC replicates), stratified by advection group"
  ) +
  theme_minimal() + btf_theme

print(p_site_rep_advgroup)


# save plot
foldername_plot <- paste0(im_folder, "swg/omsev/qq_plots/")

ggsave(
  filename = paste0(foldername_plot, "qq_cumuls_by_site_allreps_by_advgroup.pdf"),
  plot = p_site_rep_advgroup,
  width = 14, height = 12, dpi = 300
)


for (g in c("still", "weak", "significant")) {
  p <- qq_by_site_rep_adv %>%
    filter(speed_class == g) %>%
    ggplot(aes(q_obs, q_sim)) +
    geom_point(alpha=0.5, size=0.55, color=btfgreen) +
    geom_abline(slope=1, intercept=0, linetype=2, color="red", alpha=0.5) +
    facet_wrap(~ site, scales="free") +
    labs(
         x="Observed cumuls (quantiles)", y="Simulated cumuls (quantiles)") +
    theme_minimal() + btf_theme
  print(p)
  ggsave(
    filename = paste0(foldername_plot, "qq_cumuls_by_site_allreps_", g, ".pdf"),
    plot = p,
    width = 10, height = 8, dpi = 300
  )
}

# number of simulations
n_sims <- length(sims_by_ep) * length(sims_by_ep[[1]])
n_sims