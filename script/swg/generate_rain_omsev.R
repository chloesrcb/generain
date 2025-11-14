
# remove objects from the environment
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
grid_omsev <- as.data.frame(coords_m / 1000)
colnames(grid_omsev) <- c("Longitude", "Latitude")
rownames(grid_omsev) <- rownames(sites_coords)



q <- 0.95
delta <- 12
step_min <- 5


params_vario <- list(
  beta1 = 0.0826,
  beta2 = 0.8884,
  alpha1 = 0.4344,
  alpha2 = 0.6261
)

p0_vect <- rep(0.1, ncol(rain))  # proportion of zeros
u_s <- apply(rain, 2, function(x) quantile(x[x > 0], probs = q, na.rm = TRUE))
kappa_vect <- egpd_params[match(colnames(rain), egpd_params$Site), "kappa"]
xi_vect    <- egpd_params[match(colnames(rain), egpd_params$Site), "xi"]
sigma_vect <- egpd_params[match(colnames(rain), egpd_params$Site), "sigma"]
params_margins <- list(
  xi    = xi_vect,
  sigma = sigma_vect,
  kappa = kappa_vect,
  p0    = as.numeric(p0_vect)
)

qEGPD_marg <- function(u, p0, xi, sigma, kappa) {
  u <- pmin(pmax(u, 1e-12), 1 - 1e-12)
  out <- numeric(length(u))
  idx_pos <- (u > p0)
  out[!idx_pos] <- 0.0
  if (any(idx_pos)) {
    v <- (u[idx_pos] - p0) / (1 - p0)
    out[idx_pos] <- qextgp(
      v, type = 1,
      xi = xi, sigma = sigma, kappa = kappa
    )
  }
  out
}

F_s <- function(x, p0, xi, sigma, kappa) {
  p0 + (1 - p0) * pextgp(x, type = 1, xi = xi, sigma = sigma, kappa = kappa)
}


F_s_single <- function(x, p0, xi, sigma, kappa) {
  p0 + (1 - p0) * pextgp(x, type=1, xi=xi, sigma=sigma, kappa=kappa)
}

qEGPD_full <- function(u, p0, xi, sigma, kappa) {
  # clamp du vecteur
  u <- pmin(pmax(u, 1e-12), 1 - 1e-12)

  out <- numeric(length(u))

  dry  <- (u <= p0)
  wet  <- (u > p0)

  # pluie = 0 si U ≤ p0
  out[dry] <- 0

  # seulement si on a des valeurs > p0
  if (any(wet)) {
    v <- (u[wet] - p0) / (1 - p0)  # prob dans la queue
    out[wet] <- qextgp(v, type = 1, xi = xi, sigma = sigma, kappa = kappa)
  }

  return(out)
}

sim_episode_coords <- function(params_vario, params_margins, coords, times, adv, t0, s0, u_s,
                               plot_debug = FALSE) {
  
  s0_name <- rownames(coords)[coords$Longitude == s0[1] &
                              coords$Latitude == s0[2]]
  u_s0    <- u_s[s0_name]

  sim <- sim_rpareto_coords(
    beta1 = params_vario$beta1,
    beta2 = params_vario$beta2,
    alpha1 = params_vario$alpha1,
    alpha2 = params_vario$alpha2,
    adv = adv,
    coords = coords,
    t = times,
    t0 = t0,
    s0 = s0,
    threshold = u_s0
  )

  Z <- sim$Z[,,1, drop = TRUE]

  nS <- nrow(coords)
  nT <- length(times)
  X  <- matrix(NA_real_, nS, nT, dimnames = list(rownames(coords), NULL))

  # EXTREMES
  for (k in seq_len(nS)) {
    p0    <- params_margins$p0[k]
    xi    <- params_margins$xi[k]
    sigma <- params_margins$sigma[k]
    kappa <- params_margins$kappa[k]
    u_sk <- u_s[k]
    ext <- Z[k, ] > u_sk

    if (any(ext)) {
      V_tail <- 1 - 1 / Z[k, ext]
      V_tail <- pmin(pmax(V_tail, 1e-12), 1 - 1e-12)

      F_u <- F_s(u_sk, p0, xi, sigma, kappa)
      V_tail_rescaled <- F_u + (1 - F_u) *
        (V_tail - (1 - 1/u_sk)) / (1/u_sk)
      V_tail_rescaled <- pmin(pmax(V_tail_rescaled, 1e-12), 1 - 1e-12)

      X[k, ext] <- qEGPD_full(V_tail_rescaled, p0, xi, sigma, kappa)
    }
  }

  # NON EXTREMES
  for (k in seq_len(nS)) {
    p0    <- params_margins$p0[k]
    xi    <- params_margins$xi[k]
    sigma <- params_margins$sigma[k]
    kappa <- params_margins$kappa[k]
    u_sk <- u_s[k]

    Zk <- Z[k, ]
    idx <- which(Zk <= u_sk)
    if (length(idx) == 0) next

    Z_nonext <- Zk[idx]

    N <- length(Zk)
    n_s <- length(Z_nonext)
    k_s <- floor(p0 * N)

    r_s <- rank(Z_nonext, ties.method="first")

    F_u <- F_s(u_sk, p0, xi, sigma, kappa)

    V_nonext <- numeric(n_s)

    if (k_s > 0)
      V_nonext[r_s <= k_s] <- (r_s[r_s <= k_s] / k_s) * p0

    if (k_s < n_s)
      V_nonext[r_s > k_s] <- p0 +
        (F_u - p0) * ((r_s[r_s > k_s] - k_s) / (N - k_s))

    V_nonext <- pmin(pmax(V_nonext, 1e-12), 1 - 1e-12)

    X[k, idx] <- qEGPD_full(V_nonext, p0, xi, sigma, kappa)
  }

  if (plot_debug) {
    for (s in rownames(Z)) {
      plot_transformation_gg(
        Z, X, u_s, site_name = s, u_s0 = u_s0,
        save_plot = TRUE,
        filename = paste0(im_folder, "swg/omsev/plot_transformation_", s, ".png")
      )
    }
  }
  return(X)
}





sim_episode_coords <- function(params_vario, params_margins, coords, times, adv, t0, s0, u_s,
                               plot_debug = FALSE) {

  s0_name <- rownames(coords)[coords$Longitude == s0[1] &
                              coords$Latitude == s0[2]]
  u_s0    <- u_s[s0_name]

  sim <- sim_rpareto_coords(
    beta1 = params_vario$beta1,
    beta2 = params_vario$beta2,
    alpha1 = params_vario$alpha1,
    alpha2 = params_vario$alpha2,
    adv = adv,
    coords = coords,
    t = times,
    t0 = t0,
    s0 = s0,
    threshold = u_s0
  )

  Z <- sim$Z[,,1, drop = TRUE]

  nS <- nrow(coords)
  nT <- length(times)

  X_prob <- matrix(NA_real_, nS, nT, dimnames = list(rownames(coords), NULL))

  for (k in seq_len(nS)) {
    p0    <- params_margins$p0[k]
    xi    <- params_margins$xi[k]
    sigma <- params_margins$sigma[k]
    kappa <- params_margins$kappa[k]
    u_sk  <- u_s[k]

    Zk <- Z[k, ]

    ext_idx <- which(Zk > u_sk)
    nonext_idx <- which(Zk <= u_sk)

    # EXTREMES
    if (length(ext_idx) > 0) {
      V_tail <- 1 - 1 / Zk[ext_idx]
      V_tail <- pmin(pmax(V_tail, 1e-12), 1 - 1e-12)

      F_u <- F_s(u_sk, p0, xi, sigma, kappa)
      V_tail_rescaled <- F_u + (1 - F_u) * (V_tail - (1 - 1/u_sk)) / (1/u_sk)
      V_tail_rescaled <- pmin(pmax(V_tail_rescaled, 1e-12), 1 - 1e-12)

      X_prob[k, ext_idx] <- V_tail_rescaled
    }

    # NON-EXTREMES
    if (length(nonext_idx) > 0) {
      Z_nonext <- Zk[nonext_idx]
      n_s <- length(Z_nonext)
      k_s <- floor(p0 * n_s)
      r_s <- rank(Z_nonext, ties.method = "first")

      F_u <- F_s(u_sk, p0, xi, sigma, kappa)
      V_nonext <- numeric(n_s)

      if (k_s > 0)
        V_nonext[r_s <= k_s] <- (r_s[r_s <= k_s] / k_s) * p0

      if (k_s < n_s)
        V_nonext[r_s > k_s] <- p0 + (F_u - p0) * ((r_s[r_s > k_s] - k_s) / (n_s - k_s))

      V_nonext <- pmin(pmax(V_nonext, 1e-12), 1 - 1e-12)

      X_prob[k, nonext_idx] <- V_nonext
    }
  }

  X <- matrix(NA_real_, nS, nT, dimnames = list(rownames(coords), NULL))
  for (k in seq_len(nS)) {
    X[k, ] <- qEGPD_full(X_prob[k, ], params_margins$p0[k],
                          params_margins$xi[k],
                          params_margins$sigma[k],
                          params_margins$kappa[k])
  }

  if (plot_debug) {
    for (s in rownames(Z)) {
      plot_transformation_gg(Z, X, u_s, site_name = s, u_s0 = u_s0,
                             save_plot = TRUE,
                             filename = paste0(im_folder, "swg/omsev/plot_transformation_", s, ".png"))
    }
  }

  return(X)
}






plot_transformation_gg <- function(Z, X, u_s, site_name, u_s0, 
                              save_plot = FALSE, filename = NULL) {
  
  s <- site_name
  Zs <- Z[s, ]
  Xs <- X[s, ]
  
  # Détecter extrêmes vs non-extrêmes
  ext     <- Zs > u_s[s]
  nonext  <- !ext
  
  df <- data.frame(
    time = seq_along(Zs),
    Z    = Zs,
    X    = Xs,
    type = ifelse(ext, "extreme", "nonextreme")
  )
  
  df_long <- df %>%
    pivot_longer(cols = c(Z, X), names_to = "variable", values_to = "value")
  
  # Ajouter type pour X seulement (Z reste noir)
  df_long <- df_long %>%
    mutate(plot_type = case_when(
      variable == "Z" ~ "Z",
      variable == "X" & type == "extreme" ~ "X extreme",
      variable == "X" & type == "nonextreme" ~ "X non-extreme"
    ))
  
  same_threshold <- abs(u_s0 - u_s[s]) < 1e-8  # tolérance numérique

  gg <- ggplot(df_long, aes(x = time, y = value, color = plot_type, shape = plot_type)) +
    geom_line(data = df_long %>% filter(variable == "Z"), size = 1.2) +
    geom_point(data = df_long %>% filter(variable == "X"), size = 2) +
    # seuils conditionnels
    {if(same_threshold) {
        geom_hline(yintercept = u_s0, color = "red", linetype = "dashed", size = 1)
    } else {
        list(
          geom_hline(yintercept = u_s[s], color = "red", linetype = "dashed", size = 1),
          geom_hline(yintercept = u_s0, color = "#a86a18", linetype = "dotted", size = 1)
        )
    }} +
    # annotations conditionnelles
    {if(same_threshold) {
        annotate("text", x = max(df_long$time)*0.8, y = u_s0,
                label = paste0("u_s0 = u_s = ", round(u_s0,3)),
                color = "red", vjust = -1)
    } else {
        list(
          annotate("text", x = max(df_long$time)*0.8, y = u_s[s], 
                  label = paste0("u_s = ", round(u_s[s],3)),
                  color = "red", vjust = -1),
          annotate("text", x = max(df_long$time)*0.8, y = u_s0, 
                  label = paste0("u_s0 = ", round(u_s0,3)),
                  color = "#a86a18", vjust = -1)
        )
    }} +
    scale_color_manual(values = c("Z" = "#686868", "X extreme" = "#a72909", "X non-extreme" = btfgreen)) +
    scale_shape_manual(values = c("Z" = NA, "X extreme" = 19, "X non-extreme" = 19)) +
    labs(
        x = "Time", y = "Value") +
    theme_minimal() +
    theme(legend.title = element_blank()) +
    btf_theme
  print(gg)
  if (save_plot & !is.null(filename)) {
    ggsave(filename, plot = gg, width = 12, height = 6, dpi = 300)
  }

}




# start from 2020-01-01 00:00:00 to 2024-12-31 23:55:00
start_all <- as.POSIXct("2020-01-01 00:00:00", tz = "UTC")
end_all   <- as.POSIXct("2024-12-31 23:55:00", tz = "UTC")
full_dates <- seq(from = start_all, to = end_all, by = paste(step_min, "mins"))

n_episodes <- 100

set.seed(123)
V_episodes <- data.frame(
    vx = rnorm(n_episodes, mean = 0, sd = 1),  # m/5min
    vy = rnorm(n_episodes, mean = 0, sd = 1)
)


final_mat <- matrix(0, nrow = length(full_dates), ncol = nrow(grid_omsev),
                    dimnames = list(format(full_dates), rownames(grid_omsev)))
times <- 0:(delta - 1)   # Local episode time indexing
start_id_candidates <- seq_len(nrow(final_mat) - (delta - 1))
episode_starts_id <- sample(start_id_candidates,
                            size = n_episodes,
                            replace = TRUE)

for (i in seq_len(n_episodes)) {

  s0 <- s0_list[i]
  coord_s0 <- as.numeric(grid_omsev[rownames(grid_omsev) == s0, ])

  adv_vect <- as.numeric(V_episodes[sample(nrow(V_episodes), 1), ])

  X_ep <- sim_episode_coords(
    params_vario   = params_vario,
    params_margins = params_margins,
    coords = grid_omsev,
    times  = times,
    adv    = adv_vect,
    t0     = 0,
    s0     = coord_s0,
    u_s    = u_s
  )

  idx0 <- episode_starts_id[i]
  idxs <- idx0 + times
  valid <- idxs <= nrow(final_mat)
  idxs  <- idxs[valid]

  X_ep_df_cut <- X_ep[valid, , drop = FALSE]

  final_mat[idxs, ] <- as.matrix(X_ep_df_cut)
}


final_df <- data.frame(
  time = full_dates,
  final_mat,
  check.names = FALSE
)
colSums(final_df[, -1] > 0)



sum(is.na(final_mat))
# ~0

compute_chi <- function(X, q = 0.95, eps = 1e-6) {
  # X = deux colonnes (stations) de données
  X <- X[complete.cases(X), ]
  X_s1_pos <- X[,1][X[,1] > eps]
  X_s2_pos <- X[,2][X[,2] > eps]
  u1 <- quantile(X_s1_pos, q)
  u2 <- quantile(X_s2_pos, q)

  exceed1 <- X[,1] > u1
  exceed2 <- X[,2] > u2

  chi <- sum(exceed1 & exceed2) / sum(exceed2)
  return(chi)
}

# over the same episode period
chi_real <- outer(colnames(rain), colnames(rain), Vectorize(function(s1, s2){
  compute_chi(rain[, c(s1, s2)])
}))

chi_sim <- outer(colnames(final_df)[-1], colnames(final_df)[-1], Vectorize(function(s1, s2){
  compute_chi(final_df[, c(s1, s2)])
}))
x11()
plot(chi_real, chi_sim,
     xlab="Real χ", ylab="Simulated χ", pch=19)
abline(0,1,col="red",lwd=2)



sum(final_df$iem)


# count na
for (s in 1:(ncol(final_df)-1)) {
  site <- colnames(final_df)[s+1]
  ggplot(final_df, aes(time, .data[[site]])) +
    geom_line(color=btfgreen) +
    labs(x = "Time", y = "Rainfall (mm/5 min)") +
    theme_minimal()

  filename_output <- paste0(im_folder, "swg/omsev/simulated_rain_", site, ".png")
  ggsave(filename = filename_output,
        width = 12, height = 6, dpi = 300)
}

site <- "iem"
id_site <- which(colnames(final_df) == site) - 1
# time zoom on one episode
episode_start <- episode_starts_id[1]
episode_end <- episode_start + delta - 1
episode_data <- final_df[episode_start:episode_end, ]

ggplot(episode_data, aes(time, .data[[site]])) +
  geom_line(color=btfgreen) +
  labs(x = "Time", y = "Rainfall (mm/5 min)") +
  theme_minimal()

filename_output <- paste0(im_folder, "swg/omsev/zoom_ep_simulated_rain_", site, ".png")
ggsave(filename = filename_output,
       width = 12, height = 6, dpi = 300)

# do boxplot for all sites
final_df_melt <- final_df %>%
  pivot_longer(cols = -time, names_to = "site", values_to = "rain")
final_df_melt <- final_df %>%
  pivot_longer(cols = -time, names_to = "site", values_to = "rain") %>%
  filter(rain > 0)

ggplot(final_df_melt, aes(x = site, y = rain)) +
  geom_boxplot(outlier.size = 0.5, fill = btfgreen, alpha = 0.7) +
  labs(x = "Site", y = "Rainfall (mm/5 min)") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
# save
filename_output <- paste0(im_folder, "swg/omsev/simulated_rain_boxplot_allsites_pos.png")
ggsave(filename = filename_output,
       width = 12, height = 6, dpi = 300)

final_df_melt_u <- final_df %>%
  pivot_longer(cols = -time, names_to = "site", values_to = "rain") %>%
  filter(rain > 0.6)

ggplot(final_df_melt_u, aes(x = site, y = rain)) +
  geom_boxplot(outlier.size = 0.5, fill = btfgreen, alpha = 0.7) +
  labs(x = "Site", y = "Rainfall (mm/5 min)") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
filename_output <- paste0(im_folder, "swg/omsev/simulated_rain_boxplot_allsites_above06mm.png")
ggsave(filename = filename_output,
       width = 12, height = 6, dpi = 300)


rain_time <- as.POSIXct(rownames(rain), tz = "UTC")
final_time <- as.POSIXct(final_df$time, tz = "UTC")

common_time <- intersect(rain_time, final_time)


real_site <- rain[, site]
idx <- which(!is.na(real_site))
rain_time_site <- rain_time[idx]
real_site <- real_site[idx]

common_time <- intersect(rain_time_site, final_time)


cmp_df <- tibble(
  time = common_time,
  real = real_site[match(common_time, rain_time_site)],
  sim  = final_df[match(common_time, final_time), site]
) %>%
  mutate(
    real = as.numeric(real),
    sim  = as.numeric(sim)
  )



rain_time <- as.POSIXct(rownames(rain), tz="UTC")
sites_names <- colnames(rain)
simulate_once <- function(
  full_dates,
  grid_omsev,
  delta,
  n_episodes,
  params_vario,
  params_margins,
  u_s,
  V_episodes = NULL,
  min_gap = NULL  # minimum gap between episodes
) {
  n_times <- length(full_dates)
  n_sites <- nrow(grid_omsev)
  
  final_mat <- matrix(0.0, nrow = n_times, ncol = n_sites,
                      dimnames = list(format(full_dates), rownames(grid_omsev)))
  
  # Generate episode start positions without overlap
  available_start <- 1:(n_times - delta + 1)
  if (is.null(min_gap)) min_gap <- delta  # default gap = episode duration
  episode_starts_id <- c()
  
  for (i in 1:n_episodes) {
    if (length(available_start) == 0) break
    start_i <- sample(available_start, 1)
    episode_starts_id <- c(episode_starts_id, start_i)
    
    # Remove window around this start to avoid overlap
    remove_range <- (start_i - min_gap + 1):(start_i + min_gap - 1)
    remove_range <- remove_range[remove_range >= 1 & remove_range <= n_times]
    available_start <- setdiff(available_start, remove_range)
  }
  
  # Randomly draw initial sites
  s0_list <- sample(rownames(grid_omsev), length(episode_starts_id), replace = TRUE)
  
  # If no wind field provided, generate random advection vectors
  if (is.null(V_episodes)) {
    V_episodes <- data.frame(
      vx = rnorm(length(episode_starts_id), 0, 1),
      vy = rnorm(length(episode_starts_id), 0, 1)
    )
  } else {
    V_episodes <- V_episodes[sample(1:nrow(V_episodes), length(episode_starts_id), replace = TRUE), ]
  }
  
  times <- 0:(delta - 1)
  
  for (i in seq_along(episode_starts_id)) {
    idx0 <- episode_starts_id[i]
    idxs <- idx0 + times
    idxs <- idxs[idxs <= n_times]  # avoid exceeding bounds
    
    s0 <- s0_list[i]
    coord_s0 <- as.numeric(grid_omsev[rownames(grid_omsev) == s0, ])
    
    adv_vect <- as.numeric(V_episodes[i, ])
    
    # Simulate episode
    X_ep <- sim_episode_coords(
      params_vario   = params_vario,
      params_margins = params_margins,
      coords = grid_omsev,
      times  = times,
      adv    = adv_vect,
      t0     = 0,
      s0     = coord_s0,
      u_s    = u_s
    )
    
    # Add slight multiplicative noise to avoid correlation = 1
    X_ep_noisy <- X_ep * matrix(runif(length(X_ep), 0.8, 1.2),
                                nrow = nrow(X_ep), ncol = ncol(X_ep))
    
    # Insert into final_mat taking the max with existing values
    final_mat[idxs, ] <- pmax(final_mat[idxs, ], X_ep_noisy)
  }
  
  final_df <- data.frame(
    time = full_dates,
    final_mat,
    check.names = FALSE
  )
  
  return(final_df)
}

# Example usage
set.seed(123)
B <- 10  # number of realizations
sims <- replicate(B, simulate_once(
  full_dates = full_dates,
  grid_omsev = grid_omsev,
  delta = delta,
  n_episodes = n_episodes,
  params_vario = params_vario,
  params_margins = params_margins,
  u_s = u_s,
  V_episodes = V_episodes
), simplify = FALSE)

# Store simulations in an array time × site × member
times <- sims[[1]]$time
sites <- setdiff(colnames(sims[[1]]), "time")

sim_arr <- array(NA_real_,
                 dim = c(length(times), length(sites), B),
                 dimnames = list(time = as.character(times),
                                 site = sites,
                                 member = paste0("sim", 1:B)))
for (b in seq_len(B)) {
  sim_arr[,,b] <- as.matrix(sims[[b]][, sites])
}

library(dplyr)
library(ggplot2)

quantiles_cond <- lapply(sites, function(s) {
  x_real <- rain[, s]
  x_real_pos <- x_real[x_real > 0]  # seuil de détection
  q_vect <- seq(0.9, 0.99, by = 0.01)
  q_real <- quantile(x_real_pos, probs = q_vect, na.rm = TRUE)
  
  x_sim <- sim_arr[, s, ]
  q_sim <- apply(x_sim, 2, function(v) {
    v_pos <- v[v > 0.21]
    quantile(v_pos, probs = q_vect, na.rm = TRUE)
  })
  
  data.frame(
    site = s,
    prob = q_vect,
    q_real = q_real,
    q_sim_mean = rowMeans(q_sim, na.rm = TRUE),
    q_sim_sd = apply(q_sim, 1, sd, na.rm = TRUE)
  )
}) %>% bind_rows()

df_qq <- quantiles_cond %>%
  transmute(
    site,
    prob,
    q_obs = q_real,
    q_sim = q_sim_mean,
    q_sim_lo = pmax(q_sim_mean - q_sim_sd, 0),
    q_sim_hi = q_sim_mean + q_sim_sd
  )

ggplot(df_qq, aes(x = q_obs, y = q_sim)) +
  geom_ribbon(aes(ymin = q_sim_lo, ymax = q_sim_hi), fill = "#2ca02c", alpha = 0.2) +
  geom_point(color = "#2ca02c") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  facet_wrap(~ site, scales = "free") +
  scale_x_log10() + scale_y_log10() +
  labs(x = "Observed quantiles", y = "Simulated quantiles") +
  theme_minimal()






cor_obs <- cor(rain, use = "pairwise.complete.obs")

library(reshape2)
library(ggplot2)

# Ordre des sites basé sur les données réelles
sites_order <- colnames(rain)

# Réordonner la matrice simulée
sim1_mat <- sim_arr[, sites_order, 1]

# Calculer la corrélation
cor_sim1 <- cor(sim1_mat, use = "pairwise.complete.obs")

# Transformer en format long pour ggplot
library(reshape2)
m_sim1 <- melt(cor_sim1, varnames = c("Site1", "Site2"), value.name = "Corr")
colnames(m_sim1) <- c("Site1", "Site2", "Corr")

# Faire de même pour la matrice observée
cor_obs_ordered <- cor_obs[sites_order, sites_order]
m_obs <- melt(cor_obs_ordered, varnames = c("Site1", "Site2"), value.name = "Corr")
colnames(m_obs) <- c("Site1", "Site2", "Corr")

# Heatmap Observée
ggplot(m_obs, aes(x = Site1, y = Site2, fill = Corr)) +
  geom_tile() +
  scale_fill_viridis_c() +
  theme_minimal() +
  labs(title = "Correlation observed") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Heatmap Simulée
ggplot(m_sim1, aes(x = Site1, y = Site2, fill = Corr)) +
  geom_tile() +
  scale_fill_viridis_c() +
  theme_minimal() +
  labs(title = "Correlation simulated (simulation 1)") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


# Multiple pairs scatter comparison (observed vs one simulation member)
pairs_list <- list(
  c("iem","mse"),
  c("poly","um"),
  c("cefe","cnrs"),
  c("crbm","archie"),
  c("chu1","chu2"),
  c("chu3","chu4")
)

sim_idx <- 1  # which simulation member to use

make_pair_df <- function(p, rain, sim_arr, sim_idx) {
  s1 <- p[1]; s2 <- p[2]
  obs <- rain[, c(s1, s2)]
  obs_pos <- obs[obs[,1] > 0 & obs[,2] > 0, , drop = FALSE]
  
  sim <- sim_arr[, c(s1, s2), sim_idx]
  sim_pos <- sim[sim[,1] > 0 & sim[,2] > 0, , drop = FALSE]
  
  rbind(
    data.frame(site1 = s1, site2 = s2, x = obs_pos[,1], y = obs_pos[,2], type = "obs",
               pair = paste(s1, s2, sep = "_")),
    data.frame(site1 = s1, site2 = s2, x = sim_pos[,1], y = sim_pos[,2], type = "sim",
               pair = paste(s1, s2, sep = "_"))
  )
}

df_dep_all <- do.call(rbind, lapply(pairs_list, make_pair_df, rain = rain, sim_arr = sim_arr, sim_idx = sim_idx))

ggplot(df_dep_all, aes(x = x, y = y, color = type)) +
  geom_point(alpha = 0.35, size = 1) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  facet_wrap(~ pair, scales = "free") +
  labs(x = "Site 1", y = "Site 2") +
  theme_minimal()
# save
filename_output <- paste0(im_folder, "swg/omsev/dependence_scatter_obs_vs_sim.png")
ggsave(filename = filename_output,
       width = 12, height = 6, dpi = 300)
# Optional: compute tail dependence chi for each pair (same function as before)
compute_chi <- function(X, q = 0.95, eps = 1e-6) {
  X <- X[complete.cases(X), ]
  X1 <- X[,1][X[,1] > eps]; X2 <- X[,2][X[,2] > eps]
  if (length(X1) < 5 || length(X2) < 5) return(NA_real_)
  u1 <- quantile(X1, q); u2 <- quantile(X2, q)
  exceed1 <- X[,1] > u1; exceed2 <- X[,2] > u2
  if (sum(exceed2) == 0) return(NA_real_)
  sum(exceed1 & exceed2) / sum(exceed2)
}

chi_df <- lapply(pairs_list, function(p) {
  s1 <- p[1]; s2 <- p[2]
  chi_obs <- compute_chi(rain[, c(s1, s2)])
  chi_sim <- compute_chi(sim_arr[, c(s1, s2), sim_idx])
  data.frame(pair = paste(s1, s2, sep = "_"), chi_obs = chi_obs, chi_sim = chi_sim)
}) %>% bind_rows()

print(chi_df)

# plot
ggplot(chi_df, aes(x = chi_obs, y = chi_sim)) +
  geom_point(color = btfgreen, size = 3) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  labs(x = "Tail dependence χ observed", y = "Tail dependence χ simulated") +
  theme_minimal()




# i <- 1
# simu_df <- sims[[i]]
# rownames(simu_df) <- simu_df$time
# simu_df <- simu_df[ , -which(colnames(simu_df) == "time")]


# create_generator_gif_points(
#   simu_df = simu_df[3200:3500, ],
#   coords = location_gauges,
#   outfile = paste0("simulation_omsev_ep", i, ".gif"),
#   interval = 0.5,
#   forcedtemp = 1000
# )




# extract_sim_episodes <- function(final_df, s0, delta, n_ep = 5) {
  
#   times <- final_df$time
#   site_series <- final_df[, s0]
#   s0_excesses <- which(site_series > 0.65)
#   idx_peaks <- sample(s0_excesses, n_ep, replace = FALSE)
#   idx_peaks <- sort(idx_peaks)
  
#   all_ep <- list()
  
#   for (idx0 in idx_peaks) {
#     idxs <- idx0:(idx0 + delta - 1)
#     idxs <- idxs[idxs <= nrow(final_df)]
    
#     ep <- final_df[idxs, ]
#     rownames(ep) <- ep$time
#     ep$time <- NULL
    
#     all_ep[[length(all_ep) + 1]] <- ep
#   }
  
#   ep_total <- do.call(rbind, all_ep)
  
#   return(ep_total)  # matrice [ (n_ep*delta) × sites ]
# }


# s0 <- "poly"         # site de conditionnement visuel
# delta <- 12          # taille épisode
# n_ep <- 5            # nombre d'épisodes à animer

# ep_sim_multi <- extract_sim_episodes(final_df, s0, delta, n_ep)

# coords_s0 <- location_gauges[location_gauges$Station == s0, c("Longitude","Latitude")]
# colnames(coords_s0) <- c("x","y")

# create_generator_gif_points(
#   simu_df = ep_sim_multi,
#   coords = location_gauges,
#   outfile = paste0("sim_multi_", s0, ".gif"),
#   interval = 0.5,
#   forcedtemp = delta * n_ep,
#   s0 = coords_s0
# )




# extract_sim_episodes_multi_s0 <- function(final_df, s0_list, delta) {
  
#   times <- final_df$time
#   all_ep <- list()
  
#   for (s0 in s0_list) {
    
#     site_series <- final_df[, s0]
    
#     # point d’ancrage : le maximum simulé pour ce site
#     idx0 <- which.max(site_series)
    
#     # index temporels de l’épisode
#     idxs <- idx0:(idx0 + delta - 1)
#     idxs <- idxs[idxs <= nrow(final_df)]  # sécurité bord fin de série
    
#     ep <- final_df[idxs, ]
#     rownames(ep) <- ep$time
#     ep$time <- NULL
    
#     all_ep[[length(all_ep) + 1]] <- ep
#   }
  
#   # empiler tous les épisodes verticalement
#   ep_total <- do.call(rbind, all_ep)
#   return(ep_total)
# }


# delta <- 12
# s0_list_multi <- sample(sites, 8)  # par ex: choisir 8 sites au hasard

# ep_sim_multi <- extract_sim_episodes_multi_s0(final_df, s0_list_multi, delta)


# # On place un marqueur sur le site courant (facultatif)
# create_generator_gif_points(
#   simu_df = ep_sim_multi,
#   coords = location_gauges,
#   outfile = "simulation_multi_s0.gif",
#   interval = 0.5,
#   forcedtemp = delta * length(s0_list_multi)
# )


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

params_vario <- list(
  beta1 = 0.3,
  beta2 = 0.78,
  alpha1 = 0.22,
  alpha2 = 0.70
)


q <- 0.95
delta <- 12
dmin <- 5  # km

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

list_episodes_points <- get_extreme_episodes(selected_points, rain,
                                     episode_size = delta, unif = FALSE)
list_episodes <- list_episodes_points$episodes

adv_filename <- paste(data_folder, "/omsev/adv_estim/combined_comephore_omsev/episode_advection_q",
                          q * 100, "_delta", delta, "_dmin", min_spatial_dist,
                          ".csv", sep = "")
adv_df <- read.csv(adv_filename, sep = ",")
head(adv_df)
episode_idx <- 50
s0 <- s0_list[episode_idx]
t0 <- t0_list[episode_idx]
adv_vect <- as.numeric(adv_df[episode_idx, c("vx_final","vy_final")])

# count proportion of zeros per site in episode
p0_vect <- sapply(sites_names, function(s) {
  ep_data <- list_episodes[[episode_idx]][, s]
  p0 <- mean(ep_data == 0, na.rm = TRUE)
  if (is.nan(p0)) p0 <- 0
  return(p0)
})


u_s <- apply(rain, 2, function(x) quantile(x[x > 0], probs = q, na.rm = TRUE))
kappa_vect <- egpd_params[match(colnames(rain), egpd_params$Site), "kappa"]
xi_vect    <- egpd_params[match(colnames(rain), egpd_params$Site), "xi"]
sigma_vect <- egpd_params[match(colnames(rain), egpd_params$Site), "sigma"]
params_margins <- list(
  xi    = xi_vect,
  sigma = sigma_vect,
  kappa = kappa_vect,
  p0    = as.numeric(p0_vect)
)
# statio quand je ferais sur le tout

sim_ep <- sim_episode_coords(
  params_vario   = params_vario,
  params_margins = params_margins,
  coords  = grid_omsev,
  times   = 0:(delta-1),
  adv     = adv_vect,
  t0      = 0,
  s0      = as.numeric(grid_omsev[rownames(grid_omsev) == s0, ]),
  u_s     = u_s
)

obs_ep <- list_episodes[[episode_idx]]
times_ep <- 0:(delta-1)

library(ggplot2)
site <- "iem" 
df_plot <- data.frame(
  time = times_ep,
  real = obs_ep[, site],
  sim  = sim_ep[site, ]
)

ggplot(df_plot, aes(x = time)) +
  geom_line(aes(y = real, color = "observed")) +
  geom_line(aes(y = sim, color = "simulated")) +
  labs(x = "Time step", y = "Rainfall (mm/5min)",
       color = "Legend") +
  theme_minimal()

# save plot
filename_output <- paste0(im_folder, "swg/omsev/episode_", episode_idx, "_site_", site, ".png")
ggsave(filename = filename_output,
       width = 12, height = 6, dpi = 300)
