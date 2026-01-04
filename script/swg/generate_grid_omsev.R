
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
library(viridis)
library(magick)


# Get all files in the folder "R"
functions_folder <- "./R"
files <- list.files(functions_folder, full.names = TRUE)
# load all functions in files
invisible(lapply(files, function(f) source(f, echo = FALSE)))
library(latex2exp)


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


pEGPD_full <- function(x, p0, xi, sigma, kappa) {
    v <- numeric(length(x))
    v[x <= 0] <- 0
    v[x == 0] <- p0
    idx <- which(x > 0)
    if(length(idx) > 0) {
      v[idx] <- p0 + (1 - p0) * pextgp(x[idx], xi = xi, sigma = sigma, kappa = kappa)
    }
    v <- pmin(pmax(v, 1e-12), 1 - 1e-12)
    return(v)
}

qEGPD_full <- function(v, p0, xi, sigma, kappa) {
  v <- pmin(pmax(v, 1e-12), 1 - 1e-12)

  out <- numeric(length(v))
  
  na_idx <- is.na(v)
  out[na_idx] <- NA
  
  dry  <- (v <= p0) & !na_idx
  wet  <- (v > p0) & !na_idx

  # rain = 0 if U ≤ p0
  out[dry] <- 0

  # if > p0
  if (any(wet)) {
    x <- (v[wet] - p0) / (1 - p0)  # tail prob
    out[wet] <- qextgp(x, type = 1, xi = xi, sigma = sigma, kappa = kappa)
  }

  return(out)
}

# qEGPD_full <- function(v, p0, xi, sigma, kappa) {
#   # v : matrix n_sites x n_times
#   nS <- nrow(v)
#   nT <- ncol(v)
#   out <- matrix(0, nS, nT)
  
#   for (k in seq_len(nS)) {
#     vk <- v[k, ]
#     pk <- p0[k]
#     xi_k <- xi[k]
#     sigma_k <- sigma[k]
#     kappa_k <- kappa[k]
    
#     dry  <- vk <= pk
#     wet  <- vk > pk
    
#     out[k, dry] <- 0
#     if(any(wet)) {
#       x <- (vk[wet] - pk) / (1 - pk)
#       out[k, wet] <- qextgp(x, type = 1, xi = xi_k, sigma = sigma_k, kappa = kappa_k)
#     }
#   }
  
#   return(out)
# }

# pEGPD_full <- function(x, p0, xi, sigma, kappa) {
#   nS <- nrow(x)
#   nT <- ncol(x)
#   v <- matrix(0, nS, nT)
  
#   for (k in seq_len(nS)) {
#     xk <- x[k, ]
#     pk <- p0[k]
#     xi_k <- xi[k]
#     sigma_k <- sigma[k]
#     kappa_k <- kappa[k]
    
#     idx_pos <- xk > 0
#     v[k, !idx_pos] <- 0
#     if(any(idx_pos)) {
#       v[k, idx_pos] <- pk + (1 - pk) * pextgp(xk[idx_pos], type = 1, xi = xi_k, sigma = sigma_k, kappa = kappa_k)
#     }
#   }
  
#   v <- pmin(pmax(v, 1e-12), 1 - 1e-12)
#   return(v)
# }


plot_transformation_gg <- function(Z, X, u, site_name,
                              save_plot = FALSE, filename = NULL) {
  
  s <- site_name
  Zs <- Z[s, ]
  Xs <- X[s, ]
  
  ext     <- Zs > u
  interm <- (Zs > 0) & (Zs <= u)
  low  <- !ext & !interm
  
  df <- data.frame(
    time = seq_along(Zs),
    Z    = Zs,
    X    = Xs,
    type = case_when(
      ext      ~ "extreme",
      interm   ~ "intermediate",
      low      ~ "low"
    )
  )
  
  df_long <- df %>%
    pivot_longer(cols = c(X), names_to = "variable", values_to = "value")
  
  df_long <- df_long %>%
    mutate(plot_type = case_when(
      variable == "X" & type == "extreme" ~ "Z extreme",
      variable == "X" & type == "intermediate" ~ "Z intermediate",
      variable == "X" & type == "low" ~ "Z low"
    ))
  

  gg <- ggplot(df_long, aes(x = time, y = value, color = plot_type, shape = plot_type)) +
    geom_point(data = df_long %>% filter(variable == "X"), size = 2) +
    scale_color_manual(values = c("Z extreme" = "#a72909", "Z intermediate" = "#f4a261", "Z low" = btfgreen)) +
    scale_shape_manual(values = c("Z extreme" = 19, "Z intermediate" = 19, "Z low" = 19)) +
    labs(
        x = "Time", y = "X value") +
    theme_minimal() +
    theme(legend.title = element_blank()) +
    btf_theme
  
  print(gg)
  if (save_plot & !is.null(filename)) {
    ggsave(filename, plot = gg, width = 12, height = 6, dpi = 300)
  }

}

# G_std <- function(z, p0, u) {
#   v <- numeric(length(z))

#   v[z < 0] <- 0
#   x0 <- 2/(1 - p0)

#   idx1 <- which(z > 0 & z <= x0)
#   v[idx1] <- p0

#   idx2 <- which(z > x0 & z <= u)
#   v[idx2] <- p0 + ( (1 - 1/u) - p0 ) / (u - x0) * (z[idx2] - x0)

#   idx3 <- which(z > u)
#   v[idx3] <- 1 - 1/z[idx3]

#   v <- pmin(pmax(v, 1e-12), 1 - 1e-12)

#   return(v)
# }

# G_std_inv <- function(v, p0, u) {
#   v <- pmin(pmax(v, 1e-12), 1 - 1e-12)
  
#   z <- numeric(length(v))
#   z[v <= 0] <- -Inf
  
#   idx0 <- which(v > 0 & v <= p0)
#   if (length(idx0) > 0) z[idx0] <- 0
  
#   v_switch <- 1 - 1 / u
  
#   x0 <- 2/(1 - p0)
  
#   idx_mid <- which(v > p0 & v <= v_switch)
#   if (length(idx_mid) > 0) {
#     z[idx_mid] <- (4 / (1 - p0)^2) * (v[idx_mid] - p0)
#   }
  
#   idx_high <- which(v > v_switch)
#   if (length(idx_high) > 0) z[idx_high] <- 1 / (1 - v[idx_high])
  
#   return(z)
# }



G_std <- function(z, p0, u) {
  v <- numeric(length(z))
  v[] <- NA_real_
  z <- as.numeric(z)
  x0 <- 2 / (1 - p0)
  
  # z < 0
  idx_neg <- which(z < 0)
  if (length(idx_neg)) v[idx_neg] <- 0
  
  # 0 <= z <= x0 : G(z) = p0 + a*z
  a <- (1 - p0)^2 / 4
  idx1 <- which(z >= 0 & z <= x0)
  if (length(idx1)) v[idx1] <- p0 + a * z[idx1]
  
  # x0 < z <= u : linear interpolation between v_x0 and 1 - 1/u
  idx2 <- which(z > x0 & z <= u)
  if (length(idx2)) {
    v_x0 <- p0 + a * x0   # = (1 + p0)/2
    v_u  <- 1 - 1 / u
    # if u == x0 then no interval; treat numerically
    if (abs(u - x0) < .Machine$double.eps) {
      v[idx2] <- v_u
    } else {
      slope2 <- (v_u - v_x0) / (u - x0)
      v[idx2] <- v_x0 + slope2 * (z[idx2] - x0)
    }
  }
  
  # region z > u : tail
  idx3 <- which(z > u)
  if (length(idx3)) v[idx3] <- 1 - 1 / z[idx3]
  
  # clamps numeriques
  v <- pmin(pmax(v, 1e-12), 1 - 1e-12)
  return(v)
}


G_std_inv <- function(v, p0, u) {
  v <- pmin(pmax(v, 1e-12), 1 - 1e-12)
  z <- numeric(length(v))
  x0 <- 2 / (1 - p0)
  a <- (1 - p0)^2 / 4
  v_x0 <- p0 + a * x0  # = (1 + p0)/2
  v_switch <- 1 - 1 / u
  
  z[v <= 0] <- -Inf
  
  idx_mass0 <- which(v > 0 & v <= p0)
  if (length(idx_mass0)) z[idx_mass0] <- 0
  
  idx_low <- which(v > p0 & v <= v_x0)
  if (length(idx_low)) z[idx_low] <- (v[idx_low] - p0) / a
  
  idx_mid <- which(v > v_x0 & v <= v_switch)
  if (length(idx_mid)) {
    if (abs(u - x0) < .Machine$double.eps) {
      z[idx_mid] <- u
    } else {
      slope2 <- (v_switch - v_x0) / (u - x0)
      z[idx_mid] <- x0 + (v[idx_mid] - v_x0) / slope2
    }
  }
  
  idx_high <- which(v > v_switch)
  if (length(idx_high)) z[idx_high] <- 1 / (1 - v[idx_high])
  
  return(z)
}


sim_episode_coords <- function(params_vario, params_margins,
                               coords, times, adv, t0, s0,
                               u, u_emp,
                               plot_debug = FALSE, filename = NULL) {
  idx_s0 <- which(names(params_margins$p0) == s0)
  p0_s0 <- params_margins$p0[idx_s0]
  xi_s0 <- params_margins$xi[idx_s0]
  sigma_s0 <- params_margins$sigma[idx_s0]
  kappa_s0 <- params_margins$kappa[idx_s0]
  x_s0 <- pEGPD_full(u_emp,
                      p0 = p0_s0,
                      xi = xi_s0,
                      sigma = sigma_s0,
                      kappa = kappa_s0)
  u <- G_std_inv(x_s0, p0 = p0_s0, u = u)
  print(u)
  s0_coords <- as.numeric(coords[rownames(coords) == s0, ])
  sim <- sim_rpareto_coords(
    beta1 = params_vario$beta1,
    beta2 = params_vario$beta2,
    alpha1 = params_vario$alpha1,
    alpha2 = params_vario$alpha2,
    adv    = adv,
    coords = coords,
    t      = times,
    t0     = t0,
    s0     =  s0_coords,
    threshold = u
  )
  
  Z <- sim$Z[,,1, drop = TRUE]
  nS <- nrow(coords)
  nT <- length(times)
  
  X <- matrix(NA_real_, nS, nT,
              dimnames = list(rownames(coords), NULL))
  V <- matrix(NA_real_, nS, nT,
              dimnames = list(rownames(coords), NULL))
  extreme_idx <- matrix(FALSE, nS, nT)
  for (k in seq_len(nS)) {

    p0    <- params_margins$p0[k]
    xi    <- params_margins$xi[k]
    sigma <- params_margins$sigma[k]
    kappa <- params_margins$kappa[k]

    Zk <- Z[k, ]
    
    V[k, ] <- G_std(Zk, p0 = p0, u = u)

    X[k, ] <- qEGPD_full(V[k, ], p0, xi, sigma, kappa)

  }

  if (plot_debug) {
    for (s in rownames(Z)) {
      if(is.null(filename)) {
        filename <- "plot_transformation_"
      } 
      plot_transformation_gg(Z, X, u, site_name = s,
                             save_plot = TRUE,
                             filename = paste0(im_folder, "swg/omsev/", filename, s, ".png"))
    }
  }
  
  return(X)
}

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
  beta1 = 0.3125,
  beta2 = 0.7805,
  alpha1 = 0.2293,
  alpha2 = 0.7016
)

# u_s <- apply(rain, 2, function(x) quantile(x[x > 0], probs = q, na.rm = TRUE))

q <- 0.95
delta <- 12
dmin <- 1200  # m

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

adv_filename <- paste(data_folder, "/omsev/adv_estim/combined_comephore_omsev/episode_advection_q",
                          q * 100, "_delta", delta, "_dmin", dmin,
                          ".csv", sep = "")
adv_df <- read.csv(adv_filename, sep = ",")
head(adv_df)

# convert adv to m/5min from km/h
adv_df$vx_final <- adv_df$vx_final * 1000 / 60 / 12
adv_df$vy_final <- adv_df$vy_final * 1000 / 60 / 12


episode_idx <- 230
real_episode <- list_episodes[[episode_idx]]
# remove columns with all NA
real_episode <- real_episode[, colSums(is.na(real_episode)) < nrow(real_episode)]
s0_real <- s0_list[episode_idx]
t0_real <- t0_list[episode_idx]
u_s0_real <- u_list[episode_idx]
adv_real <- as.numeric(adv_df[episode_idx, c("vx_final", "vy_final")])

coord_s0_real <- as.numeric(grid_omsev[rownames(grid_omsev) == s0_real, ])

ep_matrix <- list_episodes[[episode_idx]]
ep_values <- unlist(ep_matrix, use.names = FALSE)
ep_values <- ep_values[!is.na(ep_values)]


episode <- list_episodes[[episode_idx]]

compute_p0_episode <- function(episode) {
  out <- numeric(ncol(episode))
  names(out) <- colnames(episode)

  for(s in colnames(episode)) {

    x <- episode[, s]

    if (all(is.na(x))) {
      out[s] <- NA
    } else {
      out[s] <- mean(x == 0, na.rm = TRUE)
    }
  }

  return(out)
}

compute_p0_all_episodes <- function(list_episodes) {

  list_p0 <- lapply(list_episodes, compute_p0_episode)
  mat_p0  <- do.call(rbind, list_p0)

  p0_mean <- colMeans(mat_p0, na.rm = TRUE)

  return(list(
    p0_by_episode = mat_p0,
    p0_mean = p0_mean
  ))
}

res <- compute_p0_all_episodes(list_episodes)

p0_episode_obs <- res$p0_by_episode
p0_mean_by_site <- res$p0_mean
p0_target <- p0_mean_by_site

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
# count number of 0 by sites in the real episode
p0_real_episode <- compute_p0_episode(real_episode)
params_margins <- list(
  xi    = xi_vect,
  sigma = sigma_vect,
  kappa = kappa_vect,
  p0    = p0_values
)


simulate_many_episodes <- function(N, u, u_emp, params_vario, params_margins,
                                   coords, times, adv, t0, s0, plot_debug = FALSE) {
  sims <- vector("list", N)
  for (i in seq_len(N)) {
    # print progress every 50
    if (i %% 50 == 0) cat("Sim", i, " / ", N, "\n")
    Xsim <- sim_episode_coords(
      params_vario   = params_vario,
      params_margins = params_margins,
      coords         = coords,
      times          = times,
      adv            = adv,
      t0             = t0,
      s0             = s0,
      u              = u,
      u_emp           = u_emp,
      plot_debug     = FALSE
    )
    sims[[i]] <- Xsim
  }
  return(sims)
}

# plot simu_1 at site "cefe"
site <- "cefe"
plot(simu_1[, site], type = "l", col = "red", lwd = 2,
     xlab = "Time step", ylab = "Rainfall (mm/5min)",
     main = paste0("Overall simulated series at site ", site))


sim_episode_grid <- function(params_vario, params_margins_common,
                        coords, times, adv, t0, s0_pixel_id,
                        u, u_emp,
                        plot_debug = FALSE, filename = NULL) {

  s0_coords <- as.numeric(coords[rownames(coords) == s0_pixel_id, ])

  x_s0 <- pEGPD_full(u_emp,
                     p0    = params_margins_common$p0,
                     xi    = params_margins_common$xi,
                     sigma = params_margins_common$sigma,
                     kappa = params_margins_common$kappa)
  u <- G_std_inv(x_s0, p0 = params_margins_common$p0, u = u)

  sim <- sim_rpareto_coords(
    beta1 = params_vario$beta1,
    beta2 = params_vario$beta2,
    alpha1 = params_vario$alpha1,
    alpha2 = params_vario$alpha2,
    adv    = adv,
    coords = coords,
    t      = times,
    t0     = t0,
    s0     = s0_coords,
    threshold = u
  )

  Z <- sim$Z[,,1, drop = TRUE]
  nS <- nrow(coords)
  nT <- length(times)

  X <- matrix(NA_real_, nS, nT, dimnames = list(rownames(coords), NULL))
  V <- matrix(NA_real_, nS, nT, dimnames = list(rownames(coords), NULL))

  for (k in seq_len(nS)) {
    Zk <- Z[k, ]
    V[k, ] <- G_std(Zk, p0 = params_margins_common$p0, u = u)
    X[k, ] <- qEGPD_full(V[k, ],
                         p0    = params_margins_common$p0,
                         xi    = params_margins_common$xi,
                         sigma = params_margins_common$sigma,
                         kappa = params_margins_common$kappa)
  }

  return(X)
}


library(sf)

# Convertir sites en sf et définir bounding box
sites_sf <- st_as_sf(sites_coords, coords = c("Longitude", "Latitude"), crs = 2154)
bbox <- st_bbox(sites_sf)
bbox_sf <- st_as_sfc(st_bbox(c(xmin = 3.848, ymin = 43.628,
                               xmax = 3.870, ymax = 43.639), crs = 4326))


# Transformer en projection métrique (Lambert 93)
bbox_m <- st_transform(bbox_sf, crs = 2154)
bbox <- st_bbox(bbox_m)

# Créer grille régulière
x_seq <- seq(bbox$xmin, bbox$xmax, by = 100)  # 100 m
y_seq <- seq(bbox$ymin, bbox$ymax, by = 100)

grid_sf <- expand.grid(x = x_seq, y = y_seq) %>%
  st_as_sf(coords = c("x", "y"), crs = 2154)
grid_coords <- st_coordinates(grid_sf)
grid_df <- as.data.frame(grid_coords)
rownames(grid_df) <- paste0("pixel_", seq_len(nrow(grid_df)))
colnames(grid_df) <- c("Longitude", "Latitude")


grid_sf <- st_as_sf(grid_df, coords = c("Longitude", "Latitude"), crs = 2154)
# add pixels id
grid_sf$pixel_id <- rownames(grid_df)
grid_latlon <- st_transform(grid_sf, crs = 4326)
grid_latlon_coords <- st_coordinates(grid_latlon)

grid_df <- as.data.frame(grid_latlon_coords)
rownames(grid_df) <- grid_sf$pixel_id
colnames(grid_df) <- c("Longitude", "Latitude")

# plot on a map pixel and rain gauges
library(ggplot2)
library(ggspatial)  

sites_sf <- st_as_sf(sites_coords, coords = c("Longitude", "Latitude"), crs = 4326)
ggplot() +
  geom_sf(data = bbox_sf, fill = "lightblue", alpha = 0.3) +
  geom_sf(data = grid_latlon, color = "grey80", size = 0.5) +
  geom_sf(data = sites_sf, color = "red", size = 2) +
  theme_minimal() +
  btf_theme

# Grille en polygones
res <- 100  # 100 m

grid_l93_poly <- st_make_grid(bbox, cellsize = res, what = "polygons")
grid_l93_poly <- st_sf(pixel_id = paste0("pixel_", seq_along(grid_l93_poly)),
                       geometry = grid_l93_poly)

# Transformer en lat/lon pour ggplot
grid_latlon_poly <- st_transform(grid_l93_poly, crs = 4326)

# Plot
ggplot() +
  geom_sf(data = grid_latlon_poly, fill = NA, color = "grey80", size = 0.3) +
  geom_sf(data = sites_sf, color = "red", size = 2) +
  theme_minimal() +
  theme(
    axis.title = element_blank(),
    axis.ticks = element_blank(),
    panel.grid = element_blank()
  ) 

# save this plot
ggsave(paste0(im_folder, "swg/omsev/grid_pixels_rain_gauges.png"),
       width = 8, height = 6, dpi = 300)


p0_common    <- mean(p0_values)
xi_common    <- mean(xi_vect)
sigma_common <- mean(sigma_vect)
kappa_common <- mean(kappa_vect)

params_margins_common <- list(
  xi    = xi_common,
  sigma = sigma_common,
  kappa = kappa_common,
  p0    = p0_common
)
s0_real <- s0_list[episode_idx]
idx_s0 <- which(rownames(grid_omsev) == s0_real)
s0_coords <- grid_omsev[idx_s0, ]
s0_coords <- data.frame(
  Longitude = s0_coords[1],
  Latitude = s0_coords[2]
)
s0_sf <- st_as_sf(s0_coords, coords = c("Longitude", "Latitude"), crs = 2154)

# Calculer distances entre s0 et chaque pixel
distances <- st_distance(s0_sf, grid_sf)

# Identifier le pixel le plus proche
closest_idx <- which.min(distances)
closest_pixel <- grid_sf[closest_idx, ]

print(closest_pixel)

s0_pixel_id <- closest_pixel$pixel_id

simu_df <- sim_episode_grid(
  params_vario   = params_vario,
  params_margins_common = params_margins_common,
  coords = grid_df,
  times  = 0:(delta - 1),
  adv    = as.numeric(adv_real),
  t0     = 0,
  s0_pixel_id = closest_pixel$pixel_id,
  u = 100,
  u_emp = u_s0_real
)


rain_min <- min(simu_df, na.rm = TRUE)
rain_max <- max(simu_df, na.rm = TRUE)

transfo <- "sqrt"

img_files <- c()

obs_list <- lapply(1:nrow(real_episode), function(t) {
  data.frame(
    site = colnames(real_episode),
    lon  = sites_coords[colnames(real_episode), "Longitude"],
    lat  = sites_coords[colnames(real_episode), "Latitude"],
    rainfall = as.numeric(real_episode[t, ]), 
    time = t
  )
})

obs_df <- do.call(rbind, obs_list)

head(obs_df)


for (t in 1:ncol(simu_df)) {

  rain_t <- simu_df[, t]
  df_plot <- data.frame(
    pixel_id = rownames(simu_df),
    X = grid_df$Longitude,
    Y = grid_df$Latitude,
    rainfall = rain_t
  )
  
  simu_poly <- merge(grid_latlon_poly, df_plot, by = "pixel_id")
  obs_t <- obs_df[obs_df$time == t, ]  # sous-ensemble du temps t

  p <- ggplot() +
    geom_sf(data = simu_poly, aes(fill = rainfall), color = NA) +
    geom_point(data = obs_t, aes(x = lon, y = lat, color = rainfall), 
               size = 3) +  # couleur selon la pluie observée
    scale_fill_viridis_c(
      option = "C",
      limits = c(rain_min, rain_max),
      trans = transfo,
      oob = scales::squish
    ) +
    scale_color_viridis_c(option = "C", limits = c(rain_min, rain_max)) +
    labs(
      title = paste0("Simulated episode – time step ", t),
      fill  = "Rainfall (mm/5min)",
      color = "Observed rain"
    ) +
    theme_minimal() +
    theme(
      axis.text = element_blank(),
      axis.ticks = element_blank()
    )

  filename <- paste0(im_folder, "swg/omsev/grid_simulated_episode_t", t, ".png")
  ggsave(filename, plot = p, width = 8, height = 6, dpi = 200)
  img_files <- c(img_files, filename)
}



for (t in 1:ncol(simu_df)) {

  rain_t <- simu_df[, t]

  df_plot <- data.frame(
    pixel_id = rownames(simu_df),
    X = grid_df$Longitude,
    Y = grid_df$Latitude,
    rainfall = rain_t
  )

  simu_poly <- merge(grid_latlon_poly, df_plot, by = "pixel_id")

  p <- ggplot() +
    geom_sf(data = simu_poly, aes(fill = rainfall), color = NA) +
    scale_fill_viridis_c(
      option = "C",
      limits = c(rain_min, rain_max),
      trans = transfo,
      oob = scales::squish
    ) +
    labs(
      title = paste0("Simulated episode – time step ", t),
      fill  = "Rainfall (mm/5min)"
    ) +
    theme_minimal() +
    theme(
      axis.text = element_blank(),
      axis.ticks = element_blank()
    )

  filename <- paste0(im_folder, "swg/omsev/grid_simulated_episode_t", t, ".png")
  ggsave(filename, plot = p, width = 8, height = 6, dpi = 200)
  img_files <- c(img_files, filename)
}

imgs <- image_read(img_files)
gif  <- image_animate(imgs, fps = 2)  # ajuster la vitesse
gif_path <- paste0(im_folder, "swg/omsev/simulated_episode.gif")
image_write(gif, path = gif_path)

cat("GIF créé :", gif_path, "\n")



plot(adv_df$vx_final, adv_df$vy_final, pch = 19,
     xlab = "vx", ylab = "vy")
speed <- sqrt(adv_df$vx_final^2 + adv_df$vy_final^2)
quantile(speed)
speed <- sqrt(adv_df$vx_final^2 + adv_df$vy_final^2)
plot(speed)

calm_threshold <- 0.001

is_calm <- speed < calm_threshold

selected_points$adv_speed <- speed
angle <- atan2(adv_df$vy_final, adv_df$vx_final)
quantile(angle)
selected_points$adv_angle <- angle
speed_class <- ifelse(
  is_calm,
  "calm",
  as.character(cut(
    speed,
    breaks = c(0.001, 0.2, 1.1, Inf),
    labels = c("low", "moderate", "high")
  ))
)

speed_class <- factor(speed_class, levels = c("calm","low","moderate","high"))

angle <- atan2(adv_df$vy_final, adv_df$vx_final)
angle_deg <- (angle * 180/pi + 360) %% 360

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

adv_class <- data.frame(speed_class, direction_class)
table(adv_class)
adv_class$group <- ifelse(
  speed_class == "calm",
  "calm",
  paste(speed_class, direction_class, sep = ".")
)
table(adv_class$group)


selected_points$adv_group <- adv_class$group
head(selected_points)


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

# concat adv_class and init_class
adv_init_class <- data.frame(
  adv_group = adv_class$group,
  init_class = selected_points$init_class
)
adv_init_class$group <- paste(adv_init_class$adv_group,
                              adv_init_class$init_class,
                              sep = "_")

table(adv_init_class$group)

selected_points$adv_init_group <- adv_init_class$group
head(selected_points)

# Clustering of episodes with similar advection in list_episodes

adv_matrix <- as.matrix(adv_df[, c("vx_final", "vy_final")])
all_group_names <- unique(selected_points$adv_group)
group_adv <- all_group_names[10]  # choose one group to simulate
#get list_episodes and s0_list for this group
indices_group <- which(selected_points$adv_group == group_adv)
list_episodes_group <- list_episodes[indices_group]
s0_list_group <- s0_list[indices_group]
u_list_group <- u_list[indices_group]
adv_group <- adv_matrix[indices_group, ]

Nsim <- 1000
sims_group <- vector("list", Nsim)
s0_sim <- integer(Nsim)

for (i in seq_len(Nsim)) {
  idx <- sample(seq_along(list_episodes_group), 1)
  episode_i <- list_episodes_group[[idx]]
  s0_i <- s0_list_group[idx]
  s0_sim[i] <- s0_i
  # X_s0_t0 <- selected_points$X_s0_t0[indices_group[idx]]
  # cat("Simulating ", i, "/", Nsim,
  #     " (episode idx:", indices_group[idx],
  #     " s0:", s0_i,
  #     " X_s0_t0:", round(X_s0_t0,2), ")\n")
  
  sims_group[[i]] <- simulate_many_episodes(
    N = 1,
    u = 1000,
    u_emp = u_list_group[idx],
    params_vario = params_vario,
    params_margins = params_margins,
    coords = grid_omsev,
    times = times,
    adv = adv_group,
    t0 = 0,
    s0 = s0_i
  )[[1]]
}

Nsim <- 1000
m_sim_ep <- ceiling(Nsim / length(list_episodes_group))
for (idx in seq_along(list_episodes_group)) {
  # idx <- sample(seq_along(list_episodes_group), 1)
  episode_i <- list_episodes_group[[idx]]
  s0_i <- s0_list_group[idx]
  s0_sim[i] <- s0_i
  X_s0_t0 <- selected_points$X_s0_t0[indices_group[idx]]
  sims_group[[i]] <- simulate_many_episodes(
    N = m_sim_ep,
    u = 1000,
    u_emp = X_s0_t0,
    params_vario = params_vario,
    params_margins = params_margins,
    coords = grid_omsev,
    times = times,
    adv = adv_group,
    t0 = 0,
    s0 = s0_i
  )[[1]]
}


df_all_sims_group <- do.call(rbind, lapply(seq_len(Nsim), function(i) {
  data.frame(
    time = seq_len(ncol(sims_group[[i]])),
    rainfall = sims_group[[i]][s0_sim[i], ],
    sim_id = i
  )
}))


df_real_group <- do.call(rbind, lapply(indices_group, function(ep_id) {
  ep <- list_episodes[[ep_id]]
  s0_i <- s0_list[ep_id]
  data.frame(
    time = seq_len(nrow(ep)),
    rainfall = ep[, s0_i],
    episode_id = ep_id
  )
}))


ggplot() +
  geom_line(
    data = df_all_sims_group,
    aes(x = time, y = rainfall, group = sim_id, color = "Simulations"),
    alpha = 0.15
  ) +
  geom_line(
    data = df_real_group,
    aes(x = time, y = rainfall, group = episode_id, color = "Observed episodes"),
    alpha = 0.7, size = 0.8
  ) +
  scale_color_manual(
    values = c(
      "Simulations" = "#979595",
      "Observed episodes" = "#c73535"
    )
  ) +
  labs(
    x = "Time step",
    y = "Rainfall (mm/5min)",
    color = ""
  ) +
  theme_minimal() +
  btf_theme

# save plot
foldername <- paste0(im_folder, "swg/omsev/group_adv_", group_adv, "/")
if(!dir.exists(foldername)) {
  dir.create(foldername, recursive = TRUE)
}

ggsave(paste0(foldername, "all_sims_group_", group_adv, "_s0_varied.png"),
       width = 12, height = 6, dpi = 300)

df_summ_group <- df_all_sims_group %>%
  group_by(time) %>%
  summarize(
    median = median(rainfall),
    q25 = quantile(rainfall, 0.25),
    q95 = quantile(rainfall, 0.95)
  )

ggplot() +
  geom_ribbon(
    data = df_summ_group,
    aes(x = time, ymin = q25, ymax = q95),
    fill = "#979595", alpha = 0.3
  ) +
  geom_line(
    data = df_summ_group,
    aes(x = time, y = median),
    color = "#979595", size = 1
  ) +
  geom_line(
    data = df_real_group,
    aes(x = time, y = rainfall, group = episode_id),
    color = "#c73535", alpha = 0.7, size = 0.8
  ) +
  labs(x = "Time step", y = "Rainfall (mm/5min)") +
  theme_minimal() +
  btf_theme

# save plot
ggsave(paste0(foldername, "CI_95_all_sims_group_", group_adv, "_s0_varied_summary.png"),
       width = 12, height = 6, dpi = 300)


df_summ_group <- df_all_sims_group %>%
  group_by(time) %>%
  summarize(
    median = median(rainfall),
    q25 = quantile(rainfall, 0.025),
    q95 = quantile(rainfall, 0.975)
  )

ggplot() +
  geom_ribbon(
    data = df_summ_group,
    aes(x = time, ymin = q25, ymax = q95),
    fill = "#979595", alpha = 0.3
  ) +
  geom_line(
    data = df_summ_group,
    aes(x = time, y = median),
    color = "#979595", size = 1
  ) +
  geom_line(
    data = df_real_group,
    aes(x = time, y = rainfall, group = episode_id),
    color = "#c73535", alpha = 0.7, size = 0.8
  ) +
  labs(x = "Time step", y = "Rainfall (mm/5min)") +
  theme_minimal() +
  btf_theme

# save plot
ggsave(paste0(foldername, "CI_99_all_sims_group_", group_adv, "_s0_varied_summary.png"),
       width = 12, height = 6, dpi = 300)


# do CRPS calculation for each simulation in the group
library(scoringRules)

# max of all episodes of all groups
max_Xs0_t0 <- max(selected_points$X_s0_t0)
max_idx <- which(selected_points$X_s0_t0 == max_Xs0_t0)
episode_idx <- max_idx[1]
ep_idx <- 105
real_episode <- list_episodes[[ep_idx]]
s0_real <- s0_list[ep_idx]
u_s0_real <- u_list[ep_idx]
adv_real <- as.numeric(adv_df[ep_idx, c("vx_final", "vy_final")])
df_long_real <- data.frame(
  time = seq_len(nrow(real_episode)),
  observed = real_episode[, s0_real],
  site = s0_real
)

# Do 1000 simulations for this episode
X_s0_t0 <- real_episode[1, s0_real]
# deltaXs0t0 <- 0.5*X_s0_t0

Nsim <- 1000
sims_list <- simulate_many_episodes(
  N = Nsim,
  u = 1000,
  u_emp = X_s0_t0,
  params_vario = params_vario,
  params_margins = params_margins,
  coords = grid_omsev,
  times = times,
  adv = as.numeric(adv_real),
  t0 = 0,
  s0 = s0_real
)

# plot all simu with ggplot at site s0_real

# plot all simu with ggplot at site s0_real
df_all_sims <- data.frame(
  time = rep(seq_len(nrow(real_episode)), Nsim),
  rainfall = unlist(lapply(sims_list, function(x) x[s0_real, ])),
  sim_id = rep(seq_len(Nsim), each = nrow(real_episode))
)

s <- s0_real

df_all_sims <- data.frame(
  time = rep(seq_len(nrow(real_episode)), Nsim),
  rainfall = unlist(lapply(sims_list, function(x) x[s, ])),
  sim_id = rep(seq_len(Nsim), each = nrow(real_episode))
)

ggplot() +
  geom_line(data = df_all_sims, aes(x = time, y = rainfall, group = sim_id, color = "Simulations"),
            alpha = 0.2) +
  geom_line(data = df_long_real %>% filter(site == s),
            aes(x = time, y = observed, color = "Observed"),
            size = 1) +
  scale_color_manual(values = c("Observed" = "#c73535", "Simulations" = "#979595")) +
  labs(
    x = "Time step",
    y = "Rainfall (mm/5min)",
    color = ""
  ) +
  theme_minimal() +
  btf_theme

# save this plot
ggsave(paste0(im_folder, "swg/omsev/all_sims_episode_", episode_idx, "_", s, "_s0_", s0_real, ".png"),
       width = 12, height = 6, dpi = 300)

library(dplyr)

df_summ <- df_all_sims %>%
  group_by(time) %>%
  summarize(
    median = median(rainfall),
    q05 = quantile(rainfall, 0.25),
    q95 = quantile(rainfall, 0.975)
  )


df_long_real <- data.frame(
  time = seq_len(nrow(real_episode)),
  observed = real_episode[, s0_real],
  site = s0_real
)
ggplot() +
  geom_ribbon(data = df_summ, aes(x = time, ymin = q05, ymax = q95),
              fill = "#979595", alpha = 0.3) +
  geom_line(data = df_summ, aes(x = time, y = median), color = "#979595", size = 1) +
  geom_line(data = df_long_real %>% filter(site == s),
            aes(x = time, y = observed),
            color = "#c73535", size = 1) +
  labs(x = "Time step", y = "Rainfall (mm/5min)") +
  theme_minimal() +
  btf_theme

ggsave(paste0(im_folder, "swg/omsev/sim_vs_real/95CI_episode_", episode_idx, "_", s, "_s0_", s0_real, ".png"),
         width = 12, height = 6, dpi = 300)

# plot one 
ggplot() +
  geom_line(data = df_all_sims %>% filter(sim_id == 100),
            aes(x = time, y = rainfall, color = "Simulation"),
            alpha = 0.8, size = 1.5) +
  geom_line(data = df_long_real %>% filter(site == s0_real),
            aes(x = time, y = observed, color = "Observed"),
            size = 1) +
  scale_color_manual(values = c("Simulation" = btfgreen, "Observed" = "#c73535")) +
  labs(
    x = "Time step",
    y = "Rainfall (mm/5min)",
    color = ""
  ) +
  theme_minimal() +
  btf_theme

# save this plot
ggsave(paste0(im_folder, "swg/omsev/1sim_episode_", episode_idx, "_", s, "_s0_", s0_real, ".png"),
       width = 12, height = 6, dpi = 300)



library(scoringRules)

# Nombre de pas de temps
n_time <- max(df_real_group$time)

# Pour chaque pas de temps, récupérer toutes les simulations et l'observation (moyenne si plusieurs épisodes)
crps_group <- sapply(seq_len(n_time), function(t) {
  sims_t <- df_all_sims_group %>%
    filter(time == t) %>%
    pull(rainfall)
  
  obs_t <- df_real_group %>%
    filter(time == t) %>%
    pull(rainfall)
  
  # Si plusieurs épisodes réels, prendre la moyenne
  obs_t_mean <- mean(obs_t)
  
  crps_sample(y = obs_t_mean, dat = sims_t)
})

# Résumé du CRPS
mean_crps_group <- mean(crps_group)
median_crps_group <- median(crps_group)

mean_crps_group
median_crps_group


library(scoringRules)
library(dplyr)

n_sims <- length(unique(df_all_sims_group$sim_id))
n_time <- max(df_all_sims_group$time)

crps_within_sims <- matrix(NA, nrow = n_time, ncol = n_sims)

for (i in 1:n_sims) {
  for (t in 1:n_time) {
    sim_i <- df_all_sims_group %>% filter(sim_id == i, time == t) %>% pull(rainfall)
    sims_other <- df_all_sims_group %>% filter(sim_id != i, time == t) %>% pull(rainfall)
    crps_within_sims[t, i] <- crps_sample(y = sim_i, dat = sims_other)
  }
}

# Moyenne sur les simulations et le temps
mean_crps_within <- mean(crps_within_sims)
median_crps_within <- median(crps_within_sims)

mean_crps_within
median_crps_within

n_time <- max(df_real_group$time)

crps_vs_obs <- sapply(seq_len(n_time), function(t) {
  sims_t <- df_all_sims_group %>%
    filter(time == t) %>%
    pull(rainfall)
  
  obs_t <- df_real_group %>%
    filter(time == t) %>%
    pull(rainfall)
  
  # Moyenne des épisodes réels si plusieurs
  obs_t_mean <- mean(obs_t)
  
  crps_sample(y = obs_t_mean, dat = sims_t)
})

# Résumé
mean_crps_vs_obs <- mean(crps_vs_obs)
median_crps_vs_obs <- median(crps_vs_obs)

mean_crps_vs_obs
median_crps_vs_obs


library(ggplot2)
library(dplyr)

# Résumé des simulations
df_summ_group <- df_all_sims_group %>%
  group_by(time) %>%
  summarize(
    median = median(rainfall),
    q025 = quantile(rainfall, 0.025),
    q975 = quantile(rainfall, 0.975)
  )

# Observations moyennes par pas de temps
df_obs_group <- df_real_group %>%
  group_by(time) %>%
  summarize(observed = mean(rainfall))

# Calcul CRPS vs observations
library(scoringRules)
crps_group <- sapply(seq_len(nrow(df_summ_group)), function(t) {
  sims_t <- df_all_sims_group %>% filter(time == t) %>% pull(rainfall)
  obs_t <- df_obs_group$observed[t]
  crps_sample(y = obs_t, dat = sims_t)
})

# Graphique
ggplot() +
  geom_ribbon(data = df_summ_group, aes(x = time, ymin = q025, ymax = q975),
              fill = "#979595", alpha = 0.3) +
  geom_line(data = df_summ_group, aes(x = time, y = median),
            color = "#979595", size = 1) +
  geom_line(data = df_obs_group, aes(x = time, y = observed),
            color = "#c73535", size = 1) +
  geom_line(aes(x = seq_len(nrow(df_summ_group)), y = crps_group*max(df_summ_group$median)/max(crps_group)),
            color = "blue", linetype = "dashed") +
  scale_y_continuous(
    name = "Rainfall (mm/5min)",
    sec.axis = sec_axis(~ . * max(crps_group) / max(df_summ_group$median),
                        name = "CRPS")
  ) +
  labs(
    x = "Time step",
    y = "Rainfall (mm/5min)",
    title = "Simulations vs Observations avec CRPS",
    subtitle = "Ribben: 95% CI simulations, rouge: observations, ligne bleue: CRPS"
  ) +
  theme_minimal()
