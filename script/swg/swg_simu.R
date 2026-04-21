library(sf)
library(dplyr)
library(ggplot2)
library(magick)
library(grid)

sim_rpareto_coords_m5 <- function(coords,
                                  steps,
                                  beta1, beta2, alpha1, alpha2,
                                  adv_m5 = c(0, 0), threshold = 1,
                                  s0_index = 1, t0_index = 1, seed = NULL) {

  if (!is.null(seed)) set.seed(seed)

  RandomFields::RFoptions(
    spConform = FALSE,
    allow_duplicated_locations = TRUE,
    install = "no"
  )

  coords <- as.data.frame(coords)
  n_sites <- nrow(coords)
  lt <- length(steps)

  # coords must be in meters (Lambert-93)
  x_coords <- coords$x
  y_coords <- coords$y

  # temporal index in number of 5-min steps
  step_grid <- expand.grid(site = seq_len(n_sites), it = seq_len(lt))

  x <- x_coords[step_grid$site]
  y <- y_coords[step_grid$site]
  step_t <- steps[step_grid$it]

  # advection in meters per 5-min step
  x_shift <- x - step_t * adv_m5[1]
  y_shift <- y - step_t * adv_m5[2]

  # Conditioning point
  ind0 <- which(step_grid$site == s0_index & step_grid$it == t0_index)
  x0s <- x_shift[ind0]
  y0s <- y_shift[ind0]
  t0_step <- step_t[ind0]

  # Models
  modelSpace <- RandomFields::RMfbm(alpha = alpha1, var = 2 * beta1, scale = 1)
  modelTime  <- RandomFields::RMfbm(alpha = alpha2, var = 2 * beta2, scale = 1)

  # gamma0(s,t) = gamma_space + gamma_time
  gamma_space_vec <- RandomFields::RFvariogram(
    modelSpace,
    x = x_shift - x0s,
    y = y_shift - y0s
  )

  gamma_time_vec <- RandomFields::RFvariogram(
    modelTime,
    x = step_t - t0_step
  )

  gamma0_vec <- gamma_space_vec + gamma_time_vec
  gamma0 <- matrix(gamma0_vec, nrow = n_sites, ncol = lt, byrow = FALSE)

  # simulate W_s on unique shifted coordinates
  key <- paste0(sprintf("%.8f", x_shift), "_", sprintf("%.8f", y_shift))
  key_u <- unique(key)
  map_idx <- match(key, key_u)

  xy_u <- do.call(rbind, strsplit(key_u, "_", fixed = TRUE))
  x_u <- as.numeric(xy_u[, 1])
  y_u <- as.numeric(xy_u[, 2])

  W_s_u <- RandomFields::RFsimulate(modelSpace, x = x_u, y = y_u, grid = FALSE)
  W_s <- W_s_u[map_idx]

  # simulate W_t on temporal steps
  W_t <- RandomFields::RFsimulate(modelTime, steps, grid = TRUE)

  W <- W_s + W_t[step_grid$it]
  W_mat <- matrix(W, nrow = n_sites, ncol = lt, byrow = FALSE)

  # r-Pareto transform
  W0 <- W_mat[s0_index, t0_index]
  Y  <- exp(W_mat - W0 - gamma0)
  R  <- evd::rgpd(n = 1, loc = 1, scale = 1, shape = 1)

  Z <- threshold * R * Y
  rownames(Z) <- rownames(coords)

  list(Z = Z, W = W_mat, gamma0 = gamma0, R = R)
}


sim_episode_grid_m5 <- function(params_vario, params_margins_common,
                                coords, steps, adv_m5,
                                t0, s0_pixel_id, u_emp) {

  s0_index <- which(rownames(coords) == s0_pixel_id)

  x_s0 <- pEGPD_full(
    u_emp,
    p0    = params_margins_common$p0,
    xi    = params_margins_common$xi,
    sigma = params_margins_common$sigma,
    kappa = params_margins_common$kappa
  )

  u <- G_std_inv(x_s0, p0 = params_margins_common$p0)

  sim <- sim_rpareto_coords_m5(
    beta1 = params_vario$beta1,
    beta2 = params_vario$beta2,
    alpha1 = params_vario$alpha1,
    alpha2 = params_vario$alpha2,
    adv_m5 = adv_m5,
    coords = coords,
    steps = steps,
    t0_index = t0 + 1,
    s0_index = s0_index,
    threshold = u
  )

  Z <- sim$Z
  nS <- nrow(coords)
  nT <- length(steps)

  X <- matrix(NA_real_, nS, nT, dimnames = list(rownames(coords), NULL))
  V <- matrix(NA_real_, nS, nT, dimnames = list(rownames(coords), NULL))

  for (k in seq_len(nS)) {
    Zk <- Z[k, ]
    V[k, ] <- G_std(Zk, p0 = params_margins_common$p0)
    X[k, ] <- qEGPD_full(
      V[k, ],
      p0    = params_margins_common$p0,
      xi    = params_margins_common$xi,
      sigma = params_margins_common$sigma,
      kappa = params_margins_common$kappa
    )
  }

  X
}


library(sf)
library(dplyr)
library(ggplot2)
library(magick)
library(grid)

s0_pixel_id <- "pixel_200"
u_emp <- 1


# pour un demi-pixel de 1 km, prendre plutôt qqch comme c(300, 400) ou c(500, 0)

# number of time steps of 5 min
nT <- 12
steps <- 0:(nT - 1)

# grid_df assumed initially in lon/lat
grid_pts_l93 <- st_as_sf(grid_df, coords = c("Longitude", "Latitude"), crs = 4326) %>%
  st_transform(2154)

grid_pts_xy <- st_coordinates(grid_pts_l93)
grid_df_l93 <- data.frame(
  x = grid_pts_xy[, 1],
  y = grid_pts_xy[, 2]
)
rownames(grid_df_l93) <- rownames(grid_df)

# polygons and gauges in Lambert-93
grid_poly_l93 <- st_transform(grid_latlon_poly, 2154)
sites_l93 <- st_transform(sites_sf, 2154)
s0_poly_l93 <- grid_poly_l93 %>% filter(pixel_id == s0_pixel_id)

# optional static plot
p_grid <- ggplot() +
  geom_sf(data = grid_poly_l93, fill = NA, color = "grey80", linewidth = 0.3) +
  geom_sf(data = sites_l93, color = "red", size = 2) +
  geom_sf(data = s0_poly_l93, fill = NA, color = "blue", linewidth = 1.2) +
  theme_void() +
  theme(
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA)
  )



ggsave(
  filename = paste0(im_folder, "swg/omsev/grid_pixels_rain_gauges_l93.png"),
  plot = p_grid,
  width = 8, height = 6, dpi = 300, bg = "white"
)

# advection in meters per 5-min step
adv_m5 <- c(50, -30)
sim_episode <- sim_episode_grid_m5(
  params_vario = params_m5min,
  params_margins_common = params_margins_common,
  coords = grid_df_l93,
  steps = steps,
  adv_m5 = adv_m5,
  t0 = 0,
  s0_pixel_id = s0_pixel_id,
  u_emp = u_emp
)

threshold <- u_emp

positive_vals <- sim_episode[is.finite(sim_episode) & sim_episode > 0]
fill_max <- if (length(positive_vals) > 0) quantile(positive_vals, 0.995, na.rm = TRUE) else threshold
fill_max <- max(fill_max, 2 * threshold)
fill_limits <- c(0, fill_max)

episode_idx <- 1
dir_frames <- file.path(im_folder, paste0("swg/omsev/frames_episode_", episode_idx))
dir.create(dir_frames, recursive = TRUE, showWarnings = FALSE)

# centroids in L93
grid_cent_l93 <- st_centroid(grid_poly_l93)
cent_xy <- st_coordinates(grid_cent_l93)

cent_df <- data.frame(
  pixel_id = grid_poly_l93$pixel_id,
  x = cent_xy[, 1],
  y = cent_xy[, 2]
)

dx_m <- adv_m5[1]
dy_m <- adv_m5[2]

for (tt in seq_len(nT)) {

  df_t <- data.frame(
    pixel_id = rownames(sim_episode),
    rain = sim_episode[, tt]
  )

  map_t <- grid_poly_l93 %>%
    left_join(df_t, by = "pixel_id") %>%
    mutate(
      extreme = rain > threshold,
      rain_plot = rain
    )

  tmp <- df_t %>%
    left_join(cent_df, by = "pixel_id") %>%
    mutate(w = pmax(rain, 0))

  has_rain <- any(tmp$w > threshold / 4, na.rm = TRUE)

  p <- ggplot() +
    geom_sf(data = map_t, aes(fill = rain_plot), color = NA) +
    geom_sf(data = s0_poly_l93, fill = NA, color = "#ee8686", linewidth = 1.2) +
    geom_sf(data = sites_l93, color = "white", size = 3) +
    geom_sf(data = sites_l93, color = "grey40", size = 2) +
    scale_fill_gradientn(
      colours = c("white", "#dbe9f6", "#9ecae1", "#4a90d9", "#08519c"),
      values = scales::rescale(c(0, 0.2, threshold, 2 * threshold, fill_limits[2])),
      limits = fill_limits,
      oob = scales::squish,
      na.value = "transparent",
      name = paste0("Rainfall (mm/5min)\nThreshold u = ", round(threshold, 2))
    ) +
    labs(
      title = paste0(
        "t = ", tt,
        "   |   advection per 5 min = (", dx_m, ", ", dy_m, ") m"
      )
    ) +
    theme_void() +
    theme(
      plot.title = element_text(size = 18, hjust = 0.5),
      plot.background = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = NA),
      legend.background = element_rect(fill = "white", color = NA),
      legend.key = element_rect(fill = "white", color = NA),
      legend.text = element_text(size = 14)
    ) +
    coord_sf(expand = FALSE)

  if (has_rain) {
    bx <- sum(tmp$x * tmp$w, na.rm = TRUE) / sum(tmp$w, na.rm = TRUE)
    by <- sum(tmp$y * tmp$w, na.rm = TRUE) / sum(tmp$w, na.rm = TRUE)

    bary_pt_l93 <- st_sfc(st_point(c(bx, by)), crs = 2154) %>% st_sf()

    arrow_line_l93 <- st_sfc(
      st_linestring(rbind(
        c(bx, by),
        c(bx + dx_m, by + dy_m)
      )),
      crs = 2154
    ) %>% st_sf()

    p <- p +
      geom_sf(data = bary_pt_l93, color = "white", size = 2) +
      geom_sf(data = bary_pt_l93, color = "red", size = 1) +
      geom_sf(
        data = arrow_line_l93,
        arrow = arrow(length = unit(0.35, "cm")),
        color = "white",
        linewidth = 2
      ) +
      geom_sf(
        data = arrow_line_l93,
        arrow = arrow(length = unit(0.3, "cm")),
        color = "red",
        linewidth = 1
      )
  }

  out_png <- file.path(dir_frames, sprintf("frame_%03d.png", tt))
  ggsave(out_png, plot = p, width = 7, height = 6, dpi = 150, bg = "transparent")
}

options(str = NULL)
files <- list.files(dir_frames, pattern = "frame_\\d+\\.png$", full.names = TRUE)
files <- sort(files)

img <- image_read(files)
gif <- image_animate(img, fps = 1)

out_gif <- file.path(im_folder, paste0("swg/omsev/simulated_episode_", episode_idx, ".gif"))
image_write(gif, path = out_gif)

out_gif
