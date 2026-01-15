# remove objects from the environment
rm(list = ls())

# clear console
cat("\014")

muse <- FALSE

if (muse) {
  folder_muse <- "/home/serrec/work_rainstsimu/rpareto/"
  setwd(folder_muse)
  source("load_libraries.R")
  im_folder <- "./images"
  source("config_com.R")
  data_folder <- "./data/"
  ncores <- 27
} else {
  source("./script/load_libraries.R")
  source("./script/optimisation/config_com.R")
  ncores <- detectCores() - 1
}

# Get all files in the folder "R"
functions_folder <- "./R"
files <- list.files(functions_folder, full.names = TRUE)
invisible(lapply(files, function(f) source(f, echo = FALSE)))

library(latex2exp)
library(sf)
library(dplyr)
library(ggplot2)
library(gridExtra)
library(parallel)

# LOAD DATA ####################################################################
filename_com <- paste0(data_folder, "comephore/comephore_2008_2024_5km.csv")
comephore_raw <- read.csv(filename_com, sep = ",")
filename_loc <- paste0(data_folder, "comephore/loc_px_zoom_5km.csv")
loc_px <- read.csv(filename_loc, sep = ",")

loc_px <- loc_px[loc_px$pixel_name %in% colnames(comephore_raw)[-1], ]
rownames(loc_px) <- NULL

df_comephore <- as.data.frame(comephore_raw)
df_comephore$date <- as.POSIXct(df_comephore$date,
                                format = "%Y-%m-%d %H:%M:%S", tz = "UTC")
rownames(df_comephore) <- format(df_comephore$date, "%Y-%m-%d %H:%M:%S")
comephore <- df_comephore[-1]

# DISTANCE AND COORDS ##########################################################
nsites <- nrow(loc_px)
sites_coords <- loc_px[, c("Longitude", "Latitude")]
rownames(sites_coords) <- loc_px$pixel_name
sites_coords_sf <- st_as_sf(sites_coords, coords = c("Longitude", "Latitude"),
                            crs = 4326)
sites_coords_sf <- st_transform(sites_coords_sf, crs = 2154)
coords_m <- st_coordinates(sites_coords_sf)
grid_coords_km <- as.data.frame(coords_m / 1000)
colnames(grid_coords_km) <- c("Longitude", "Latitude")
rownames(grid_coords_km) <- rownames(sites_coords)

# Spatial chi WLSE #############################################################
foldername <- paste0(data_folder, "/comephore/WLSE/")
df_spa_result <- read.csv(paste0(foldername, "wlse_spa_summary.csv"))
df_temp_result <- read.csv(paste0(foldername, "wlse_temp_summary.csv"))
df_result_all <- merge(df_spa_result, df_temp_result, by = NULL)

df_result_all <- df_result_all[, c("q_spa", "q_temp",
                                   "beta1", "alpha1", "beta2", "alpha2")]

df_result_all$beta1 <- round(df_result_all$beta1, 4)
df_result_all$alpha1 <- round(df_result_all$alpha1, 4)
df_result_all$beta2 <- round(as.numeric(df_result_all$beta2), 4)
df_result_all$alpha2 <- round(as.numeric(df_result_all$alpha2), 4)

# # compute U_taylor for each quantile
# df_result_all$L <- (1 / df_result_all$beta1) ^ (1 / df_result_all$alpha1)
# df_result_all$D <- (1 / df_result_all$beta2) ^ (1 / df_result_all$alpha2)
# df_result_all$U_Taylor <- df_result_all$L / df_result_all$D

df_result <- df_result_all[df_result_all$q_spa == q &
                           df_result_all$q_temp == q, ]

beta1 <- df_result$beta1
beta2 <- df_result$beta2
alpha1 <- df_result$alpha1
alpha2 <- df_result$alpha2


# CHOOSE EXTREME EPISODE FOR R-PARETO ##########################################
set_st_excess <- get_spatiotemp_excess(data = comephore, quantile = q,
                                      remove_zeros = TRUE)
first_ts <- as.POSIXct(rownames(comephore)[1], tz = "UTC")
starting_year <- year(first_ts)
starting_year_suffix <- if (starting_year == 2008) "" else paste0("_from", starting_year)

list_s <- set_st_excess$list_s
list_t <- set_st_excess$list_t
list_u <- set_st_excess$list_u

for (i in seq_along(list_s)) {
  s0 <- list_s[[i]]
  t0 <- list_t[[i]][1]
  u_s0 <- list_u[[i]][1]
  if (comephore[t0, s0] <= u_s0) {
    stop("Excess is not above threshold for s =", s0, " and t =", t0)
  }
}




library(data.table)
library(dplyr)
library(ggplot2)

dmins <- c(3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15)  # in km

episode_pair_stats <- function(sel, sites_coords, dmin_km,
                 delta_steps,
                 latlon = FALSE, beta = 0) {

  if (nrow(sel) < 2) {
  return(list(
    n_episodes = nrow(sel),
    n_pairs = 0,
    n_pairs_dtlt_delta = 0,
    p_close_space_50km = NA_real_,
    p_close_space_120km = NA_real_,
    p_close_space_50km_given_dtlt_delta = NA_real_,
    p_conflict_dmin_delta = NA_real_,
    spatial_d10 = NA_real_,
    temporal_d10 = NA_real_
  ))
  }

  dist_matrix <- get_dist_mat(sites_coords, latlon = latlon)

  s <- sel$s0
  t <- sel$t0

  idx <- combn(seq_along(s), 2)
  s1 <- s[idx[1, ]]; s2 <- s[idx[2, ]]
  t1 <- t[idx[1, ]]; t2 <- t[idx[2, ]]

  spatial_d <- dist_matrix[cbind(s1, s2)]

  temporal_steps <- abs(t1 - t2)
  temporal_hours <- temporal_steps

  delta_eff_steps <- delta_steps + 2 * beta
  is_dt_close <- temporal_steps < delta_eff_steps

  p_close_50_all <- mean(spatial_d < 50, na.rm = TRUE)
  p_close_120_all <- mean(spatial_d < 120, na.rm = TRUE)

  n_pairs_dt <- sum(is_dt_close, na.rm = TRUE)
  p_close_50_given_dt <- if (n_pairs_dt == 0) NA_real_ else
  mean((spatial_d < 50)[is_dt_close], na.rm = TRUE)

  p_conflict <- mean((spatial_d < dmin_km) & is_dt_close, na.rm = TRUE)

  list(
  n_episodes = nrow(sel),
  n_pairs = length(spatial_d),
  n_pairs_dtlt_delta = n_pairs_dt,
  p_close_space_50km = p_close_50_all,
  p_close_space_120km = p_close_120_all,
  p_close_space_50km_given_dtlt_delta = p_close_50_given_dt,
  p_conflict_dmin_delta = p_conflict,
  spatial_d10 = as.numeric(quantile(spatial_d, 0.10, na.rm = TRUE)),
  temporal_d10 = as.numeric(quantile(temporal_hours, 0.10, na.rm = TRUE))
  )
}

episode_size <- 24
set_st_excess <- get_spatiotemp_excess(comephore, quantile = q, 
                     remove_zeros = TRUE)
res <- lapply(dmins, function(dm) {
  sel <- get_s0t0_pairs(
  sites_coords = grid_coords_km,
  data = comephore,
  min_spatial_dist = dm,
  episode_size = episode_size,
  set_st_excess = set_st_excess,
  n_max_episodes = 10000,
  latlon = FALSE,
  beta = 0
  )

  st <- episode_pair_stats(
  sel = sel,
  sites_coords = grid_coords_km,
  dmin_km = dm,
  delta_steps = episode_size,
  latlon = FALSE,
  beta = 0
  )

  data.frame(
  dmin_km = dm,
  n_episodes = st$n_episodes,
  n_pairs_dtlt_delta = st$n_pairs_dtlt_delta,
  p_close_space_50km_all = st$p_close_space_50km,
  p_close_space_50km_given_dtlt_delta = st$p_close_space_50km_given_dtlt_delta,
  p_conflict_dmin_delta = st$p_conflict_dmin_delta,
  spatial_d10_km = st$spatial_d10,
  temporal_d10_hours = st$temporal_d10
  )
})

df_tradeoff <- bind_rows(res)

pA <- ggplot(df_tradeoff, aes(x = dmin_km, y = n_episodes)) +
  geom_line(size = 1.1, color = btfgreen) +
  geom_point(size = 2, color = btfgreen) +
  geom_vline(xintercept = 5, linetype = "dashed", color = "red") +
  theme_minimal() +
  labs(
  x = expression(d[min]~"(km)"),
  y = "Number of selected episodes"
  ) 

foldername <- paste0(im_folder, "/optim/comephore/choice_config/")
if (!dir.exists(foldername)) {
  dir.create(foldername, recursive = TRUE)
}
filename <- paste0(foldername, "tradeoff_dmin_episodes.png")
ggsave(filename, plot = pA, width = 7, height = 5, units = "in", dpi = 300)

dmin_fixed <- 5
delta_grid <- c(10, 12, 15, 20, 24, 30, 36, 38, 40, 45, 48)  # in "steps" of 5 minutes

set_st_excess <- get_spatiotemp_excess(comephore, quantile = q,
                     remove_zeros = TRUE)

res_delta <- lapply(delta_grid, function(delta_steps) {

  sel <- get_s0t0_pairs(
  sites_coords = grid_coords_km,
  data = comephore,
  min_spatial_dist = dmin_fixed,
  episode_size = delta_steps,
  set_st_excess = set_st_excess,
  n_max_episodes = 10000,
  latlon = FALSE,
  beta = 0
  )

  st <- episode_pair_stats(
  sel = sel,
  sites_coords = grid_coords_km,
  dmin_km = dmin_fixed,
  delta_steps = delta_steps,
  latlon = FALSE,
  beta = 0
  )

  data.frame(
  delta_steps = delta_steps,
  delta_hours = delta_steps,
  n_episodes = st$n_episodes,
  n_pairs_dtlt_delta = st$n_pairs_dtlt_delta,
  p_close_space_50km_given_dtlt_delta = st$p_close_space_50km_given_dtlt_delta,
  p_conflict_dmin_delta = st$p_conflict_dmin_delta,
  spatial_d10_km = st$spatial_d10,
  temporal_d10_hours = st$temporal_d10
  )
})

df_tradeoff_delta <- bind_rows(res_delta)


p_deltaA <- ggplot(df_tradeoff_delta, aes(x = delta_steps, y = n_episodes)) +
  geom_line(size = 1.1, color = btfgreen) +
  geom_point(size = 2, color = btfgreen) +
  geom_vline(xintercept = 24, linetype = "dashed", color = "red") +
  theme_minimal() +
  labs(
  x = expression(delta~"(hours)"),
  y = "Number of selected episodes"
  )

print(p_deltaA)

ggsave(paste0(foldername, "tradeoff_delta_episodes_dmin5.png"),
     plot = p_deltaA, width = 7, height = 5, units = "in", dpi = 300)


















min_spatial_dist <- 5
s0t0_set <- get_s0t0_pairs(grid_coords_km, comephore,
                            min_spatial_dist = min_spatial_dist,
                            episode_size = delta,
                            set_st_excess = set_st_excess,
                            n_max_episodes = 10000,
                            latlon = FALSE,
                            beta = 0)

selected_points <- s0t0_set
n_episodes <- length(selected_points$s0)
list_episodes_points <- get_extreme_episodes(selected_points, comephore,
                                     episode_size = delta, unif = FALSE)

list_episodes <- list_episodes_points$episodes
s0_list <- selected_points$s0
t0_list <- selected_points$t0
u_list <- selected_points$u_s0

# Compute lags and excesses #####################################################
df_coords <- as.data.frame(grid_coords_km)
tau_vect <- 0:10

list_results <- parallel::mclapply(1:length(s0_list), function(i) {
  s0 <- s0_list[i]
  s0_coords <- df_coords[s0, ]
  u <- u_list[i]
  episode <- list_episodes[[i]]
  ind_t0_ep <- 0
  lags <- get_conditional_lag_vectors(df_coords, s0_coords, ind_t0_ep,
                                tau_vect, latlon = FALSE)
  excesses <- empirical_excesses_rpar(episode, threshold = u,
                                  df_lags = lags, t0 = ind_t0_ep)
  list(lags = lags, excesses = excesses)
}, mc.cores = ncores)

list_lags <- lapply(list_results, `[[`, "lags")
list_excesses <- lapply(list_results, `[[`, "excesses")

# ADD ADVECTION ESTIMATES ######################################################
adv_filename <- paste0(
  data_folder, "comephore/adv_estim/timeframe_t0_plus_minus_2/advection_results_q",
  q * 100, "_delta", delta, "_dmin", min_spatial_dist,
  starting_year_suffix, ".csv"
)
adv_df <- read.csv(adv_filename, sep = ",")
adv_df$t0 <- as.POSIXct(adv_df$t0, tz = "UTC")
selected_points$t0_date <- as.POSIXct(selected_points$t0_date, tz = "UTC")

matching_indices <- match(selected_points$t0_date, adv_df$t0)
matching_indices <- matching_indices[!is.na(matching_indices)]
adv_df <- adv_df[matching_indices, ]
rownames(adv_df) <- NULL

selected_episodes <- selected_points
selected_episodes$adv_x <- NA
selected_episodes$adv_y <- NA

for (i in 1:nrow(selected_episodes)) {
  t0_date <- selected_episodes$t0_date[i]
  adv_row <- adv_df[adv_df$t0 == t0_date, ]
  if (nrow(adv_row) > 0) {
    selected_episodes$adv_x[i] <- adv_row$mean_dx_kmh[1]
    selected_episodes$adv_y[i] <- adv_row$mean_dy_kmh[1]
  }
}



# COMPUTE ADVECTION CLASSES ####################################################
selected_episodes$adv_speed <- sqrt(selected_episodes$adv_x^2 +
                                    selected_episodes$adv_y^2)

selected_episodes$adv_direction <- atan2(selected_episodes$adv_y,
                                          selected_episodes$adv_x) * (180 / pi)
selected_episodes$adv_direction[selected_episodes$adv_direction < 0] <-
  selected_episodes$adv_direction[selected_episodes$adv_direction < 0] + 360

# max_speed_bins <- 5
# n_speed_bins <- min(max_speed_bins,
#                     max(1, floor(n_advective / min_per_class)))

# quants_speed <- quantile(selected_episodes$adv_speed,
#                               probs = c(seq(0, 1, length.out = n_speed_bins + 1), 0.95))
speed_class <- as.character(cut(
  selected_episodes$adv_speed,
  breaks = c(0, 1, Inf),
  include.lowest = TRUE,
  labels = c("still", "significant")
))
speed_class <- factor(speed_class, levels = c("still", "significant"))
table(speed_class)

adv_group <- as.character(speed_class)
table(adv_group)
selected_episodes$adv_group <- adv_group
selected_episodes$speed_class <- speed_class
# quick diagnostic per class
adv_summary <- selected_episodes %>%
  group_by(adv_group) %>%
  summarise(
    n = n(),
    median_speed = median(adv_speed, na.rm = TRUE)
  ) %>%
  arrange(desc(n))
print(adv_summary)

# plot of histogram of advection speeds
ggplot(selected_episodes, aes(x = adv_speed)) +
  geom_histogram(binwidth = 0.2, fill = btfgreen, color = "black", alpha = 0.6) +
  labs(title = "", x = "Advection speed (km/h)", y = "Episodes count") +
  theme_minimal() +
  btf_theme

angle_adv <- atan2(selected_episodes$adv_y,
                     selected_episodes$adv_x) %% (2 * pi)
angle_math_deg <- (atan2(selected_episodes$adv_y, selected_episodes$adv_x) * 180/pi) %% 360
angle_geo_deg  <- (90 - angle_math_deg) %% 360

head(angle_adv)
head(angle_math_deg)
head(angle_geo_deg)
head(direction_class)
direction_class <- cut(
  angle_geo_deg,
  breaks = c(-45, 45, 135, 225, 315, 405),
  labels = c("N", "E", "S", "W", "N"),
  include.lowest = TRUE,
  right = FALSE
)


selected_episodes$direction_class <- direction_class

# Plot advection classes
plot_dir_df <- selected_episodes
ggplot(plot_dir_df, aes(x = direction_class, fill = speed_class)) +
  geom_bar(stat = "count", width = 1, color = "black") +
  coord_polar(start = -45 * (pi/180)) +
  scale_fill_brewer(palette = "YlGnBu", name = "Speed class") +
  labs(title = "", x = "Direction class", y = "Episodes count (advective only)") +
  theme_minimal() +
  btf_theme

ggsave(filename = paste0(im_folder, "/optim/comephore/adv_classes_",
  q*100, "q", min_spatial_dist, "dmin", delta, "delta_NSEW.png"))

ggplot(plot_dir_df, aes(x = speed_class)) +
  geom_bar(stat = "count", width = 0.5, fill = btfgreen, alpha = 0.6) +
  scale_x_discrete(drop = TRUE) +
  labs(title = "", x = "Speed class (km/h)", y = "Episodes count") +
  theme_minimal() +
  btf_theme

ggsave(filename = paste0(im_folder, "/optim/comephore/adv_classes_",
  q*100, "q", min_spatial_dist, "dmin", delta, "delta_speedonly.png"))

# OPTIMIZATION FOR EACH ADVECTION CLASS #######################################
hmax <- 10
adv_groups <- unique(selected_episodes$speed_class)
# adv_groups <- adv_groups[!is.na(adv_groups) &
#                          adv_groups != "significant"]  # exclude significant for now
results_all_classes <- list()
min_episodes_per_group <- 30

for (adv_group_selected in adv_groups) {

  selected_indices <- which(selected_episodes$adv_group == adv_group_selected)
  number_episodes  <- length(selected_indices)

  if (number_episodes < min_episodes_per_group) {
    cat("Skipping group", adv_group_selected,
        "with only", number_episodes,
        "episodes (<", min_episodes_per_group, ")\n")
    next
  }

  cat("Optimizing group:", adv_group_selected,
      "with", number_episodes, "episodes\n")

  list_episodes_filtered  <- list_episodes[selected_indices]
  list_lags_filtered      <- list_lags[selected_indices]
  list_excesses_filtered  <- list_excesses[selected_indices]
  wind_df_filtered        <- selected_episodes[selected_indices, c("adv_x", "adv_y")]
  colnames(wind_df_filtered) <- c("vx", "vy")

  stopifnot(
    length(list_episodes_filtered) == nrow(wind_df_filtered),
    length(list_lags_filtered)     == nrow(wind_df_filtered),
    length(list_excesses_filtered) == nrow(wind_df_filtered)
  )

  init_param <- c(beta1, beta2, alpha1, alpha2, 1, 1)

  result <- optim(
    par   = init_param,
    fn    = neg_ll_composite,
    list_lags      = list_lags_filtered,
    list_episodes  = list_episodes_filtered,
    list_excesses  = list_excesses_filtered,
    hmax           = hmax,
    wind_df        = wind_df_filtered,
    latlon         = FALSE,
    distance       = distance_type,
    method         = "L-BFGS-B",
    lower          = c(1e-08, 1e-08, 1e-08, 1e-08, 1e-08, 1e-08),
    upper          = c(10, 10, 1.999, 1.999, 10, 10),
    control        = list(maxit = 10000, trace = 0)
  )

  if (result$convergence != 0) {
    print("Optimization did not converge for group:", adv_group_selected)
    next
  }

  print(result$par)

  results_all_classes[[adv_group_selected]] <- list(
    adv_group   = adv_group_selected,
    par         = result$par,
    nll         = result$value,
    n_episodes  = number_episodes,
    indices     = selected_indices
  )
}

# compute adv with eta1 and eta2 fixed
eta1_fixed <- results_all_classes[[1]]$par[5]
eta2_fixed <- results_all_classes[[1]]$par[6]

# plot the theorical chi from estimated parameters
params_est <- results_all_classes[[1]]$par
chi_th <- theoretical_chi(
  params = params_est,
  df_lags = list_lags_filtered[[3]],
  distance = distance_type,
  latlon = FALSE
)

# for fixed tau
tau_val <- 2
chi_tau <- chi_th[tau == tau_val, ]
plot(chi_tau$hnorm, chi_tau$chi,
     xlab = "Distance h (km)",
     ylab = "Theoretical chi",
     main = paste0("Theoretical chi at tau = ", tau_val))

library(data.table)

res_list <- vector("list", length(list_lags_filtered))

for (i in seq_along(list_lags_filtered)) {

  lags_i     <- list_lags_filtered[[i]]
  excess_i  <- list_excesses_filtered[[i]]
  adv_x     <- wind_df_filtered$vx[i]
  adv_y     <- wind_df_filtered$vy[i]

  # paramètres estimés (exemple)
  params_i <- results_all_classes[[1]]$par
  eta1_i <- params_i[5]
  eta2_i <- params_i[6]
  adv_norm <- sqrt(adv_x^2 + adv_y^2)
  adv_norm_transformed <- eta1_i * adv_norm^eta2_i
  if (adv_norm > 0) {
    adv_x <- adv_x / adv_norm * adv_norm_transformed
    adv_y <- adv_y / adv_norm * adv_norm_transformed
  } else {
    adv_x <- 0
    adv_y <- 0
  }

  # chi théorique pour cet épisode
  chi_th_i <- theoretical_chi(
    params   = params_i,
    df_lags  = lags_i,
    latlon   = FALSE
  )

  # assemblage
  res_i <- data.table(
    episode  = i,
    s2       = lags_i$s2,
    tau      = lags_i$tau,
    h        = lags_i$hnorm,
    chi_emp  = excess_i$kij,   # 0 ou 1
    chi_theo = chi_th_i$chi,
    adv_x    = adv_x,
    adv_y    = adv_y
  )

  res_list[[i]] <- res_i
}

res <- rbindlist(res_list, fill = TRUE)

res$hbin <- cut(res$h, breaks = seq(0, 10, by = 1), include.lowest = TRUE)


library(dplyr)

res_cmp <- res %>%
  filter(tau >= 0) %>%
  group_by(tau, hbin) %>%
  summarise(
    chi_emp_bar  = mean(chi_emp, na.rm = TRUE),
    chi_theo_bar = mean(chi_theo, na.rm = TRUE),
    n_pairs = n(),
    .groups = "drop"
  )


# plot chi empirical vs theoretical
ggplot(res_cmp, aes(x = chi_theo_bar, y = chi_emp_bar)) +
  geom_point(alpha = 0.6, color = btfgreen) +
  geom_abline(slope = 1, intercept = 0, color = "red", linetype = "dashed") +
  labs(title = "Empirical vs Theoretical Chi",
       x = "Theoretical Chi",
       y = "Empirical Chi") +
  theme_minimal() +
  btf_theme


# fixed tau
tau_fixed <- 0
res_tau <- res_cmp %>% filter(tau == tau_fixed) 

ggplot(res_tau, aes(x = chi_theo_bar, y = chi_emp_bar)) +
  geom_point(alpha = 0.6, color = btfgreen) +
  geom_abline(slope = 1, intercept = 0, color = "red", linetype = "dashed") +
  labs(title = "Empirical vs Theoretical Chi",
       x = "Theoretical Chi",
       y = "Empirical Chi") +
  theme_minimal() +
  btf_theme



adv_df_significant <- selected_episodes[selected_episodes$speed_class == "significant", ]
adv_df_significant$adv_x_transformed <- eta1_fixed * sign(adv_df_significant$adv_x) * abs(adv_df_significant$adv_x)^eta2_fixed
adv_df_significant$adv_y_transformed <- eta1_fixed * sign(adv_df_significant$adv_y) * abs(adv_df_significant$adv_y)^eta2_fixed

plot(adv_df_significant$adv_x,
     adv_df_significant$adv_y,
     xlab = "Transformed adv_x",
     ylab = "Transformed adv_y",
     main = "Transformed Advection Vectors (Significant Speed Class)")

x <- adv_df_significant$adv_x
y <- adv_df_significant$adv_y

r <- sqrt(x^2 + y^2)
theta <- atan2(y, x)

r_new <- eta1_fixed * r^eta2_fixed

adv_df_significant$adv_x_transformed <- r_new * cos(theta)
adv_df_significant$adv_y_transformed <- r_new * sin(theta)

plot(adv_df_significant$adv_x_transformed,
     adv_df_significant$adv_y_transformed,
     xlab = "Transformed adv_x",
     ylab = "Transformed adv_y",
     main = "Transformed Advection Vectors (Significant Speed Class)")




# get a dataframe with all results to save it
if (length(results_all_classes) == 0) {
  stop("No advection class produced a valid optim result.")
}

results_df <- do.call(rbind, lapply(results_all_classes, function(res) {
  data.frame(
    adv_group = res$adv_group,
    beta1 = res$par[1],
    beta2 = res$par[2],
    alpha1 = res$par[3],
    alpha2 = res$par[4],
    eta1 = res$par[5],
    eta2 = res$par[6],
    nll = res$nll,
    n_episodes  = res$n_episodes
  )
}))


# SAVE RESULTS #################################################################
eta_type <- if (!is.na(fixed_eta1) && !is.na(fixed_eta2)) {
  "fixed_eta"
} else if (!is.na(fixed_eta1)) {
  "fixed_eta1"
} else if (!is.na(fixed_eta2)) {
  "fixed_eta2"
} else {
  "free_eta"
}

foldername_res <- file.path(
  paste0(data_folder, "comephore/optim_results"),
  distance_type,
  eta_type
)

if (!dir.exists(foldername_res)) {
  dir.create(foldername_res, recursive = TRUE)
}

# Save optimization results for all classes
filename <- paste0(foldername_res, "/results_com_classes_q",
           q * 100, "_delta", delta, "_dmin", min_spatial_dist, ".csv")
write.csv(results_df, filename, row.names = FALSE)
cat("Results saved to", filename, "\n")
print(results_df)





# JACKKNIFE ANALYSIS ##########################################################

jackknife_all_classes <- list()

get_season <- function(date) {
  yday <- as.integer(format(date, "%j"))
  if (yday >= 80  & yday <= 171) "Spring"
  else if (yday >= 172 & yday <= 263) "Summer"
  else if (yday >= 264 & yday <= 354) "Autumn"
  else "Winter"
}

for (adv_group_selected in names(results_all_classes)) {

  cat("Jackknife for group:", adv_group_selected, "\n")

  mle_obj  <- results_all_classes[[adv_group_selected]]
  mle_par  <- mle_obj$par
  indices  <- mle_obj$indices

  list_episodes_filtered <- list_episodes[indices]
  list_lags_filtered     <- list_lags[indices]
  list_excesses_filtered <- list_excesses[indices]
  wind_df_filtered       <- selected_episodes[indices, c("adv_x", "adv_y")]

  selected_episodes_filtered <- selected_episodes[indices, ]
  selected_episodes_filtered$t0_date <- as.POSIXct(
    selected_episodes_filtered$t0_date,
    format = "%Y-%m-%d %H:%M:%S",
    tz = "UTC"
  )

  season_vec <- sapply(selected_episodes_filtered$t0_date, get_season)
  year_vec   <- format(selected_episodes_filtered$t0_date, "%Y")
  season_year_vec <- paste(season_vec, year_vec)
  unique_season_years <- sort(unique(season_year_vec))

  jack_estimates_list <- parallel::mclapply(
    unique_season_years,
    function(season_year) {

      exclude_idx <- which(season_year_vec == season_year)

      jack_episodes  <- list_episodes_filtered[-exclude_idx]
      jack_lags      <- list_lags_filtered[-exclude_idx]
      jack_excesses  <- list_excesses_filtered[-exclude_idx]
      jack_wind      <- wind_df_filtered[-exclude_idx, , drop = FALSE]
      colnames(jack_wind) <- c("vx", "vy")
      res <- tryCatch({
        optim(
          par   = mle_par,
          fn    = neg_ll_composite,
          list_lags      = jack_lags,
          list_episodes  = jack_episodes,
          list_excesses  = jack_excesses,
          hmax           = hmax,
          wind_df        = jack_wind,
          latlon         = FALSE,
          distance       = distance_type,
          method         = "L-BFGS-B",
          lower          = c(1e-08, 1e-08, 1e-08, 1e-08, 1e-08, 1e-08),
          upper          = c(10, 10, 1.999, 1.999, 10, 10),
          control        = list(maxit = 10000)
        )
      }, error = function(e) NULL)

      if (!is.null(res) && res$convergence == 0) {

        too_close <- any(
          abs(res$par - c(1e-08, 1e-08, 1e-08, 1e-08, 1e-08, 1e-08)) < 1e-6 |
          abs(res$par - c(10, 10, 1.999, 1.999, 10, 10)) < 1e-6
        )

        if (!too_close) return(res$par)
      }

      rep(NA, length(mle_par))
    },
    mc.cores = ncores
  )

  jack_estimates <- do.call(rbind, jack_estimates_list)
  jack_estimates <- na.omit(jack_estimates)

  if (nrow(jack_estimates) == 0) next

  n_eff <- nrow(jack_estimates)
  jack_mean <- colMeans(jack_estimates)

  pseudo_values <- matrix(NA, nrow = n_eff, ncol = length(mle_par))
  for (i in seq_len(n_eff)) {
    pseudo_values[i, ] <- n_eff * jack_mean - (n_eff - 1) * jack_estimates[i, ]
  }

  jack_mean_pseudo <- colMeans(pseudo_values)
  jack_se <- apply(pseudo_values, 2, sd) / sqrt(n_eff)

  z <- qnorm(0.975)

  jk_df <- data.frame(
    adv_group     = adv_group_selected,
    Parameter     = c("beta1", "beta2", "alpha1", "alpha2", "eta1", "eta2"),
    Estimate_full = mle_par,
    Estimate_jk   = jack_mean_pseudo,
    StdError      = jack_se,
    CI_lower      = jack_mean_pseudo - z * jack_se,
    CI_upper      = jack_mean_pseudo + z * jack_se,
    n_eff         = n_eff
  )

  jackknife_all_classes[[adv_group_selected]] <- jk_df
}



# Save jackknife results for all classes
if (length(jackknife_all_classes) > 0) {
  jackknife_combined <- do.call(rbind, jackknife_all_classes)
  filename_jk <- paste0(foldername_res, "/jackknife_results_com_classes_q",
             q * 100, "_delta", delta, "_dmin", min_spatial_dist, ".csv")
  write.csv(jackknife_combined, filename_jk, row.names = FALSE)
  cat("Jackknife results saved to", filename_jk, "\n")
}
