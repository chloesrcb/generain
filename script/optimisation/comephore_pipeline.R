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



# PREPARE EPISODES WITH AVAILABLE ADVECTION ###################################
selected_episodes$adv_speed <- sqrt(selected_episodes$adv_x^2 +
                                    selected_episodes$adv_y^2)
speed_class <- cut(
  selected_episodes$adv_speed,
  breaks = c(0, 5.6),
  labels = c("omsev-like"),
  include.lowest = TRUE
)

table(speed_class)
selected_episodes$adv_group <- speed_class

hmax <- 10
min_episodes_per_group <- 30
adv_groups <- levels(speed_class)
adv_groups <- adv_groups[!is.na(adv_groups)]
# remove 'still' group if speed threshold is applied
adv_groups <- adv_groups[adv_groups != "still"]


if (!length(adv_groups)) {
  stop("No advection classes available (check advection estimates).")
}

library(data.table)
results_all_classes <- list()
res_cmp_classes <- list()

for (adv_group_selected in adv_groups) {
  selected_indices <- which(selected_episodes$adv_group == adv_group_selected)
  number_episodes <- length(selected_indices)

  if (number_episodes < min_episodes_per_group) {
    cat("Skipping group", adv_group_selected, "with only",
        number_episodes, "episodes (<", min_episodes_per_group, ")\n")
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
    warning("Optimization did not converge for group ", adv_group_selected)
    next
  }
  print(result$par)

  res_list <- vector("list", length(list_lags_filtered))

  for (i in seq_along(list_lags_filtered)) {
    lags_i    <- list_lags_filtered[[i]]
    excess_i  <- list_excesses_filtered[[i]]
    adv_x     <- wind_df_filtered$vx[i]
    adv_y     <- wind_df_filtered$vy[i]

    params_i <- result$par
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

    chi_th_i <- theoretical_chi(
      params   = params_i,
      df_lags  = lags_i,
      latlon   = FALSE
    )

    res_list[[i]] <- data.table(
      episode  = i,
      s2       = lags_i$s2,
      tau      = lags_i$tau,
      h        = lags_i$hnorm,
      chi_emp  = excess_i$kij,
      chi_theo = chi_th_i$chi,
      adv_x    = adv_x,
      adv_y    = adv_y
    )
  }

  res <- rbindlist(res_list, fill = TRUE)
  res$hbin <- cut(res$h, breaks = seq(0, 10, by = 1), include.lowest = TRUE)

  res_cmp <- res %>%
    filter(tau >= 0) %>%
    group_by(tau, hbin) %>%
    summarise(
      chi_emp_bar  = mean(chi_emp, na.rm = TRUE),
      chi_theo_bar = mean(chi_theo, na.rm = TRUE),
      n_pairs = n(),
      .groups = "drop"
    )

  res_cmp_classes[[adv_group_selected]] <- res_cmp

  results_all_classes[[adv_group_selected]] <- list(
    adv_group = adv_group_selected,
    par = result$par,
    nll = result$value,
    n_episodes = number_episodes,
    indices = selected_indices,
    res_cmp = res_cmp
  )

  params_rds <- list(
    par = result$par,
    q = q,
    delta = delta,
    min_spatial_dist = min_spatial_dist,
    distance = distance_type,
    adv_group = adv_group_selected,
    speed_threshold = speed_threshold
  )
  params_rds_path <- file.path(
    paste0(data_folder, "comephore/optim_results"),
    distance_type,
    "free_eta",
    sprintf("params_q%d_delta%d_dmin%d_group_%s.rds",
            as.integer(q * 100), delta, min_spatial_dist, adv_group_selected)
  )
  dir.create(dirname(params_rds_path), recursive = TRUE, showWarnings = FALSE)
  saveRDS(params_rds, params_rds_path)
  cat("Saved COMEPHORE parameters for group", adv_group_selected,
      "to", params_rds_path, "\n")
}


if (!length(results_all_classes)) {
  stop("No advection class produced a valid optimization result.")
}

default_group <- names(results_all_classes)[1]
cat("Using group", default_group, "for diagnostic plots.\n")

ref_cmp <- res_cmp_classes[[default_group]]
tau_fixed <- 0
ref_tau <- ref_cmp %>% filter(tau == tau_fixed)

ggplot(ref_cmp, aes(x = chi_theo_bar, y = chi_emp_bar, size = n_pairs)) +
  geom_point(alpha = 0.6, color = btfgreen) +
  geom_abline(slope = 1, intercept = 0, color = "red", linetype = "dashed") +
  labs(
    x = "Theoretical Chi",
    y = "Empirical Chi",
    size = "Number of pairs"
  ) +
  theme_minimal() +
  btf_theme

# save plot
filename <- paste0(
  data_folder, "comephore/optim_results/",
  distance_type, "/free_eta/",
  "chi_emp_vs_chi_theo_q", q * 100,
  "_delta", delta,
  "_dmin", min_spatial_dist,
  "_group_", default_group,
  ".png"
)
ggsave(filename, width = 8, height = 7)
cat("Plot saved to", filename, "\n")


ggplot(ref_tau, aes(x = chi_theo_bar, y = chi_emp_bar)) +
  geom_point(alpha = 0.6, color = btfgreen) +
  geom_abline(slope = 1, intercept = 0, color = "red", linetype = "dashed") +
  labs(
    title = paste("Empirical vs Theoretical Chi (tau = ", tau_fixed,
                  ", ", default_group, ")", sep = ""),
    x = "Theoretical Chi",
    y = "Empirical Chi"
  ) +
  theme_minimal() +
  btf_theme

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
    n_episodes = res$n_episodes
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
