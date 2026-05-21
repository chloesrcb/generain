# remove objects from the environment
rm(list = ls())

# clear console
cat("\014")

muse <- F

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
filename_com <- paste0(data_folder, "comephore/rebuild_clean/comephore_2008_2025_within5km.csv")
comephore_raw <- read.csv(filename_com, sep = ",")
filename_loc <- paste0(data_folder, "comephore/rebuild_clean/coords_pixels_within5km.csv")
loc_px <- read.csv(filename_loc, sep = ",")

loc_px <- loc_px[loc_px$pixel_name %in% colnames(comephore_raw)[-1], ]
rownames(loc_px) <- NULL

df_comephore <- as.data.frame(comephore_raw)
df_comephore$date <- as.POSIXct(df_comephore$date,
                                format = "%Y-%m-%d %H:%M:%S", tz = "UTC")
rownames(df_comephore) <- format(df_comephore$date, "%Y-%m-%d %H:%M:%S")
tail(df_comephore)
head(df_comephore)
# remove every dates before 2008
df_comephore <- df_comephore[format(df_comephore$date, "%Y") >= "2008", ]

# plot 2025
# df_2025 <- df_comephore[format(df_comephore$date, "%Y") == "2025", ]
# df_2025_p101 <- df_2025$p101
# plot(df_2025_p101, type = "l", main = "p101 in 2025", xlab = "Time step", ylab = "Value")
# remove 2025
# df_comephore <- df_comephore[format(df_comephore$date, "%Y") != "2025", ]
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


# distance matrix in km
dist_mat <- get_dist_mat(grid_coords_km, latlon = FALSE)
# Spatial chi WLSE #############################################################
beta1 <- 0.01611951
beta2 <- 0.7944467
alpha1 <- 1.42891232
alpha2 <- 0.8049938

head(comephore)
# remove column "dates"
comephore <- comephore[, which(colnames(comephore) != "date")]
# CHOOSE EXTREME EPISODE FOR R-PARETO ##########################################
set_st_excess <- get_spatiotemp_excess(data = comephore, quantile = q,
                                      remove_zeros = TRUE)
# first_ts <- as.POSIXct(rownames(comephore), tz = "UTC")
# starting_year <- year(first_ts)[1]
# starting_year_suffix <- if (starting_year == 2008) "" else paste0("_from", starting_year)

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


######################################################################
# OPTIM
min_spatial_dist <- 5
delta <- 24
s0t0_set <- get_s0t0_pairs(grid_coords_km, comephore,
                            min_spatial_dist = min_spatial_dist,
                            episode_size = delta,
                            set_st_excess = set_st_excess,
                            n_max_episodes = 10000,
                            latlon = FALSE,
                            beta = 0)

selected_points <- s0t0_set
n_episodes <- length(selected_points$s0)
n_episodes
tail(selected_points)
list_episodes_points <- get_extreme_episodes(selected_points, comephore,
                                     episode_size = delta, unif = FALSE)

list_episodes <- list_episodes_points$episodes
episodes <- list_episodes[[130]]
plot(episodes$p101, type = "l", main = "Example episode", xlab = "Time step", ylab = "Value")

s0_list <- selected_points$s0
t0_list <- selected_points$t0
u_list <- selected_points$u_s0

selected_points$t0_date <- as.POSIXct(selected_points$t0_date, format = "%Y-%m-%d %H:%M:%S", tz = "UTC")

# Round to the next hour
selected_points$t0_date_rounded <- ceiling_date(selected_points$t0_date, "hour")

# Round to the next hour
datetimes <- unique(selected_points$t0_date)
datetimes_hour <- unique(selected_points$t0_date_rounded)

# save datetime list to csv
datetime_filename <- paste(data_folder, "/comephore/episodes/t0_episodes/t0_episodes_q", q * 100,
                           "_delta", delta, "_dmin", min_spatial_dist,
                           ".csv", sep = "")
write.csv(data.frame(t0_date = datetimes), datetime_filename, row.names = FALSE)


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
  data_folder, "comephore/adv_estim/advection_results_q",
  q * 100, "_delta", delta, "_dmin", min_spatial_dist, ".csv"
)
adv_df <- read.csv(adv_filename, sep = ",")
adv_df$t0 <- as.POSIXct(adv_df$t0, tz = "UTC")
speed <- sqrt(adv_df$mean_dx_kmh^2 + adv_df$mean_dy_kmh^2)
plot(speed, main = "Advection speed (km/h)", xlab = "Episode index", ylab = "Speed (km/h)")
hist(speed, main = "Histogram of advection speeds", xlab = "Speed (km/h)", ylab = "Frequency")

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


selected_episodes$t0_date <- as.POSIXct(
  selected_episodes$t0_date,
  format = "%Y-%m-%d %H:%M:%S",
  tz = "UTC"
)

# same_temporality <- function(t0_i, t0_j, delta, step_min = 5) {
#   abs(as.numeric(difftime(t0_i, t0_j, units = "mins"))) < delta * step_min
# }

# t0_vec <- selected_episodes$t0_date
# n_ep <- length(t0_vec)

# temporal_overlap <- matrix(FALSE, n_ep, n_ep)

# for (i in 1:n_ep) {
#   for (j in 1:n_ep) {
#     if (i != j) {
#       temporal_overlap[i, j] <- same_temporality(t0_vec[i], t0_vec[j], delta)
#     }
#   }
# }

# has_temporal_overlap <- apply(temporal_overlap, 1, any)
# table(has_temporal_overlap)

# choose local radius for episodes with temporal overlap
# IMPORTANT: adapt depending on unit of hnorm
# d_local <- 2.5  # if hnorm is in km
# # d_local <- 600 # if hnorm is in meters

# list_lags_overlap <- vector("list", length(list_lags))
# list_excesses_overlap <- vector("list", length(list_excesses))

# for (i in seq_along(list_lags)) {
#   if (has_temporal_overlap[i]) {
#     keep <- list_lags[[i]]$hnorm <= d_local
#   } else {
#     keep <- rep(TRUE, nrow(list_lags[[i]]))
#   }

#   list_lags_overlap[[i]] <- list_lags[[i]][keep, , drop = FALSE]
#   list_excesses_overlap[[i]] <- list_excesses[[i]][keep, , drop = FALSE]
# }

# # selected_episodes_filtered <- selected_episodes
# # list_episodes_filtered <- list_episodes
# list_lags <- list_lags_overlap
# list_excesses <- list_excesses_overlap


# # remove 2025 episodes
# idx_2025 <- which(format(selected_episodes$t0_date, "%Y") == "2025")
# if (length(idx_2025) > 0) {
#   selected_episodes <- selected_episodes[-idx_2025, ]
#   list_episodes <- list_episodes[-idx_2025]
#   list_lags <- list_lags[-idx_2025]
#   list_excesses <- list_excesses[-idx_2025]
#   cat("Removed", length(idx_2025), "episodes from 2025\n")
# } else {
#   cat("No episodes from 2025 found\n")
# }

# # keep only episodes with speed adv <=5.6 km/h
# index_speed_above <- which(sqrt(selected_episodes$adv_x^2 + selected_episodes$adv_y^2) > 5.6)
# length(index_speed_above)
# index_speed <- which(sqrt(selected_episodes$adv_x^2 + selected_episodes$adv_y^2) <= 5.6)
# selected_episodes <- selected_episodes[index_speed, ]
# list_episodes <- list_episodes[index_speed]
# list_lags <- list_lags[index_speed]
# list_excesses <- list_excesses[index_speed]


# COMPUTE ADVECTION CLASSES ####################################################
selected_episodes$adv_speed <- sqrt(selected_episodes$adv_x^2 +
                                    selected_episodes$adv_y^2)

selected_episodes$adv_direction <- atan2(selected_episodes$adv_y,
                                          selected_episodes$adv_x) * (180 / pi)
selected_episodes$adv_direction[selected_episodes$adv_direction < 0] <-
  selected_episodes$adv_direction[selected_episodes$adv_direction < 0] + 360

V_adv <- selected_episodes[, c("adv_x", "adv_y")]
colnames(V_adv) <- c("vx", "vy")
# OPTIMIZATION ON ALL EPISODES ################################################################
init <- c(beta1, beta2, alpha1, alpha2, 1, 1)

hmax <- 10
result <- optim(
  par = init,
  fn = neg_ll_composite,
  list_lags = list_lags,
  list_episodes = list_episodes,
  list_excesses = list_excesses,
  hmax = hmax,
  wind_df = V_adv,
  latlon = FALSE,
  distance = "euclidean",
  method = "L-BFGS-B",
  lower = c(1e-08, 1e-08, 1e-08, 1e-08, 1e-08, 1e-08),
  upper = c(10, 10, 1.999, 1.999, 10, 10),
  control = list(maxit = 20000, trace = 1)
)

if(result$convergence == 0) {
  cat("Optimization converged successfully\n")
} else {
  warning("Optimization did not converge")
}

result

results_df <- data.frame(
    beta1 = result$par[1],
    beta2 = result$par[2],
    alpha1 = result$par[3],
    alpha2 = result$par[4],
    eta1 = result$par[5],
    eta2 = result$par[6],
    nll = result$value
  )


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
results_df <- read.csv(filename, sep = ",")

# JACKKNIFE ANALYSIS ##########################################################

selected_episodes$t0_date <- as.POSIXct(
  selected_episodes$t0_date,
  format = "%Y-%m-%d %H:%M:%S",
  tz = "UTC"
)

month_vec <- format(selected_episodes$t0_date, "%m")
year_vec  <- format(selected_episodes$t0_date, "%Y")
month_year_vec <- paste(year_vec, month_vec, sep = "-")
unique_month_years <- sort(unique(month_year_vec))

param_names <- c("beta1", "beta2", "alpha1", "alpha2", "eta1", "eta2")
npar <- length(param_names)

lower_bounds <- c(1e-08, 1e-08, 1e-08, 1e-08, 1e-08, 1e-08)
upper_bounds <- c(10, 10, 1.999, 1.999, 10, 10)

theta_full <- as.numeric(result$par)
names(theta_full) <- param_names

# folder for jackknife outputs
jk_folder <- file.path(foldername_res, "jackknife_estimates")
if (!dir.exists(jk_folder)) {
  dir.create(jk_folder, recursive = TRUE)
}

jackknife_one_block <- function(month_year) {
  exclude_idx <- which(month_year_vec == month_year)

  if (length(exclude_idx) == 0) {
    return(c(rep(NA_real_, npar), convergence = NA_real_, n_removed = 0))
  }

  jack_episodes <- list_episodes[-exclude_idx]
  jack_lags <- list_lags[-exclude_idx]
  jack_excesses <- list_excesses[-exclude_idx]
  jack_wind <- V_adv[-exclude_idx, , drop = FALSE]

  if (length(jack_episodes) == 0 ||
      length(jack_lags) == 0 ||
      length(jack_excesses) == 0 ||
      nrow(jack_wind) == 0) {
    return(c(rep(NA_real_, npar), convergence = NA_real_, n_removed = length(exclude_idx)))
  }

  res <- tryCatch(
    optim(
      par = theta_full,
      fn = neg_ll_composite,
      list_lags = jack_lags,
      list_episodes = jack_episodes,
      list_excesses = jack_excesses,
      hmax = hmax,
      wind_df = jack_wind,
      latlon = FALSE,
      distance = "euclidean",
      method = "L-BFGS-B",
      lower = lower_bounds,
      upper = upper_bounds,
      control = list(maxit = 10000)
    ),
    error = function(e) NULL
  )

  if (is.null(res)) {
    return(c(rep(NA_real_, npar), convergence = NA_real_, n_removed = length(exclude_idx)))
  }

  if (is.null(res$par) || !is.numeric(res$par) || length(res$par) != npar) {
    return(c(rep(NA_real_, npar), convergence = NA_real_, n_removed = length(exclude_idx)))
  }

  out <- c(as.numeric(res$par),
           convergence = ifelse(is.null(res$convergence), NA_real_, res$convergence),
           n_removed = length(exclude_idx))
  names(out)[1:npar] <- param_names
  return(out)
}

jack_estimates_list <- parallel::mclapply(
  unique_month_years,
  jackknife_one_block,
  mc.cores = ncores
)

cat("Jackknife estimates computed for all month-year blocks\n")

jack_raw <- do.call(rbind, jack_estimates_list)
jack_raw <- as.data.frame(jack_raw)
jack_raw$block <- unique_month_years

# reorder columns
jack_raw <- jack_raw[, c("block", param_names, "convergence", "n_removed")]

filename_jack_all <- file.path(
  jk_folder,
  paste0("all_jk_q", q * 100, "_delta", delta, "_dmin", min_spatial_dist, ".csv")
)
write.csv(jack_raw, filename_jack_all, row.names = FALSE)
cat("Raw jackknife estimates saved to", filename_jack_all, "\n")

# Keep only valid converged rows
valid_rows <- complete.cases(jack_raw[, param_names]) &
  !is.na(jack_raw$convergence) &
  jack_raw$convergence == 0

jack_valid <- jack_raw[valid_rows, , drop = FALSE]

cat("Number of valid jackknife replicates:", nrow(jack_valid), "\n")
cat("Number of discarded replicates:", nrow(jack_raw) - nrow(jack_valid), "\n")

if (nrow(jack_valid) < 2) {
  cat("Not enough valid jackknife replicates to compute standard errors.\n")
} else {

  jack_estimates <- as.matrix(jack_valid[, param_names, drop = FALSE])
  storage.mode(jack_estimates) <- "double"

  G <- nrow(jack_estimates)
  theta_dot <- colMeans(jack_estimates)

  # Jackknife bias-corrected estimate
  theta_jack <- G * theta_full - (G - 1) * theta_dot

  # Equivalent pseudo-values (optional, for checking / export)
  pseudo_values <- matrix(
    NA_real_,
    nrow = G,
    ncol = npar,
    dimnames = list(jack_valid$block, param_names)
  )

  for (i in seq_len(G)) {
    pseudo_values[i, ] <- G * theta_full - (G - 1) * jack_estimates[i, ]
  }

  # Jackknife standard error from leave-one-block-out estimates
  jack_se <- sqrt(
    (G - 1) / G * colSums(
      (jack_estimates - matrix(theta_dot, G, npar, byrow = TRUE))^2
    )
  )

  z <- qnorm(0.975)

  jk_df <- data.frame(
    Parameter = param_names,
    Estimate_full = as.numeric(theta_full),
    Estimate_jackknife = as.numeric(theta_jack),
    StdError = as.numeric(jack_se),
    CI_lower = as.numeric(theta_jack - z * jack_se),
    CI_upper = as.numeric(theta_jack + z * jack_se),
    n_eff = G
  )

  filename_jk <- file.path(
    foldername_res,
    paste0("jackknife_results_com_q", q * 100,
           "_delta", delta, "_dmin", min_spatial_dist, ".csv")
  )
  write.csv(jk_df, filename_jk, row.names = FALSE)
  cat("Jackknife results saved to", filename_jk, "\n")
  print(jk_df)

  # Save pseudo-values
  pseudo_df <- as.data.frame(pseudo_values)
  pseudo_df$block <- jack_valid$block
  pseudo_df <- pseudo_df[, c("block", param_names)]

  filename_pseudo <- file.path(
    jk_folder,
    paste0("pseudo_values_q", q * 100,
           "_delta", delta, "_dmin", min_spatial_dist, ".csv")
  )
  write.csv(pseudo_df, filename_pseudo, row.names = FALSE)
  cat("Pseudo-values saved to", filename_pseudo, "\n")


  # --------------------------------------------------------------------------
  # OPTIONAL: log-scale jackknife for eta1 and eta2
  # --------------------------------------------------------------------------
  eta_names <- c("eta1", "eta2")

  if (all(jack_estimates[, eta_names] > 0) && all(theta_full[eta_names] > 0)) {

    jack_eta_log <- log(jack_estimates[, eta_names, drop = FALSE])
    theta_full_log <- log(theta_full[eta_names])
    theta_dot_log <- colMeans(jack_eta_log)

    theta_jack_log <- G * theta_full_log - (G - 1) * theta_dot_log

    jack_se_log <- sqrt(
      (G - 1) / G * colSums(
        (jack_eta_log - matrix(theta_dot_log, G, length(eta_names), byrow = TRUE))^2
      )
    )

    ci_log_lower <- theta_jack_log - z * jack_se_log
    ci_log_upper <- theta_jack_log + z * jack_se_log

    jk_eta_log_df <- data.frame(
      Parameter = eta_names,
      Estimate_full = as.numeric(theta_full[eta_names]),
      Estimate_jackknife = as.numeric(exp(theta_jack_log)),
      StdError_logscale = as.numeric(jack_se_log),
      CI_lower = as.numeric(exp(ci_log_lower)),
      CI_upper = as.numeric(exp(ci_log_upper)),
      n_eff = G
    )

    filename_jk_eta_log <- file.path(
      foldername_res,
      paste0("jackknife_results_log_eta_q", q * 100,
             "_delta", delta, "_dmin", min_spatial_dist, ".csv")
    )
    write.csv(jk_eta_log_df, filename_jk_eta_log, row.names = FALSE)
    cat("Log-scale jackknife results for eta saved to", filename_jk_eta_log, "\n")
    print(jk_eta_log_df)

  } else {
    cat("Log-scale CI for eta not computed because some eta estimates are non-positive.\n")
  }
}

# com_results <- c(0.1579947, 0.9923270, 0.5800128, 0.6627969, 3.8962660, 2.2208320)
# com_params <- com_results
com_params <- results_df[, c("beta1", "beta2", "alpha1", "alpha2", "eta1", "eta2")]
com_params <- as.numeric(com_params[1, ])
# # get distance matrix and breaks for chi estimation
df_dist <- reshape_distances(dist_mat)
n_hbins <- 30
h_all <- df_dist$value

h_breaks_com <- quantile(
  h_all,
  probs = seq(0, 1, length.out = n_hbins + 1),
  na.rm = TRUE
)
# h_breaks_com <- seq(0, 13, by=2)

table(cut(h_all, breaks = h_breaks_com, include.lowest = TRUE), useNA = "ifany")
length(h_breaks_com)
h_breaks_com

h_breaks_com <- unique(as.numeric(h_breaks_com))
h_breaks_com[length(h_breaks_com)] <- 13

chi_com <- compute_group_chi(
  list_lags = list_lags,
  list_excesses = list_excesses,
  wind_df = V_adv,
  params = com_params,
  tau_fixed = 0,
  h_breaks = h_breaks_com,
  adv_transform = TRUE
)

corr_com <- summary_correlation(chi_com$res_cmp)
message(sprintf("COM chi Spearman correlation: %.3f", corr_com))

# plot chi_come emp vs theo
print(chi_com$plots$all)

ggplot(res_com_cmp, aes(x = chi_theo_bar, y = chi_emp_bar, color = factor(tau), size = n_pairs)) +
  geom_point(alpha = 0.8) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
  theme_bw()

res_com_cmp <- chi_com$res_cmp
# ggplot(res_com_cmp, aes(x = chi_theo_bar, y = chi_emp_bar, size = n_pairs)) +
#   geom_point(alpha = 0.7, color = btfgreen) +
#   geom_abline(linetype = "dashed", color = "red") +
#   theme_bw() +
#   scale_size_continuous(range = c(1, 4)) +
#   labs(x = "Theoretical Chi", y = "Empirical Chi",
#       size = "Number of pairs") +
#   btf_theme

ggplot(res_com_cmp, aes(x = chi_theo_bar, y = chi_emp_bar, size = n_pairs)) +
  geom_point(alpha = 0.7, color = btfgreen) +
  geom_abline(linetype = "dashed", color = "red") +
  theme_bw() +
  scale_size_continuous(range = c(1, 4)) +
  labs(x = TeX("Theoretical $\\chi$"), y = TeX("Empirical $\\chi$"),
      size = "Number of pairs") +
  btf_theme +
  theme(legend.position = "right",
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 14))
# save plot
foldername <- paste0(im_folder, "workflows/full_pipeline/")
if (!dir.exists(foldername)) {
  dir.create(foldername, recursive = TRUE)
}
filename <- paste0(foldername, "com_chi_th_emp.png")
ggsave(filename,
       width = 7,
       height = 6)

