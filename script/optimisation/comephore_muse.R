# remove objects from the environment
rm(list = ls())

# clear console
cat("\014")

muse <- FALSE

if (muse) {
  folder_muse <- "/home/serrec/work_rainstsimu/rpareto/"
  setwd(folder_muse)
  .libPaths("~/R/x86_64-pc-linux-gnu-library/4.3")  
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
# beta2 <- 0.7944467
beta2 <- 1.553
alpha1 <- 1.42891232
# alpha2 <- 0.8049938
alpha2 <- 0.283

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
# init <- c(0.1579947, 0.9923211, 0.5800001, 0.6627909, 3.896299, 2.221052)
# v0 <- median(selected_episodes$adv_speed , na.rm = TRUE) # 0.8987283

hmax <- 10
result <- optim(
  par = init,
  fn = neg_ll_composite,
  list_lags = list_lags,
  list_episodes = list_episodes,
  list_excesses = list_excesses,
  hmax = hmax,
  wind_df = V_adv,
  v0 = 1,
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

year_vec <- format(selected_episodes$t0_date, "%Y")
unique_years <- sort(unique(year_vec))

param_names <- c("beta1", "beta2", "alpha1", "alpha2", "eta1", "eta2")
npar <- length(param_names)

lower_bounds <- c(1e-08, 1e-08, 1e-08, 1e-08, 1e-08, 1e-08)
upper_bounds <- c(10, 10, 1.999, 1.999, 10, 10)

theta_full <- as.numeric(result$par)
names(theta_full) <- param_names

jk_folder <- file.path(foldername_res, "jackknife_annual_estimates_new")
if (!dir.exists(jk_folder)) {
  dir.create(jk_folder, recursive = TRUE)
}

jackknife_one_year <- function(target_year) {
  cat("Exclusion de l'année :", target_year, "\n")
  exclude_idx <- which(year_vec == target_year)

  if (length(exclude_idx) == 0) {
    return(c(rep(NA_real_, npar), convergence = NA_real_, n_removed = 0))
  }

  jack_episodes <- list_episodes[-exclude_idx]
  jack_lags     <- list_lags[-exclude_idx]
  jack_excesses <- list_excesses[-exclude_idx]
  jack_wind     <- V_adv[-exclude_idx, , drop = FALSE]

  if (length(jack_episodes) == 0 || length(jack_lags) == 0 || 
      length(jack_excesses) == 0 || nrow(jack_wind) == 0) {
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
      v0 = v0,
      latlon = FALSE,
      distance = "euclidean",
      method = "L-BFGS-B",
      lower = lower_bounds,
      upper = upper_bounds,
      control = list(maxit = 15000)
    ),
    error = function(e) NULL
  )

  if (is.null(res) || is.null(res$par) || length(res$par) != npar) {
    return(c(rep(NA_real_, npar), convergence = NA_real_, n_removed = length(exclude_idx)))
  }

  out <- c(as.numeric(res$par),
           convergence = ifelse(is.null(res$convergence), NA_real_, res$convergence),
           n_removed = length(exclude_idx))
  names(out)[1:npar] <- param_names
  return(out)
}

jack_estimates_list <- parallel::mclapply(
  unique_years,
  jackknife_one_year,
  mc.cores = ncores
)

cat("Jackknife estimates computed for all annual blocks\n")

jack_raw <- do.call(rbind, jack_estimates_list)
jack_raw <- as.data.frame(jack_raw)
jack_raw$block <- unique_years
jack_raw <- jack_raw[, c("block", param_names, "convergence", "n_removed")]

filename_jack_all <- file.path(
  jk_folder,
  paste0("all_annual_jk_q", q * 100, "_delta", delta, "_dmin", min_spatial_dist, ".csv")
)
write.csv(jack_raw, filename_jack_all, row.names = FALSE)

valid_rows <- complete.cases(jack_raw[, param_names]) & 
              !is.na(jack_raw$convergence) & 
              jack_raw$convergence == 0

jack_valid <- jack_raw[valid_rows, , drop = FALSE]

G <- nrow(jack_valid)

jack_estimates <- as.matrix(jack_valid[, param_names, drop = FALSE])
storage.mode(jack_estimates) <- "double"


theta_dot <- colMeans(jack_estimates)

pseudo_values <- matrix(NA_real_, nrow = G, ncol = npar, dimnames = list(jack_valid$block, param_names))
for (i in seq_len(G)) {
  pseudo_values[i, ] <- G * theta_full - (G - 1) * jack_estimates[i, ]
}

jack_se <- sqrt(
  ((G - 1) / G) * colSums((jack_estimates - matrix(theta_dot, G, npar, byrow = TRUE))^2)
)

theta_jack_corrected <- colMeans(pseudo_values)

z <- qnorm(0.975)
ci_lower <- theta_full - z * jack_se
ci_upper <- theta_full + z * jack_se

ci_lower <- pmax(1e-08, ci_lower)

jk_df <- data.frame(
  Parameter          = param_names,
  Estimate_full      = as.numeric(theta_full),
  Estimate_jackknife = as.numeric(theta_jack_corrected), 
  StdError           = as.numeric(jack_se),
  CI_lower           = as.numeric(ci_lower),
  CI_upper           = as.numeric(ci_upper),
  n_blocks_eff       = G
)

filename_jk <- file.path(
  foldername_res,
  paste0("jackknife_annual_results_natural_q", q * 100, "_delta", delta, "_dmin", min_spatial_dist, ".csv")
)
write.csv(jk_df, filename_jk, row.names = FALSE)
print(jk_df)
pseudo_values_log[, c("eta1", "eta2")]
pseudo_natural_df <- as.data.frame(pseudo_values)
pseudo_natural_df$block <- jack_valid$block
pseudo_natural_df <- pseudo_natural_df[, c("block", param_names)]

filename_pseudo <- file.path(jk_folder, paste0("pseudo_values_natural_q", q * 100, ".csv"))
write.csv(pseudo_natural_df, filename_pseudo, row.names = FALSE)


influence_df <- jack_valid[, c("block", "eta1", "eta2", "n_removed")]

influence_df$delta_eta1 <- influence_df$eta1 - theta_full["eta1"]
influence_df$delta_eta2 <- influence_df$eta2 - theta_full["eta2"]

influence_df$abs_delta_eta1 <- abs(influence_df$delta_eta1)
influence_df$abs_delta_eta2 <- abs(influence_df$delta_eta2)

influence_df <- influence_df[order(-influence_df$abs_delta_eta1), ]

print(influence_df)


influence_df$influence_score <- sqrt(
  (influence_df$delta_eta1 / sd(jack_valid$eta1))^2 +
  (influence_df$delta_eta2 / sd(jack_valid$eta2))^2
)

influence_df <- influence_df[order(-influence_df$influence_score), ]
print(influence_df)


ggplot(influence_df, aes(x = eta2, y = eta1, label = block)) +
  geom_point(size = 3) +
  geom_text_repel() +
  geom_point(
    aes(x = theta_full["eta2"], y = theta_full["eta1"]),
    color = "red",
    size = 4
  ) +
  labs(
    x = expression(hat(eta)[2]),
    y = expression(hat(eta)[1])
  ) +
  btf_theme
 