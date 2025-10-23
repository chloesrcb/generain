################################################################################
# FULL ETA OPTIMIZATION PIPELINE — CLEAN VERSION
# Tests eta1 = 0 vs. eta1 = eta1_com for various (delta, dmin, q)
# Author: you, refactored by ChatGPT
################################################################################

library(numDeriv)
library(parallel)
library(data.table)
library(dplyr)
library(lubridate)
library(sf)

################################################################################
# 1. LOAD FUNCTIONS AND COMMON DATA
################################################################################

functions_folder <- "./R"
files <- list.files(functions_folder, full.names = TRUE)
# # load all functions in files
invisible(lapply(files, function(f) source(f, echo = FALSE)))

filename_loc <- paste0(data_folder,
                           "omsev/loc_rain_gauges.csv")
# get location of each rain gauge
location_gauges <- read.csv(filename_loc)
location_gauges$Station <- c("iem", "mse", "poly", "um", "cefe", "cnrs",
                             "crbm", "archiw", "archie", "um35", "chu1",
                             "chu2", "chu3", "chu4", "chu5", "chu6", "chu7",
                             "cines", "brives", "hydro")

# get rain data
filename_rain <- paste0(data_folder, "omsev/omsev_5min/rain_mtp_5min_2019_2024_cleaned.csv")
rain_omsev <- read.csv(filename_rain)


rain_test <- rain_omsev
rain_test$nb_sites_non_NA <- apply(rain_test[ , -1], 1, function(x) sum(!is.na(x)))
date_2_sites <- rain_test$dates[which(rain_test$nb_sites_non_NA >= 2)[1]]
date_3_sites <- rain_test$dates[which(rain_test$nb_sites_non_NA >= 3)[1]]
date_4_sites <- rain_test$dates[which(rain_test$nb_sites_non_NA >= 4)[1]]
date_5_sites <- rain_test$dates[which(rain_test$nb_sites_non_NA >= 5)[1]]

# begin in 2020-01-01
rain_omsev <- rain_omsev[rain_omsev$dates >= date_5_sites, ]
head(rain_omsev)

# put dates as rownames
rownames(rain_omsev) <- rain_omsev$dates
rain <- rain_omsev[-1] # remove dates column
# rain$mse
# Get distances matrix
dist_mat <- get_dist_mat(location_gauges)
df_dist <- reshape_distances(dist_mat)

# remove cines, hydro, brives
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
dist_mat <- get_dist_mat(grid_coords_km, latlon = FALSE)

################################################################################
# 2. BUILD EPISODE DATA FUNCTION
################################################################################

get_episode_data <- function(q, delta, dmin,
                             data,
                             grid_coords_km,
                             dist_mat = NULL,
                             data_folder,
                             tau_vect = 0:10,,
                             adv_name = c("combined_comephore_omsev"),
                             max_episodes = 10000,
                             latlon = FALSE,
                             match_tol_sec = 5*60) {

  # --- Distance matrix and hmax (km)
  df_dist <- reshape_distances(dist_mat)
  df_dist$value <- df_dist$value
  hmax <- max(df_dist$value, na.rm = TRUE)

  # --- Compute spatio-temporal exceedances
  st_excess <- get_spatiotemp_excess(data, quantile = q, remove_zeros = TRUE)

  # --- Choose conditional (s0, t0, u_s0)
  s0t0_set <- get_s0t0_pairs(
    grid_coords_km,
    data = data,
    min_spatial_dist = dmin,
    episode_size = delta,
    set_st_excess = st_excess,
    n_max_episodes = max_episodes,
    latlon = FALSE
  )

  s0t0_set <- s0t0_set %>%
    mutate(t0_date = as.POSIXct(t0_date, format = "%Y-%m-%d %H:%M:%S", tz = "UTC"))

  # --- Build raw extreme episodes
  episodes_obj <- get_extreme_episodes(
    selected_points = s0t0_set,
    data = data,
    episode_size = delta,
    unif = FALSE
  )
  list_episodes <- episodes_obj$episodes

  # --- Load advection and match to selected episodes
  adv_filename <- paste0(
    data_folder, "/omsev/adv_estim/", adv_name[1],
    "/episode_advection_q", q * 100,
    "_delta", delta,
    "_dmin", dmin, ".csv"
  )
  adv_df <- read.csv(adv_filename, sep = ",")
  adv_df$t0_omsev <- as.POSIXct(adv_df$t0_omsev, format = "%Y-%m-%d %H:%M:%S", tz = "UTC")

  setDT(s0t0_set)
  setDT(adv_df)
  setkey(adv_df, t0_omsev)
  selected_episodes <- adv_df[s0t0_set, roll = match_tol_sec, on = .(t0_omsev = t0_date)]
  colnames(selected_episodes)[colnames(selected_episodes) == "dx_comb_kmh"] <- "adv_x"
  colnames(selected_episodes)[colnames(selected_episodes) == "dy_comb_kmh"] <- "adv_y"

  selected_episodes_nona <- selected_episodes[!is.na(adv_x) & !is.na(adv_y)]
  if (nrow(selected_episodes_nona) == 0)
    stop("No episodes with matched wind — check advection file or time tolerance.")

  V_episodes <- data.frame(vx = selected_episodes_nona$adv_x,
                           vy = selected_episodes_nona$adv_y)

  # --- Build lag vectors and empirical excesses
  df_coords <- as.data.frame(grid_coords_km)
  s0_list <- selected_episodes_nona$s0
  u0_list <- selected_episodes_nona$u_s0

  episodes_obj <- get_extreme_episodes(
    selected_points = selected_episodes_nona,
    data = data,
    episode_size = delta,
    unif = FALSE
  )
  list_episodes <- episodes_obj$episodes

  list_results <- mclapply(seq_along(s0_list), function(i) {
    s0 <- s0_list[i]
    col_s0 <- which(rownames(df_coords) == s0)
    s0_coords <- df_coords[col_s0, , drop = FALSE]
    episode <- list_episodes[[i]]
    ind_t0_ep <- 0

    lags <- get_conditional_lag_vectors(df_coords, s0_coords, ind_t0_ep,
                                        tau_vect, latlon = FALSE)
    u <- u0_list[i]

    excesses <- empirical_excesses_rpar(
      episode, quantile = u, threshold = TRUE, df_lags = lags, t0 = ind_t0_ep
    )
    lags$tau <- lags$tau * 5 / 60
    list(lags = lags, excesses = excesses)
  }, mc.cores = max(1, detectCores() - 1))

  list_lags <- lapply(list_results, `[[`, "lags")
  list_excesses <- lapply(list_results, `[[`, "excesses")

  return(list(
    list_episodes = list_episodes,
    list_excesses = list_excesses,
    list_lags = list_lags,
    wind_df = V_episodes,
    hmax = hmax,
    meta = list(n_total = nrow(s0t0_set),
                n_with_wind = nrow(selected_episodes_nona))
  ))
}

################################################################################
# 3. ETA OPTIMIZATION FUNCTION (FOR ONE CASE)
################################################################################

run_eta_optim <- function(delta, dmin, q, params_com,
                          wind_df, list_lags, list_episodes, list_excesses, hmax) {

  eta1_com <- params_com[5]
  eta2_com <- params_com[6]
  init_params <- params_com[1:4]

  cases <- list(
    list(name = "eta1=0", fixed_eta1 = 0, fixed_eta2 = eta2_com),
    list(name = "eta1=com", fixed_eta1 = eta1_com, fixed_eta2 = eta2_com)
  )

  results <- list()

  for (case in cases) {
    cat("\n--------------------------------------------------\n")
    cat("Δ =", delta, "| dmin =", dmin, "| q =", q, "|", case$name, "\n")
    cat("--------------------------------------------------\n")

    objfun <- function(p) {
      neg_ll_composite_fixed_eta(
        params = p,
        list_episodes = list_episodes,
        list_excesses = list_excesses,
        list_lags = list_lags,
        wind_df = wind_df,
        hmax = hmax,
        latlon = FALSE,
        distance = "lalpha",
        fixed_eta1 = case$fixed_eta1,
        fixed_eta2 = case$fixed_eta2
      )
    }

    opt <- optim(
      par = init_params,
      fn = objfun,
      method = "L-BFGS-B",
      lower = c(1e-8, 1e-8, 1e-8, 1e-8),
      upper = c(10, 10, 1.999, 1.999),
      control = list(maxit = 20000, trace = 0)
    )

    grad_val <- grad(objfun, opt$par)
    val_check <- objfun(opt$par)
    stopifnot(isTRUE(all.equal(val_check, opt$value, tolerance = 1e-6)))

    results[[case$name]] <- list(
      delta = delta,
      dmin = dmin,
      q = q,
      case = case$name,
      params = opt$par,
      negll = opt$value,
      grad = grad_val,
      grad_norm = sqrt(sum(grad_val^2)),
      conv = opt$convergence
    )

    cat("→ NegLL:", round(opt$value, 2), "\n")
    cat("→ Gradient norm:", round(sqrt(sum(grad_val^2)), 5), "\n")
    cat("→ Convergence code:", opt$convergence, "\n")
  }

  return(results)
}

################################################################################
# 4. MAIN LOOP OVER (delta, dmin, q)
################################################################################

delta_values <- c(12, 15)
dmin_values  <- c(5, 10)
q_values     <- c(0.95, 0.99)

com_results <- read.csv(paste0(
  data_folder, "/comephore/optim_results/lalpha/free_eta/combined_optim_results.csv"
))

all_results <- list()
iter <- 1

for (delta in delta_values) {
  for (dmin in dmin_values) {
    for (q in q_values) {
      cat("\n### Running combination (Δ =", delta, ", dmin =", dmin, ", q =", q, ") ###\n")

      params_com <- com_results[com_results$q == q * 100 &
                                  com_results$delta == 30 &
                                  com_results$dmin == 5, ]

      if (nrow(params_com) == 0) {
        warning("No comephore params for this triplet, skipping.")
        next
      }

      params_com <- c(params_com$beta1, params_com$beta2,
                      params_com$alpha1, params_com$alpha2,
                      params_com$eta1, params_com$eta2)

      edata <- get_episode_data(
        q = q,
        delta = delta,
        dmin = dmin,
        rain = rain,
        grid_coords_km = grid_coords_km,
        dist_mat_m = dist_mat,
        data_folder = data_folder,
        tau_vect = 0:10,
        adv_name = "combined_comephore_omsev",
        max_episodes = 10000,
        latlon = FALSE
      )

      res <- run_eta_optim(
        delta = delta,
        dmin = dmin,
        q = q,
        params_com = params_com,
        wind_df = edata$wind_df,
        list_lags = edata$list_lags,
        list_episodes = edata$list_episodes,
        list_excesses = edata$list_excesses,
        hmax = edata$hmax
      )

      all_results[[iter]] <- res
      iter <- iter + 1
    }
  }
}

################################################################################
# 5. ORGANIZE AND EXPORT RESULTS
################################################################################

df_res <- do.call(rbind, lapply(all_results, function(x) {
  do.call(rbind, lapply(x, function(r) {
    data.frame(
      delta = r$delta,
      dmin = r$dmin,
      q = r$q,
      case = r$case,
      negll = r$negll,
      grad_norm = r$grad_norm,
      beta1 = r$params[1],
      beta2 = r$params[2],
      alpha1 = r$params[3],
      alpha2 = r$params[4],
      conv = r$conv
    )
  }))
}))

df_summary <- df_res %>%
  group_by(delta, dmin, q) %>%
  summarize(
    NLL_eta0 = negll[case == "eta1=0"],
    NLL_etaCom = negll[case == "eta1=com"],
    Diff_NLL = NLL_etaCom - NLL_eta0,
    .groups = "drop"
  )

print(df_summary)
write.csv(df_res, file = "eta_comparison_results.csv", row.names = FALSE)
