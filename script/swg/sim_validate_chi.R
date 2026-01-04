rm(list = ls())
cat("\014")

source("./script/load_libraries.R")

library(ggplot2)
library(dplyr)
library(data.table)
library(sf)
library(parallel)

functions_folder <- "./R"
files <- list.files(functions_folder, full.names = TRUE)
invisible(lapply(files, function(f) source(f, echo = FALSE)))

# ------------------------------------------------------------------------------
# Configuration (edit to test other parameter sets)
# ------------------------------------------------------------------------------
sim_config <- list(
  q = 0.95,
  delta = 12,
  dmin = 1200,              # meters between conditional points
  tau_vect = 0:6,           # temporal lags in 5-min units
  Nsim = 1000,
  adv_vector = c(vx = 0.2, vy = 0.1),  # advection used during simulation
  beta1 = 1.435,            # km / h^alpha1
  beta2 = 4.964,            # h units
  alpha1 = 0.183,
  alpha2 = 0.962,
  eta1 = 0.596,
  eta2 = 6.734,
  h_breaks = seq(0, 10, by = 0.5),
  threshold_latent = 1000,  # latent threshold used in sim_rpareto_coords
  seed = 42
)

convert_beta_units <- function(beta1, beta2, alpha1, alpha2,
                               c_x = 1000, c_t = 12) {
  list(
    beta1 = beta1 / (c_x^alpha1),
    beta2 = beta2 / (c_t^alpha2)
  )
}

prepare_validation_inputs <- function(cfg) {
  message("Loading rainfall and site metadata...")

  rain_file <- paste0(
    data_folder,
    "omsev/omsev_5min/rain_mtp_5min_2019_2024_cleaned.csv"
  )
  rain_df <- read.csv(rain_file)
  rownames(rain_df) <- rain_df$dates
  rain <- rain_df[, -1]

  stations_to_remove <- c("cines", "hydro", "brives")
  rain <- rain[, !(colnames(rain) %in% stations_to_remove)]

  loc_file <- paste0(data_folder, "omsev/loc_rain_gauges.csv")
  locations <- read.csv(loc_file)
  locations$Station <- c(
    "iem", "mse", "poly", "um", "cefe", "cnrs", "crbm", "archiw", "archie",
    "um35", "chu1", "chu2", "chu3", "chu4", "chu5", "chu6", "chu7",
    "cines", "brives", "hydro"
  )
  locations <- locations[!(locations$Station %in% stations_to_remove), ]

  coords_sf <- st_as_sf(
    locations[, c("Longitude", "Latitude")],
    coords = c("Longitude", "Latitude"),
    crs = 4326
  )
  coords_sf <- st_transform(coords_sf, crs = 2154)
  coords_m <- st_coordinates(coords_sf)

  grid_coords_km <- as.data.frame(coords_m / 1000)
  colnames(grid_coords_km) <- c("Longitude", "Latitude")
  rownames(grid_coords_km) <- locations$Station

  message("Estimating marginal EGPD parameters...")
  egpd_file <- paste0(
    data_folder,
    "../thesis/resources/images/EGPD/OMSEV/2019_2024/egpd_results.csv"
  )
  egpd_params <- read.csv(egpd_file)

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

  message("Building extreme episode catalogue...")
  set_st_excess <- get_spatiotemp_excess(
    rain,
    quantile = cfg$q,
    remove_zeros = TRUE
  )

  s0t0_set <- get_s0t0_pairs(
    grid_coords_km,
    rain,
    min_spatial_dist = cfg$dmin / 1000,
    episode_size = cfg$delta,
    set_st_excess = set_st_excess,
    n_max_episodes = 10000,
    latlon = FALSE
  )

  selected_points <- data.table::as.data.table(s0t0_set)
  selected_points[, t0_date := as.POSIXct(
    t0_date,
    format = "%Y-%m-%d %H:%M:%S",
    tz = "UTC"
  )]

  episodes <- get_extreme_episodes(
    selected_points,
    rain,
    episode_size = cfg$delta,
    unif = FALSE
  )$episodes

  list(
    grid_coords = grid_coords_km,
    params_margins = params_margins,
    selected_points = selected_points,
    list_episodes = episodes,
    s0_list = selected_points$s0,
    u_list = selected_points$u_s0,
    times = 0:(cfg$delta - 1)
  )
}

simulate_chi_validation <- function(inputs, cfg) {
  if (!is.null(cfg$seed)) set.seed(cfg$seed)

  params_m5 <- convert_beta_units(
    cfg$beta1,
    cfg$beta2,
    cfg$alpha1,
    cfg$alpha2
  )

  params_vario_sim <- list(
    beta1 = params_m5$beta1,
    beta2 = params_m5$beta2,
    alpha1 = cfg$alpha1,
    alpha2 = cfg$alpha2
  )

  nsim <- cfg$Nsim
  tau_vect <- cfg$tau_vect
  df_coords <- inputs$grid_coords
  list_lags <- vector("list", nsim)
  list_excesses <- vector("list", nsim)

  message("Simulating ", nsim, " conditional episodes...")
  for (i in seq_len(nsim)) {
    idx <- sample(seq_along(inputs$s0_list), 1)
    s0_name <- inputs$s0_list[idx]
    u_emp   <- inputs$u_list[idx]
    s0_coords <- df_coords[s0_name, , drop = TRUE]

    sim_episode_list <- NULL
    capture.output({
      sim_episode_list <- simulate_many_episodes(
        N = 1,
        u = cfg$threshold_latent,
        u_emp = u_emp,
        params_vario = params_vario_sim,
        params_margins = inputs$params_margins,
        coords = df_coords,
        times = inputs$times,
        adv = cfg$adv_vector,
        t0 = 0,
        s0 = s0_name
      )
    })

    sim_episode <- sim_episode_list[[1]]
    Z_ep <- sim_episode$Z
    lags <- get_conditional_lag_vectors(
      df_coords = df_coords,
      s0 = s0_coords,
      t0 = 0,
      tau_vect = tau_vect,
      latlon = FALSE
    )

    u0 <- Z_ep[s0_name, 1]
    excesses <- empirical_excesses_rpar(
      data_rain = t(Z_ep),
      threshold = u0,
      df_lags = lags,
      t0 = 0
    )

    lags$tau <- lags$tau * (5 / 60)
    list_lags[[i]] <- lags
    list_excesses[[i]] <- excesses
  }

  wind_df <- data.frame(
    vx = rep(cfg$adv_vector[["vx"]], nsim),
    vy = rep(cfg$adv_vector[["vy"]], nsim)
  )

  params_for_theory <- c(
    cfg$beta1,
    cfg$beta2,
    cfg$alpha1,
    cfg$alpha2,
    cfg$eta1,
    cfg$eta2
  )

  chi_res <- plot_th_emp_chi(
    list_lags_filtered = list_lags,
    list_excesses_filtered = list_excesses,
    wind_df_filtered = wind_df,
    params_estimates = params_for_theory,
    tau_min = 0,
    tau_fixed = cfg$tau_vect[1],
    h_breaks = cfg$h_breaks,
    latlon = FALSE,
    adv_transform = TRUE
  )

  corr_val <- with(
    chi_res$res_cmp,
    cor(chi_emp_bar, chi_theo_bar, use = "complete.obs")
  )
  attr(chi_res$res_cmp, "correlation") <- corr_val

  chi_res
}

run_sim_validation <- function(cfg = sim_config) {
  inputs <- prepare_validation_inputs(cfg)
  simulate_chi_validation(inputs, cfg)
}

if (identical(environment(), globalenv())) {
  validation <- run_sim_validation(sim_config)

  cor_val <- attr(validation$res_cmp, "correlation")
  cat(
    "\nCorrelation between empirical and theoretical chi: ",
    round(cor_val, 3),
    "\n",
    sep = ""
  )

  print(validation$plots$all)

  # Save plot
  foldername <- paste0(im_folder, "swg/sim_validate_chi/")
  if (!dir.exists(foldername)) {
    dir.create(foldername, recursive = TRUE)
  }

  filename <- paste0(
    foldername,
    "simulated_vs_theoretical_chi_adv_02_01_significant.png"
  )
  ggsave(filename,
         plot = validation$plots$all,
         width = 6,
         height = 5)


  df_plot <- validation$res_cmp

  # plot chi empirical vs theoretical with point size = n_pairs
  p <- ggplot(df_plot, aes(
    x = chi_theo_bar,
    y = chi_emp_bar,
    size = n_pairs
  )) +
    labs(
      x = "Theoretical Chi",
      y = "Empirical Chi"
    ) +
    geom_abline(linetype = "dashed", color = "red") +
    geom_point(alpha = 0.8, color = btfgreen) +
    theme_minimal() +
    btf_theme
  p

  filename <- paste0(
    foldername,
    "simulated_vs_theoretical_chi_points_adv_02_01_significant.png"
  )
  ggsave(filename,
         plot = p,
         width = 6,
         height = 5)
}
