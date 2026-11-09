rm(list = ls())

# clear console
cat("\014")


# -----------------------------------------------------------------------------
# CONFIGURATIONS
# -----------------------------------------------------------------------------
qs <- c(0.95, 0.96, 0.97, 0.98, 0.99)
min_spatial_dists <- c(800, 1000, 1200, 1300)    # mètres
deltas <- c(6, 8, 12, 15, 24)       # 5 minutes


muse <- FALSE
seed <- 1000L
n_max_episodes <- 10000L
nworkers <- min(8L, max(1L, parallel::detectCores() - 1L))

# Paramètres fixés et initialisation identiques au script d'origine
com_results <- c(0.1579947, 0.9923270, 0.5800128, 0.6627969,
                 3.8962660, 2.2208320)
eta1 <- com_results[5]
eta2 <- com_results[6]
init_par <- c(1.08, 3.9, 0.12, 0.65)
par_names <- c("beta1", "beta2", "alpha1", "alpha2")

if (muse) {
  setwd("/home/serrec/work_rainstsimu/rpareto/")
  source("load_libraries.R")
  source("config_omsev.R")
} else {
  source("./script/load_libraries.R")
  source("./script/optimisation/config_omsev.R")
}

invisible(lapply(list.files("./R", full.names = TRUE, pattern = "\\.[Rr]$"),
                 function(f) source(f, echo = FALSE)))
source("./R/advection.R")

location_gauges <- read.csv(file.path(data_folder, "omsev/loc_rain_gauges.csv"))
rain_omsev <- read.csv(file.path(
  data_folder, "omsev/omsev_5min/rain_mtp_5min_2019_2025.csv"
))
rownames(rain_omsev) <- rain_omsev$dates
rain <- rain_omsev[, -1, drop = FALSE]

sites_coords <- location_gauges[, c("Longitude", "Latitude")]
rownames(sites_coords) <- location_gauges$Station
sites_coords_sf <- sf::st_as_sf(sites_coords, coords = c("Longitude", "Latitude"), crs = 4326)
sites_coords_sf <- sf::st_transform(sites_coords_sf, crs = 2154)
coords_m <- sf::st_coordinates(sites_coords_sf)
grid_coords_m <- as.data.frame(coords_m)
colnames(grid_coords_m) <- c("Longitude", "Latitude")
rownames(grid_coords_m) <- rownames(sites_coords)
grid_coords_km <- grid_coords_m / 1000
df_coords <- as.data.frame(grid_coords_km)

dist_mat <- get_dist_mat(location_gauges)
hmax <- max(dist_mat, na.rm = TRUE)

out_dir <- file.path(data_folder, "omsev/optim_results/configurations")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
random_dir <- file.path(data_folder, "omsev/random")
dir.create(random_dir, recursive = TRUE, showWarnings = FALSE)
summary_file <- file.path(out_dir, "variogram_all_configurations.csv")

excess_cache <- new.env(parent = emptyenv())

empty_result <- function(q, dmin, delta) {
  data.frame(q = q, min_spatial_dist = dmin, delta = delta,
             n_selected = NA_integer_, n_episodes = NA_integer_,
             beta1 = NA_real_, beta2 = NA_real_,
             alpha1 = NA_real_, alpha2 = NA_real_,
             eta1 = eta1, eta2 = eta2,
             neg_loglik = NA_real_, convergence = NA_integer_,
             message = NA_character_, elapsed_seconds = NA_real_)
}

estimate_one <- function(q, dmin, delta) {
  start_time <- Sys.time()
  ans <- empty_result(q, dmin, delta)
  cat(sprintf("\n--- q=%s | dmin=%sm | delta=%s ---\n", q, dmin, delta))

  tryCatch({
    q_key <- as.character(q)
    if (!exists(q_key, envir = excess_cache, inherits = FALSE)) {
      assign(q_key, get_spatiotemp_excess(
        rain, quantile = q, remove_zeros = TRUE
      ), envir = excess_cache)
    }
    set_st_excess <- get(q_key, envir = excess_cache)

    selected_points <- get_s0t0_pairs(
      grid_coords_m, rain,
      min_spatial_dist = dmin,
      episode_size = delta,
      set_st_excess = set_st_excess,
      n_max_episodes = n_max_episodes,
      latlon = FALSE
    )
    selected_points <- data.table::as.data.table(selected_points)
    ans$n_selected <- nrow(selected_points)
    if (!nrow(selected_points)) stop("Aucun épisode sélectionné")
    selected_points[, t0_date := as.POSIXct(t0_date,
       format = "%Y-%m-%d %H:%M:%S", tz = "UTC")]

    # get_advection lit les listes des t0 enregistrées sur disque.
    t0_unique <- unique(selected_points$t0_date)
    prefix <- sprintf("_q%s_delta%s_dmin%s.csv", q * 100, delta, dmin)
    write.csv(data.frame(t0_date = unique(lubridate::ceiling_date(t0_unique, "hour"))),
              file.path(random_dir, paste0("t0_episodes", prefix)), row.names = FALSE)
    write.csv(data.frame(t0_date = t0_unique),
              file.path(random_dir, paste0("t0_5min_episodes", prefix)), row.names = FALSE)

    future::plan(future::multisession, workers = nworkers)
    adv_df <- tryCatch(get_advection(q = q, delta = delta,
                         dmin = dmin, data_folder = data_folder, force = TRUE),
                       finally = future::plan(future::sequential))
    adv_df <- data.table::as.data.table(adv_df)
    adv_df[, t0_omsev := as.POSIXct(t0_omsev,
       format = "%Y-%m-%d %H:%M:%S", tz = "UTC")]

    # Association exacte des dates (comme la jointure finale du script source).
    adv_df_u <- adv_df[, .(vx_final = vx_final[1], vy_final = vy_final[1]),
                       by = t0_omsev]
    selected_episodes <- merge(
      selected_points, adv_df_u,
      by.x = "t0_date", by.y = "t0_omsev", all.x = TRUE, sort = FALSE
    )
    data.table::setnames(selected_episodes,
                         c("vx_final", "vy_final"), c("adv_x", "adv_y"))
    selected_episodes <- selected_episodes[
      is.finite(adv_x) & is.finite(adv_y)]
    ans$n_episodes <- nrow(selected_episodes)
    if (!nrow(selected_episodes)) stop("Aucun épisode avec advection valide")

    # Épisodes, lags et excès doivent rester alignés dans le même ordre.
    list_episodes <- get_extreme_episodes(
      selected_episodes, rain, episode_size = delta, unif = FALSE
    )$episodes
    if (length(list_episodes) != nrow(selected_episodes))
      stop("Nombre d'épisodes incohérent")

    u_by_site <- unique(data.frame(
      site = unlist(set_st_excess$list_s),
      threshold = unlist(set_st_excess$list_u)
    ))
    u_vec <- stats::setNames(u_by_site$threshold, u_by_site$site)

    # Le lag temporel ne doit pas dépasser la taille d'un épisode.
    tau_vect <- 0:min(10L, delta - 1L)
    n_cores_lags <- max(1L, parallel::detectCores() - 1L)
    process_episode <- function(i) {
      s0 <- selected_episodes$s0[i]
      s0_coords <- df_coords[s0, , drop = FALSE]
      if (!nrow(s0_coords)) stop("Station absente des coordonnées : ", s0)
      lags <- get_conditional_lag_vectors(
        df_coords, s0_coords, 0, tau_vect, latlon = FALSE
      )
      excesses <- empirical_excesses_rpar(
        data_rain = list_episodes[[i]], df_lags = lags,
        thresholds_by_site = u_vec, t0 = 0
      )
      lags$tau <- lags$tau * 5 / 60  # heures
      list(lags = lags, excesses = excesses)
    }
    if (.Platform$OS.type == "windows") {
      precomputed <- lapply(seq_along(list_episodes), process_episode)
    } else {
      precomputed <- parallel::mclapply(seq_along(list_episodes), process_episode,
                                       mc.cores = n_cores_lags)
    }
    if (any(vapply(precomputed, inherits, logical(1), "try-error")))
      stop("Erreur dans le calcul d'un épisode")
    list_lags <- lapply(precomputed, `[[`, "lags")
    list_excesses <- lapply(precomputed, `[[`, "excesses")
    V_episodes <- data.frame(vx = selected_episodes$adv_x,
                             vy = selected_episodes$adv_y)

    result <- optim(
      par = init_par, fn = neg_ll_composite_fixed_eta,
      list_lags = list_lags, list_episodes = list_episodes,
      list_excesses = list_excesses, hmax = hmax,
      wind_df = V_episodes, latlon = FALSE, distance = "euclidean",
      fixed_eta1 = eta1, fixed_eta2 = eta2,
      method = "L-BFGS-B",
      lower = c(1e-8, 1e-8, 1e-8, 1e-8),
      upper = c(10, 10, 1.999, 1.999),
      control = list(maxit = 20000, trace = 0)
    )
    for (j in seq_along(par_names)) ans[[par_names[j]]] <- result$par[j]
    ans$neg_loglik <- result$value
    ans$convergence <- result$convergence
    ans$message <- if (result$convergence == 0) "OK" else "Non-convergence"
  }, error = function(e) {
    ans$message <<- conditionMessage(e)
    message("Échec : ", conditionMessage(e))
  })

  ans$elapsed_seconds <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))
  ans
}


configurations <- expand.grid(
  q = qs, min_spatial_dist = min_spatial_dists, delta = deltas,
  KEEP.OUT.ATTRS = FALSE
)

# if (file.exists(summary_file)) {
#   all_results <- read.csv(summary_file, stringsAsFactors = FALSE)
# } else {
all_results <- empty_result(qs[1], min_spatial_dists[1], deltas[1])[0, ]
# }
# all_results <- empty_result(qs[1], min_spatial_dists[1], deltas[1])[0, ]

for (i in seq_len(nrow(configurations))) {
  cfg <- configurations[i, ]
  # Ignore les configurations déjà réussies ; les échecs sont retentés.
  already_done <- with(all_results,
    q == cfg$q & min_spatial_dist == cfg$min_spatial_dist &
      delta == cfg$delta & !is.na(convergence) & convergence == 0)
  if (any(already_done, na.rm = TRUE)) {
    message("Already estimated : q=", cfg$q, " dmin=", cfg$min_spatial_dist,
            " delta=", cfg$delta)
    next
  }

  result_row <- estimate_one(cfg$q, cfg$min_spatial_dist, cfg$delta)
  existing <- with(all_results,
    q == cfg$q & min_spatial_dist == cfg$min_spatial_dist & delta == cfg$delta)
  all_results <- rbind(all_results[!existing, , drop = FALSE], result_row)
  all_results <- all_results[order(all_results$q,
                      all_results$min_spatial_dist, all_results$delta), ]
  write.csv(all_results, summary_file, row.names = FALSE)
  print(result_row)
}



all_results <- read.csv(summary_file, stringsAsFactors = FALSE)
library(knitr)

all_results_table <- all_results %>%
  dplyr::select(q, min_spatial_dist, delta, n_episodes,
                beta1, beta2, alpha1, alpha2,
                neg_loglik) 
latex_table <- kable(
  all_results_table,
  format = "latex",
  booktabs = TRUE,
  digits = 3,
  row.names = FALSE,
  caption = "Variogram parameter estimates for different episode configurations",
  label = "tab:variogram_results"
)

cat(latex_table)




library(dplyr)
library(knitr)
library(kableExtra)

table_results <- all_results %>%
  filter(q == 0.95) %>%
  arrange(delta, min_spatial_dist) %>%
  select(
    delta, min_spatial_dist, n_episodes,
    beta1, beta2, alpha1, alpha2
  )

latex_table <- table_results %>%
  kbl(
    format = "latex",
    booktabs = TRUE,
    digits = 3,
    row.names = FALSE,
    escape = FALSE,
    col.names = c(
      "$\\delta$",
      "$d_{\\min}$ (m)",
      "$N$",
      "$\\beta_1$",
      "$\\beta_2$",
      "$\\alpha_1$",
      "$\\alpha_2$"
    ),
    caption = "Variogram parameter estimates for $q=0.95$.",
    label = "variogram_sensitivity_delta"
  ) %>%
  kable_styling(
    latex_options = "hold_position",
    font_size = 9
  ) %>%
  collapse_rows(
    columns = 1,
    valign = "middle",
    latex_hline = "major",
    custom_latex_hline = 1
  )

cat(latex_table)


# same by fixing delta to 12 and varying q and dmin

table_results_q <- all_results %>%
  filter(delta == 12) %>%
  arrange(q, min_spatial_dist) %>%
  select(
    q, min_spatial_dist, n_episodes,
    beta1, beta2, alpha1, alpha2
  )

latex_table_q <- table_results_q %>%
  kbl(
    format = "latex",
    booktabs = TRUE,
    digits = 3,
    row.names = FALSE,
    escape = FALSE,
    col.names = c(
      "$q$",
      "$d_{\\min}$ (m)",
      "$N$",
      "$\\beta_1$",
      "$\\beta_2$",
      "$\\alpha_1$",
      "$\\alpha_2$"
    ),
    caption = "Variogram parameter estimates for $\\delta=12$.",
    label = "variogram_sensitivity_q"
  ) %>%
  kable_styling(
    latex_options = "hold_position",
    font_size = 9
  ) %>%
  collapse_rows(
    columns = 1,
    valign = "middle",
    latex_hline = "major",
    custom_latex_hline = 1
  ) 

cat(latex_table_q)




library(dplyr)
library(tidyr)
library(ggplot2)

params_long <- all_results %>%
  filter(delta == 12, convergence == 0) %>%
  mutate(dmin = factor(min_spatial_dist)) %>%
  pivot_longer(
    cols = c(beta1, beta2, alpha1, alpha2),
    names_to = "parameter",
    values_to = "estimate"
  ) %>%
  mutate(
    parameter = factor(
      parameter,
      levels = c("beta1", "beta2", "alpha1", "alpha2"),
      labels = c("beta[1]", "beta[2]",
                 "alpha[1]", "alpha[2]")
    )
  )


params_long$parameter <- factor(
  params_long$parameter,
  levels = c("beta[1]", "beta[2]", "alpha[1]", "alpha[2]")
)

math_labels <- as_labeller(
  c(
    "beta[1]" = "hat(beta)[1]",
    "beta[2]" = "hat(beta)[2]",
    "alpha[1]" = "hat(alpha)[1]",
    "alpha[2]" = "hat(alpha)[2]"
  ),
  default = label_parsed
)

p_q <- ggplot(
  params_long,
  aes(x = q, y = estimate,
      color = dmin, group = dmin)
) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 2) +
  facet_wrap(
    ~ parameter,
    scales = "free_y",
    labeller = math_labels,
    ncol = 2
  ) +
  scale_x_continuous(breaks = seq(0.95, 0.99, 0.01)) +
  labs(
    x = "Threshold quantile q",
    y = "Parameter estimate",
    color = expression(d[min] ~ "(m)")
  ) +
  btf_theme +
  theme(
    legend.position = "bottom",
    strip.background = element_blank(),
    strip.text = element_text(
      size = 18,
      face = "plain",
      margin = margin(t = 12, b = 12)
    ),
    strip.placement = "outside",
    panel.spacing = unit(1.5, "lines"),
    panel.grid.minor = element_blank()
  )

print(p_q)

ggsave(
  "sensitivity_q.pdf",
  p_q,
  width = 10,
  height = 8
)


# same plot but fixing q to 0.95 and varying delta and dmin
params_long <- all_results %>%
  filter(q == 0.95, convergence == 0) %>%
  mutate(dmin = factor(min_spatial_dist)) %>%
  pivot_longer(
    cols = c(beta1, beta2, alpha1, alpha2),
    names_to = "parameter",
    values_to = "estimate"
  ) %>%
  mutate(
    parameter = factor(
      parameter,
      levels = c("beta1", "beta2", "alpha1", "alpha2"),
      labels = c("beta[1]", "beta[2]",
                 "alpha[1]", "alpha[2]")
    )
  )

math_labels <- as_labeller(
  c(
    "beta[1]" = "hat(beta)[1]",
    "beta[2]" = "hat(beta)[2]",
    "alpha[1]" = "hat(alpha)[1]",
    "alpha[2]" = "hat(alpha)[2]"
  ),
  default = label_parsed
)

ylims <- params_long %>%
  group_by(parameter) %>%
  summarise(ymin = min(estimate, na.rm = TRUE),
            ymax = max(estimate, na.rm = TRUE)) %>%
  ungroup()

library(ggh4x)
library(dplyr)

# Limites communes calculées sur toutes les configurations
ylims <- params_long %>%
  group_by(parameter) %>%
  summarise(
    ymin = min(estimate, na.rm = TRUE),
    ymax = max(estimate, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    padding = 0.2 * (ymax - ymin),
    ymin = ymin - padding,
    ymax = ymax + padding
  )

# Une échelle Y par paramètre
y_scales <- lapply(seq_len(nrow(ylims)), function(i) {
  scale_y_continuous(
    limits = c(ylims$ymin[i], ylims$ymax[i])
  )
})

p_delta <- ggplot(
  params_long %>% filter(q == 0.95),
  aes(x = delta, y = estimate,
      color = dmin, group = dmin)
) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 2) +
  facet_wrap(
    ~ parameter,
    scales = "free_y",
    labeller = math_labels,
    ncol = 2
  ) +
  ggh4x::facetted_pos_scales(
    y = y_scales
  ) +
  scale_x_continuous(breaks = deltas) +
  labs(
    x = "Episode size (5-min intervals)",
    y = "Parameter estimate",
    color = expression(d[min] ~ "(m)")
  ) +
  btf_theme +
  theme(
    legend.position = "bottom",
    strip.background = element_blank(),
    strip.text = element_text(
      size = 18,
      face = "plain",
      margin = margin(t = 12, b = 12)
    ),
    strip.placement = "outside",
    panel.spacing = unit(1.5, "lines"),
    panel.grid.minor = element_blank()
  )

print(p_delta)

print(p_delta)

foldername <- paste0(im_folder, "/optim/episode_config_omsev")
if (!dir.exists(foldername)) {
  dir.create(foldername, recursive = TRUE)
}
ggsave(
    file.path(foldername, "sensitivity_delta.pdf"),
    p_delta,
    width = 10,
    height = 8
    )
