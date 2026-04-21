################################################################################
# Simulations
################################################################################
grid_omsev <- grid_coords_km
adv_matrix <- as.matrix(adv_df_transfo[, c("vx_t", "vy_t")])
M <- 5 # repetitions per episode
Nsim <- length(list_episodes)  # number of episodes to simulate
sims_by_ep <- vector("list", Nsim)
timesteps <- 0:11  # every 5 minutes for 1 hour
for (j in seq_len(Nsim)) {
  ep_idx <- j
  s0_j  <- s0_list[ep_idx]
  u_j   <- u_list[ep_idx]
  adv_j <- adv_matrix[ep_idx, ]

  sims_by_ep[[j]] <- vector("list", M)

  for (m in seq_len(M)) {
    sims_by_ep[[j]][[m]] <- simulate_many_episodes(
      N = 1,
      u_emp = u_j,
      params_vario = params_kmh,
      params_margins = params_margins,
      coords = grid_omsev,
      times = timesteps * 5 / 60,
      adv = adv_j,
      t0 = 0,
      s0 = s0_j
    )[[1]]
  }
}

sites_names <- colnames(rain)

# Cumulatives observed and simulated
cum_obs_mat <- t(vapply(
  seq_along(list_episodes),
  function(i) sum_over_time_by_site(list_episodes[[i]], sites_names),
  FUN.VALUE = setNames(numeric(length(sites_names)), sites_names)
))

cum_obs_long <- as.data.frame(cum_obs_mat) |>
  mutate(ep_id = seq_len(nrow(cum_obs_mat))) |>
  pivot_longer(-ep_id, names_to = "site", values_to = "cum_obs")

cum_sim_long <- purrr::map_dfr(seq_along(sims_by_ep), function(j) {
  purrr::map_dfr(seq_along(sims_by_ep[[j]]), function(m) {
    simsX_jm <-sims_by_ep[[j]][[m]]$X
    simsX_jm[simsX_jm < 0.22] <- 0  # thresholding at 0.22 mm
    v <- sum_over_time_by_site(simsX_jm, sites_names)
    tibble(
      ep_id = j,
      rep   = m,
      site  = names(v),
      cum_sim = as.numeric(v)
    )
  })
})

# Get information about episodes
meta_ep <- as.data.frame(selected_points) |>
  mutate(ep_id = seq_len(nrow(selected_points))) |>
  select(ep_id, speed_class)
head(meta_ep)

# Join OBS + SIM at the (episode, site) level
df_cum <- cum_obs_long |>
  left_join(meta_ep, by = "ep_id") |>
  left_join(cum_sim_long, by = c("ep_id", "site"))

# Cumuls observed and simulated per site (median over repetitions)
df_cum_site <- df_cum |>
  group_by(ep_id, site) |>
  summarise(
    cum_obs = first(cum_obs),
    cum_sim = median(cum_sim, na.rm = TRUE),
    .groups = "drop"
  )

# Scatter plot observed vs simulated cumuls
plot(df_cum_site$cum_obs, df_cum_site$cum_sim,
     xlab = "Cumul observed (mm)",
     ylab = "Cumul simulated (mm)",
     main = "Observed vs simulated cumuls at site level")
abline(0, 1, col = "red", lty = 2)  
