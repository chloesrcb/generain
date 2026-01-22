
################################################################################
# Simulations
################################################################################
grid_omsev <- grid_coords_km
adv_matrix <- as.matrix(adv_df_transfo[, c("vx_t", "vy_t")])
M <- 5 # repetitions per episode
Nsim <- length(list_episodes)  # number of episodes to simulate
sims_by_ep <- vector("list", Nsim)

for (j in seq_len(Nsim)) {
  ep_idx <- j
  s0_j  <- s0_list[ep_idx]
  u_j   <- u_list[ep_idx]
  adv_j <- adv_matrix[ep_idx, ]

  sims_by_ep[[j]] <- vector("list", M)

  for (m in seq_len(M)) {
    sims_by_ep[[j]][[m]] <- simulate_many_episodes(
      N = 1,
      u = 1000,
      u_emp = u_j,
      params_vario = params_kmh,
      params_margins = params_margins,
      coords = grid_omsev,
      times = times * 5 / 60,
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
    v <- sum_over_time_by_site(sims_by_ep[[j]][[m]]$X, sites_names)
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

# quantiles by site and replicate
qq_by_site_rep <- df_cum |>
  filter(is.finite(cum_obs), is.finite(cum_sim)) |>
  group_by(site, rep) |>
  summarise(qq = list(make_qq_df(cum_obs, cum_sim, min_n = 100)), .groups="drop") |>
  unnest(qq)

ggplot(qq_by_site_rep, aes(q_obs, q_sim)) +
  geom_point(alpha = 0.12, size = 0.6, color=btfgreen) +
  geom_abline(slope = 1, intercept = 0, linetype = 2, color="red", alpha=0.5) +
  facet_wrap(~ site, scales = "free") +
  labs(
    x = "Quantiles observed cumuls",
    y = "Quantiles simulated cumuls",
  ) +
  theme_minimal()

# save plot
foldername_plot <- paste0(
  im_folder,
  "swg/omsev/qq_plots/"
)

filename_plot <- paste0(
  foldername_plot,
  "qq_plots_cumuls_by_site_allreps.png"
)

ggsave(
  filename = filename_plot,
  plot = last_plot(),
  width = 10,
  height = 8,
  dpi = 300
)


# Get cumuls observed and simulated per episode x site
cum_obs_mat <- t(vapply(
  seq_along(list_episodes),
  function(i) sum_over_time_by_site(list_episodes[[i]], sites_names),
  FUN.VALUE = setNames(numeric(length(sites_names)), sites_names)
))

cum_obs_long <- as.data.frame(cum_obs_mat) %>%
  mutate(ep_id = seq_len(nrow(cum_obs_mat))) %>%
  pivot_longer(-ep_id, names_to = "site", values_to = "cum_obs")

# one value per (episode, replicate, site)
cum_sim_long_rep <- purrr::map_dfr(seq_along(sims_by_ep), function(j) {
  purrr::map_dfr(seq_along(sims_by_ep[[j]]), function(m) {
    epX <- sims_by_ep[[j]][[m]]$X
    epX[epX < 0.22] <- 0  # thresholding at 0.22 mm (measurement precision)
    v <- sum_over_time_by_site(epX, sites_names)

    tibble(
      ep_id   = j,
      rep     = m,
      site    = names(v),
      cum_sim = as.numeric(v)
    )
  })
})

# Aggregate over replicates: median cum_sim per (episode, site)
cum_sim_long <- cum_sim_long_rep %>%
  group_by(ep_id, site) %>%
  summarise(
    cum_sim = median(cum_sim, na.rm = TRUE),
    .groups = "drop"
  )

meta_ep <- as.data.frame(selected_points) %>%
  mutate(ep_id = seq_len(nrow(selected_points))) %>%
  select(ep_id, s0, speed_class)

# Join OBS + SIM at the (episode, site) level
df_cum <- cum_obs_long %>%
  left_join(meta_ep, by = "ep_id") %>%
  left_join(cum_sim_long, by = c("ep_id", "site"))

cum_obs_mat <- t(vapply(
  seq_along(list_episodes),
  function(i) sum_over_time_by_site(list_episodes[[i]], sites_names),
  FUN.VALUE = setNames(numeric(length(sites_names)), sites_names)
))

cum_obs_long <- as.data.frame(cum_obs_mat) |>
  mutate(ep_id = seq_len(nrow(cum_obs_mat))) |>
  pivot_longer(-ep_id, names_to = "site", values_to = "cum_obs")

# Simulated cumuls per episode x site x replicate
cum_sim_long_rep <- purrr::map_dfr(seq_along(sims_by_ep), function(j) {
  purrr::map_dfr(seq_along(sims_by_ep[[j]]), function(m) {
    v <- sum_over_time_by_site(sims_by_ep[[j]][[m]]$X, sites_names)
    tibble(ep_id = j, rep = m, site = names(v), cum_sim = as.numeric(v))
  })
})

cum_sim_long <- cum_sim_long_rep |>
  group_by(ep_id, site) |>
  summarise(cum_sim = median(cum_sim, na.rm = TRUE), .groups = "drop")

meta_ep <- as.data.frame(selected_points) |>
  mutate(ep_id = seq_len(nrow(selected_points))) |>
  select(ep_id, speed_class)

df_cum <- cum_obs_long |>
  left_join(meta_ep, by = "ep_id") |>
  left_join(cum_sim_long, by = c("ep_id", "site"))


df_cum_rep <- cum_obs_long %>%
  left_join(meta_ep, by = "ep_id") %>%
  left_join(cum_sim_long_rep, by = c("ep_id", "site"))

qq_by_site_rep <- df_cum_rep %>%
  filter(is.finite(cum_obs), is.finite(cum_sim)) %>%
  group_by(site, rep) %>%
  summarise(
    qq = list(make_qq_df(cum_obs, cum_sim, min_n = 100)),
    .groups = "drop"
  ) %>%
  unnest(qq)

p_site_rep <- ggplot(qq_by_site_rep, aes(q_obs, q_sim)) +
  geom_point(alpha = 0.5, size = 0.6, color = btfgreen) +
  geom_abline(slope = 1, intercept = 0, linetype = 2, color = "red", alpha = 0.5) +
  facet_wrap(~ site, scales = "free") +
  labs(
    x = "Observed cumuls (quantiles)",
    y = "Simulated cumuls (quantiles)",
  ) +
  theme_minimal() +
  btf_theme

print(p_site_rep)


foldername_plot <- paste0(im_folder, "swg/omsev/qq_plots/")
if (!dir.exists(foldername_plot)) dir.create(foldername_plot, recursive = TRUE)

ggsave(
  filename = paste0(foldername_plot, "qq_cumuls_by_site_allreps.pdf"),
  plot = p_site_rep,
  width = 10, height = 8, dpi = 300
)
