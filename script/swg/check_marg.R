
################################################################################
# Plot function for split violin plots
################################################################################

# Split-violin geom
GeomSplitViolin <- ggproto(
  "GeomSplitViolin", GeomViolin,
  draw_group = function(self, data, panel_scales, coord, draw_quantiles = NULL) {
    data <- transform(data,
                      xminv = x - violinwidth * (x - xmin),
                      xmaxv = x + violinwidth * (xmax - x))
    grp <- data[1, "group"]
    newdata <- data
    newdata$x <- if (grp %% 2 == 1) data$xminv else data$xmaxv
    newdata <- newdata[order(newdata$y), ]
    newdata <- rbind(newdata[1, ], newdata, newdata[nrow(newdata), ], newdata[1, ])
    newdata[c(1, nrow(newdata)-1, nrow(newdata)), "x"] <- newdata[1, "x"]
    GeomPolygon$draw_panel(newdata, panel_scales, coord)
  }
)

geom_split_violin <- function(mapping = NULL, data = NULL, stat = "ydensity",
                              position = "identity", ..., trim = TRUE, scale = "area",
                              na.rm = FALSE, show.legend = NA, inherit.aes = TRUE) {
  layer(
    data = data, mapping = mapping, stat = stat, geom = GeomSplitViolin,
    position = position, show.legend = show.legend, inherit.aes = inherit.aes,
    params = list(trim = trim, scale = scale, na.rm = na.rm, ...)
  )
}


################################################################################
# Simulations
################################################################################
grid_omsev <- grid_coords_km
adv_matrix <- as.matrix(adv_df_transfo[, c("vx_t", "vy_t")])
Nsim <- 100
s0_sim <- integer(Nsim)
sims_all <- vector("list", Nsim)
adv_sim <- matrix(0, nrow = Nsim, ncol = 2)
for (i in seq_len(Nsim)) {
  idx <- sample(seq_along(list_episodes), 1)
  s0_i <- s0_list[idx]
  s0_sim[i] <- s0_i
  u_i <- u_list[idx]
  adv_i <- adv_matrix[idx, ]
  adv_sim[i, ] <- adv_i
  sims_all[[i]] <- simulate_many_episodes(
    N = 1,
    u = 1000,
    u_emp = u_i,
    params_vario = params_kmh,
    params_margins = params_margins,
    coords = grid_omsev, # km
    times = times * 5 / 60,  # in hours
    adv = adv_i,
    t0 = 0,
    s0 = s0_i
  )[[1]]
}

Xsim_all <- do.call(cbind, lapply(sims_all, function(sim) sim$X))
site <- rownames(sims_all[[1]]$X)[4]
lt <- ncol(sims_all[[1]]$X)

Xsim_site <- do.call(rbind, lapply(sims_all, function(sim) as.numeric(sim$X[site, ])))

Xobs_site <- do.call(rbind, lapply(list_episodes, function(ep) as.numeric(ep[, site])))

df <- data.frame(
  value = c(Xobs_site, Xsim_site),
  type = rep(c("Observed episodes", "Simulated episodes"))
) %>%
  filter(is.finite(value) & value > 0.2) %>%
  mutate(yplot = value)

ggplot(df, aes(x = type, y = yplot, fill = type)) +
  geom_violin(alpha = 0.7, trim = FALSE, color = NA) +
  labs(x = "", y = "Rain") +
  theme_minimal() + btf_theme +
  theme(legend.position = "none")

df_log <- data.frame(
  value = c(Xobs_site, Xsim_site),
  type = rep(c("Observed episodes", "Simulated episodes"))
) %>%
  filter(is.finite(value) & value > 0.1) %>%
  mutate(yplot = log1p(value))

ggplot(df_log, aes(x = type, y = yplot, fill = type)) +
  geom_violin(alpha = 0.7, trim = FALSE, color = NA) +
  labs(x = "", y = "log(1 + rain)") +
  theme_minimal() + btf_theme +
  theme(legend.position = "none")


sites <- rownames(sims_all[[1]]$X)
lt <- ncol(sims_all[[1]]$X)

Xsim_all <- lapply(sims_all, function(sim) sim$X)

sites <- rownames(sims_all[[1]]$X)
sites_plot <- sites_names[13:17]  # sites to plot

df_all <- do.call(rbind, lapply(sites_plot, function(site) {

  Xobs_site <- unlist(lapply(list_episodes, function(ep) ep[, site]))
  Xobs_site <- Xobs_site[is.finite(Xobs_site)]

  Xsim_site <- unlist(lapply(sims_all, function(sim) sim$X[site, ]))
  Xsim_site <- Xsim_site[is.finite(Xsim_site)]

  data.frame(
    site  = site,
    value = c(Xobs_site, Xsim_site),
    type  = rep(c("Observed", "Simulated"))
  )
}))


df_all <- df_all %>%
  filter(is.finite(value) & value > 0.22) %>%
  mutate(yplot = value)

# Plot
ggplot(df_all, aes(x = type, y = yplot, fill = type)) +
  geom_violin(alpha = 1, trim = FALSE, color = NA) +
  facet_wrap(~site, scales = "free_y", ncol = 3) +
  labs(
    x = "", y = "Rainfall (mm/5min)"
  ) +
  theme_minimal() + btf_theme +
  theme(legend.position = "none")

ggsave(
  filename = paste0(
    foldername_plot,
    "marginal_distributions_observed_vs_simulated_episodes_",
    "_n", Nsim, "_3",
    ".png"
  ),
  plot = last_plot(),
  width = 10,
  height = 8,
  dpi = 300
)