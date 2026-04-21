
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
grid_omsev <- grid_coords_m
adv_matrix <- as.matrix(adv_df_transfo[, c("vx_t", "vy_t")])
Nsim <- length(list_episodes)
s0_sim <- integer(Nsim)
sims_all <- vector("list", Nsim)
adv_sim <- matrix(0, nrow = Nsim, ncol = 2)
for (i in seq_len(Nsim)) {
  # idx <- sample(seq_along(list_episodes), 1)
  idx <- i
  s0_i <- s0_list[idx]
  s0_sim[i] <- s0_i
  u_i <- u_list[idx]
  adv_i <- adv_matrix[idx, ]
  adv_sim[i, ] <- adv_i
  sims_all[[i]] <- simulate_many_episodes(
    N = 1,
    u_emp = u_i,
    params_vario = params_m5min,
    params_margins = params_margins,
    coords = grid_omsev, # km
    times = times,  # in hours
    adv = adv_i,
    t0 = 0,
    s0 = s0_i
  )[[1]]
}




Xsim_all <- do.call(cbind, lapply(sims_all, function(sim) sim$X))
# q <- c(0.1,0.25,0.5,0.75,0.9,0.95,0.99)
# apply_measurement <- function(x, p) {
#   x[x < 1e-4] <- 0
#   x[x > 0 & x < p] <- p
#   x
# }

# Qobs <- quantile(Xobs_site[Xobs_site>0], q, na.rm=TRUE)
# Qsim <- quantile(Xsim_site[Xsim_site>0.06], q, na.rm=TRUE)
# Xsim_meas <- apply_measurement(Xsim_site, p = 0.2152)
# Xsim_pos <-  Xsim_meas[Xsim_meas > 0.001]
# Qsim_m <- quantile(Xsim_pos, q, na.rm=TRUE)

# cbind(q=q, obs=Qobs, sim=Qsim, sim_measured=Qsim_m,
#       diff_before=Qsim-Qobs, diff_after=Qsim_m-Qobs)
apply_measurement <- function(x, p) {
  # x[x < 1e-4] <- 0
  x[x > 0 & x < p] <- p
  x
}

apply_disc <- function(x, p) {
  x[x < 1e-4] <- 0
  x[x > 0 & x < p] <- p
  x[x > p & x < 2*p] <- 2*p
  x
}

site_names <- colnames(sims_all[[1]]$X)

for (site in site_names) {
  # sim <- sims_all[[1]]
  Xsim_site <- unlist(lapply(sims_all, function(sim) as.numeric(sim$X[, site])))
  Xobs_site <- unlist(lapply(list_episodes, function(ep) as.numeric(ep[, site])))
  Xobs_site <- Xobs_site[!is.na(Xobs_site)]
  xr <- range(c(Xobs_site, Xsim_site), finite = TRUE)
  dens     <- density(Xobs_site, from = xr[1], to = xr[2], na.rm = TRUE)
  dens_sim <- density(Xsim_site, from = xr[1], to = xr[2], na.rm = TRUE)

  Xobs_site <- Xobs_site[Xobs_site > 0]
  p <- min(Xobs_site[Xobs_site > 0])
  # # put all values between 0 and 0.001 to 0 
  # Xsim_site[Xsim_site < 0.1] <- 0
  # and all values between 0.001 and 0.2 to 0.2
  # Xsim_site[Xsim_site >= 0.001 & Xsim_site < p] <- p
  Xsim_site <- Xsim_site[Xsim_site >  1e-3]

  Xsim_meas <- apply_disc(Xsim_site, p = 0.2152)
  Xsim_meas <- Xsim_site
  Xobs_meas <- Xobs_site

  # éventuellement retirer les 0 pour comparer les positives
  Xsim_site <- Xsim_meas[Xsim_meas > 0]
  Xobs_site <- Xobs_meas[Xobs_meas > 0]


  config <- "above0_simdisc_p02152"
  df <- rbind(
    data.frame(type = "Observed episodes", value = Xobs_site),
    data.frame(type = "Simulated episodes", value = Xsim_site)
  )
  bw_common <- bw.nrd0(df$value)

  # psite <- ggplot(df, aes(x = type, y = value, fill = type)) +
  #   geom_violin(alpha = 0.7, trim = F, scale = "width",
  #               bw = bw_common) +
  #   labs(x = "", y = "Rainfall (mm/5min)") +
  #   theme_minimal() + btf_theme +
  #   theme(legend.position = "none") 

  # # save plot
  # ggsave(paste0(im_folder, "swg/omsev/margins/violin_all", site, ".png"), psite, width = 6, height = 4)
  df_log <- df %>% mutate(yplot = log1p(value))
  # df_log <- df_log[df_log$yplot > 0,]
  psite_log <- ggplot(df_log, aes(x = type, y = yplot, fill = type)) +
    geom_violin(alpha = 0.7, trim = F, color = NA, bw = bw_common, scale = "width", bounds= c(0, Inf)) +
    labs(x = "", y = "log(1 + X)") +
    # scale_y_log10() +
    theme_minimal() + btf_theme +
    geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA) +
    theme(legend.position = "none")
  
  # save plot
  folder_site <- paste0(im_folder, "swg/omsev/margins/", config, "/")
  if (!dir.exists(folder_site)) {
    dir.create(folder_site, recursive = TRUE)
  }
  ggsave(paste0(folder_site, "violin_log_", config, "_", site, ".png"), psite_log, width = 6, height = 4)
}



# remove episode with more than 80% of active sites ie less than 20% of NA
list_episodes_filtered <- list_episodes[sapply(list_episodes, function(ep) mean(is.na(ep)) < 0.8)]
length(list_episodes_filtered) # 17 episodes left


cum_obs <- sapply(list_episodes_filtered, function(ep) sum(ep, na.rm = TRUE))
cum_sim <- sapply(sims_all, function(sim) sum(sim$X, na.rm = TRUE))

dens_cum_obs <- density(cum_obs, from = 0, to = max(cum_obs), na.rm = TRUE)
dens_cum_sim <- density(cum_sim, from = 0, to = max(cum_sim), na.rm = TRUE)
plot(dens_cum_obs, main = "Cumulative rainfall distribution", xlab = "Cumulative rainfall (mm)", ylab = "Density", col = "blue")
lines(dens_cum_sim, add = TRUE, col = "red")


dens_cum_obs$y <- dens_cum_obs$y / max(dens_cum_obs$y)
dens_cum_sim$y <- dens_cum_sim$y / max(dens_cum_sim$y)

plot(dens_cum_obs, col = "blue",
     main = "Normalized cumulative rainfall distribution",
     xlab = "Cumulative rainfall (mm)", ylab = "Relative density")
lines(dens_cum_sim, col = "red")
legend("topright", legend = c("Observed", "Simulated"),
       col = c("blue", "red"), lwd = 2)

# plot cumulative distribution
df_cum <- rbind(
  data.frame(type = "Observed episodes", cumul = cum_obs),
  data.frame(type = "Simulated episodes", cumul = cum_sim)
)

ggplot(df_cum, aes(x = type, y = cumul, fill = type)) +
  geom_violin(alpha = 0.7, trim = FALSE) +
  labs(x = "", y = "Total cumulative rainfall (mm)") +
  theme_minimal() + btf_theme +
  theme(legend.position = "none")
