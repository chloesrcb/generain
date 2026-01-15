
################################################################################

# look at marginals
sims_all_X <- do.call(cbind, lapply(sims_group, function(sim) sim$X))

# crps
library(scoringRules)

# sims_group : list of simulations
# each sim$X : matrix (n_sites x lt) with rownames = sites

sites <- rownames(sims_group[[1]]$X)
lt <- ncol(sims_group[[1]]$X)
Nsim <- length(sims_group)

# Stack all values for each site across sims and times
# Result: matrix n_sites x (Nsim*lt)
X_by_site <- do.call(cbind, lapply(sims_group, function(sim) sim$X))
# X_by_site rows = sites, cols = concatenated times across sims
dim(X_by_site)


Hgpd <- function(x, sigma, xi) {
  x <- pmax(x, 0)
  if (abs(xi) < 1e-10) return(1 - exp(-x / sigma))
  1 - (1 + xi * x / sigma)^(-1/xi)
}

Fegpd_pos <- function(x, kappa, sigma, xi) {
  stopifnot(all(x > 0))
  Hgpd(x, sigma, xi)^kappa
}

Qgpd <- function(p, sigma, xi) {
  p <- pmin(pmax(p, 1e-12), 1 - 1e-12)
  if (abs(xi) < 1e-10) return(-sigma * log(1 - p))
  sigma/xi * ((1 - p)^(-xi) - 1)
}

r_egpd_pos <- function(n, kappa, sigma, xi) {
  u <- runif(n)
  Qgpd(u^(1 / kappa), sigma, xi)
}

crps_from_samples <- function(y, xs) {
  mean(abs(xs - y)) - 0.5 * mean(abs(outer(xs, xs, "-")))
}

crps_egpd_mc_pos <- function(y, kappa, sigma, xi, m = 500) {
  if (y <= 0) return(NA_real_)
  xs <- r_egpd_pos(m, kappa, sigma, xi)
  crps_from_samples(y, xs)
}


X_by_site <- do.call(cbind, lapply(sims_group, function(sim) sim$X))
sites <- rownames(X_by_site)


crps_site <- sapply(seq_along(sites), function(i) {

  x <- as.numeric(X_by_site[i, ])
  x <- x[is.finite(x) & x > 0] 

  if (length(x) < 50) return(NA_real_)

  mean(vapply(
    x,
    function(y) crps_egpd_mc_pos(
      y,
      kappa = params_margins$kappa[i],
      sigma = params_margins$sigma[i],
      xi    = params_margins$xi[i],
      m = 400
    ),
    numeric(1)
  ))
})

data.frame(site = sites, crps_egpd = crps_site)


library(scoringRules)

crps_site_emp <- sapply(seq_along(sites), function(i) {
  x <- as.numeric(X_by_site[i, ])
  x <- x[is.finite(x) & x > 0]
  if (length(x) < 50) return(NA_real_)
  mean(vapply(x, function(y) crps_sample(y = y, dat = x), numeric(1)))
})

res <- data.frame(
  site = sites,
  crps_egpd = crps_site,
  crps_emp  = crps_site_emp,
  ratio     = crps_site / crps_site_emp
)

res[order(res$ratio), ]

# save in csv
foldername_crps <- paste0(
  data_folder,
  "swg/omsev/crps/"
)
if (!dir.exists(foldername_crps)) {
  dir.create(foldername_crps, recursive = TRUE)
}
filename_crps <- paste0(
  foldername_crps,
  "crps_ratio_advgroup_",
  group_adv,
  "_n", Nsim,
  ".csv"
)

write.csv(res, filename_crps, row.names = FALSE)

pit_site_pos <- function(x, kappa, sigma, xi) {
  x <- x[x > 0]
  Fegpd_pos(x, kappa, sigma, xi)
}

# Exemple site i
i <- 1
x <- as.numeric(X_by_site[i, ])
u <- pit_site_pos(
  x,
  params_margins$kappa[i],
  params_margins$sigma[i],
  params_margins$xi[i]
)


res$group <- "All sites"

ggplot(res, aes(x = group, y = ratio)) +
  geom_boxplot(
    fill = "#69b3a2",
    alpha = 0.6,
    width = 0.3
  ) +
  labs(
    y = "CRPS(EGPD) / CRPS(empirical)",
    x = ""
  ) +
  ylim(0.9, 1.2) +
  theme_minimal()


# save plot
foldername_plot_crps <- paste0(
  im_folder,
  "swg/omsev/simu/crps/"
)

if (!dir.exists(foldername_plot_crps)) {
  dir.create(foldername_plot_crps, recursive = TRUE)
}

filename_plot_crps <- paste0(
  foldername_plot_crps,
  "crps_ratio_advgroup_",
  group_adv,
  "_n", Nsim,
  ".png"
)






compute_crps_marginal_group <- function(sims_group, params_margins) {

  # sims_group : list of simulations (un seul groupe d’advection)

  sites <- rownames(sims_group[[1]]$X)

  # Stack all values for each site across sims and times
  X_by_site <- do.call(cbind, lapply(sims_group, function(sim) sim$X))

  crps_site <- sapply(seq_along(sites), function(i) {

    x <- as.numeric(X_by_site[i, ])
    x <- x[is.finite(x) & x > 0]

    if (length(x) < 50) return(NA_real_)

    mean(vapply(
      x,
      function(y) crps_egpd_mc_pos(
        y,
        kappa = params_margins$kappa[i],
        sigma = params_margins$sigma[i],
        xi    = params_margins$xi[i],
        m = 400
      ),
      numeric(1)
    ))
  })

  crps_site_emp <- sapply(seq_along(sites), function(i) {

    x <- as.numeric(X_by_site[i, ])
    x <- x[is.finite(x) & x > 0]

    if (length(x) < 50) return(NA_real_)

    mean(vapply(
      x,
      function(y) crps_sample(y = y, dat = x),
      numeric(1)
    ))
  })

  data.frame(
    site = sites,
    crps_egpd = crps_site,
    crps_emp  = crps_site_emp,
    ratio     = crps_site / crps_site_emp
  )
}


results_all_groups <- list()

for (g in all_group_names) {

  cat("Processing group:", g, "\n")

  indices_group <- which(selected_points$adv_group == g)
  if (length(indices_group) < 10) next

  list_episodes_group <- list_episodes[indices_group]
  s0_list_group <- s0_list[indices_group]
  u_list_group  <- u_list[indices_group]
  adv_group     <- adv_matrix[indices_group, ]

  Nsim <- 100
  s0_sim <- integer(Nsim)
  sims_group <- vector("list", Nsim)

  for (i in seq_len(Nsim)) {
    idx <- sample(seq_along(list_episodes_group), 1)
    s0_i <- s0_list_group[idx]
    s0_sim[i] <- s0_i

    sims_group[[i]] <- simulate_many_episodes(
      N = 1,
      u = 1000,
      u_emp = u_list_group[idx],
      params_vario = params_vario_kmh,
      params_margins = params_margins,
      coords = grid_omsev,
      times = times * 5 / 60,
      adv = adv_group[idx, ],
      t0 = 0,
      s0 = s0_i
    )[[1]]
  }

  res_group <- compute_crps_marginal_group(
    sims_group,
    params_margins
  )

  res_group$group <- g
  results_all_groups[[g]] <- res_group
}


crps_all <- do.call(rbind, results_all_groups)
crps_all <- crps_all[is.finite(crps_all$ratio), ]

library(ggplot2)

ggplot(crps_all, aes(x = group, y = ratio)) +
  geom_boxplot(fill = "#69b3a2", alpha = 0.6) +
  coord_flip() +
  labs(
    x = "Advection group",
    y = "CRPS(EGPD) / CRPS(empirical)",
    title = ""
  ) +
  ylim(0.9, 1.3) +
  theme_minimal()

# save plot
foldername_plot_crps_groups <- paste0(
  im_folder,
  "swg/omsev/simu/crps/"
)

if (!dir.exists(foldername_plot_crps_groups)) {
  dir.create(foldername_plot_crps_groups, recursive = TRUE)
}

filename_plot_crps_groups <- paste0(
  foldername_plot_crps_groups,
  "crps_ratio_all_advgroups_n", Nsim,
  ".png"
)

ggsave(
  filename = filename_plot_crps_groups,
  width = 8,
  height = 6,
  dpi = 300
)