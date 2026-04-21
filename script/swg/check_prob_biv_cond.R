


#### 
# Conditional pairs probability plots
####


grid_omsev <- grid_coords_km
adv_matrix <- as.matrix(adv_df_transfo[, c("vx_t", "vy_t")])
all_group_names <- unique(selected_points$speed_class)
# all_significant_groups <- all_group_names[grep("significant", all_group_names)]
group_adv <- all_group_names[3]  # choose one group to simulate
# get list_episodes and s0_list for this group
indices_group <- which(selected_points$speed_class == group_adv)
list_episodes_group <- list_episodes[indices_group]
length(list_episodes_group)
s0_list_group <- s0_list[indices_group]
u_list_group <- u_list[indices_group]
adv_group <- adv_matrix[indices_group, ]
Nsim <- 100
s0_sim <- integer(Nsim)
u_sim <- numeric(Nsim)
sims_group <- vector("list", Nsim)
adv_sim <- matrix(0, nrow = Nsim, ncol = 2)
for (i in seq_len(Nsim)) {
  idx <- sample(seq_along(list_episodes_group), 1)
  s0_i <- s0_list_group[idx]
  s0_sim[i] <- s0_i
  u_i <- u_list_group[idx]
  u_sim[i] <- u_i
  adv_i <- adv_group[idx, ]
  adv_sim[i, ] <- adv_i
  sims_group[[i]] <- simulate_many_episodes(
    N = 1,
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


sim_ep <- sims_group[[1]]
X_ep <- sim_ep$X

site1 <- "iem"
site2 <- "mse"
site3 <- "cnrs"
u <- u_list_group[1]

I_s3    <- as.integer(X_ep[, site3] > u)
I_s1s2s3 <- as.integer(X_ep[, site1] > u & X_ep[, site2] > u & X_ep[, site3] > u)

prob_s1s2_given_s3 <- sum(I_s1s2s3) / sum(I_s3)
prob_s1s2_given_s3

ind_dep_s3 <- X_ep[, site3] > u
n_cond <- sum(ind_dep_s3)

prob_s1s2_given_s3 <- if (n_cond == 0) NA_real_ else
  mean(X_ep[ind_dep_s3, site1] > u & X_ep[ind_dep_s3, site2] > u)

prob_s1s2_given_s3


# Computes P(X_s1>u & X_s2>u | X_s3>u) for all sites in one sim
cond_probs_all_sites <- function(X, u) {
  sites <- colnames(X)
  p <- length(sites)
  n <- nrow(X)

  # exceedances matrix: n x p (logical)
  E <- X > u

  # denominators for each s3: count(X_s3 > u)
  denom <- colSums(E)  # length p
  # if there is NA put denom to 0
  denom[is.na(denom)] <- 0
  # allocate result: s1 x s2 x s3
  P <- array(NA, dim = c(p, p, p),
             dimnames = list(s1 = sites, s2 = sites, s3 = sites))

  # loop over conditioning site s3, compute pairwise probs among exceedances of s3
  for (k in seq_len(p)) {
    nk <- denom[k] # number of exceedances at site k
    if (nk == 0) next  # stays NA

    Ek <- E[E[, k], , drop = FALSE] # rows where s3 exceeds

    P[, , k] <- crossprod(Ek) / nk
  }

  P
}


probs_by_group_obs <- lapply(seq_along(list_episodes_group), function(i) {
  X <- list_episodes_group[[i]]
  u <- u_list_group[i]
  cond_probs_all_sites(X, u)
})

names(probs_by_group_obs) <- paste0("group_", seq_along(probs_by_group_obs))

# Example: probability for group 1, (iem,mse | cnrs)
probs_by_group_obs[["group_1"]]["iem", "mse", "cnrs"]

# Average over all groups without NA
p <- dim(probs_by_group_obs[[1]])[1]
cond_probs_avg_obs <- array(NA_real_, dim = c(p, p, p),
                        dimnames = list(s1 = sites_names, s2 = sites_names, s3 = sites_names))
for (i in seq_len(p)) {
  for (j in seq_len(p)) {
    for (k in seq_len(p)) { 
      vals_ijk <- sapply(probs_by_group_obs, function(P) P[i, j, k])
      cond_probs_avg_obs[i, j, k] <- mean(vals_ijk, na.rm = TRUE)
    }
  }
}



# Apply to all sims in sims_group
# u_list_group[g] assumed to match sims_group[[g]]
probs_by_group_sim <- lapply(seq_along(sims_group), function(g) {
  X <- sims_group[[g]]$X
  u <- u_sim[g]
  cond_probs_all_sites(X, u)
})

names(probs_by_group_sim) <- paste0("group_", seq_along(probs_by_group_sim))

# Example: probability for group 1, (iem,mse | cnrs)
probs_by_group_sim[["group_1"]]["iem", "mse", "cnrs"]

# Average over all groups without NA
p <- dim(probs_by_group_sim[[1]])[1]
cond_probs_avg_sim <- array(NA_real_, dim = c(p, p, p),
                        dimnames = list(s1 = sites_names, s2 = sites_names, s3 = sites_names))
for (i in seq_len(p)) {
  for (j in seq_len(p)) {
    for (k in seq_len(p)) { 
      vals_ijk <- sapply(probs_by_group_sim, function(P) P[i, j, k])
      cond_probs_avg_sim[i, j, k] <- mean(vals_ijk, na.rm = TRUE)
    }
  }
}



################################################################################
# Conditional exceedance probabilities: plots for all groups and average
################################################################################

cond_probs_to_df <- function(P, group = NA_character_) {
  dn <- dimnames(P)
  if (is.null(dn)) {
    stop("cond_probs_to_df: array must have dimnames for s1, s2, s3")
  }
  df <- expand.grid(
    s1 = dn$s1,
    s2 = dn$s2,
    s3 = dn$s3,
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )
  df$prob <- as.vector(P)
  if (!is.na(group)) df$group <- group
  df
}

# Long data for all groups
df_cond_probs_all_sim <- do.call(
  rbind,
  lapply(names(probs_by_group_sim), function(g) {
    cond_probs_to_df(probs_by_group_sim[[g]], g)
  })
)

df_cond_probs_all_obs <- do.call(
  rbind,
  lapply(names(probs_by_group_obs), function(g) {
    cond_probs_to_df(probs_by_group_obs[[g]], g)
  })
)

# Plot: all groups x all conditioning sites

foldername_plot <- file.path(im_folder, "swg/conditional_probs")
if (!dir.exists(foldername_plot)) {
  dir.create(foldername_plot, recursive = TRUE)
}

# Plot: average over groups (all conditioning sites)
df_cond_probs_avg_sim <- cond_probs_to_df(cond_probs_avg_sim)
df_cond_probs_avg_obs <- cond_probs_to_df(cond_probs_avg_obs)

ggplot(
  df_cond_probs_avg_sim,
  aes(x = s1, y = s2, fill = prob)
) +
  geom_tile(color = NA) +
  facet_wrap(~ s3) +
  scale_fill_gradient(
    low = "white",
    high = "steelblue",
    limits = c(0, 1),
    na.value = "grey90"
  ) +
  coord_equal() +
  labs(
    x = "Site 1 (s1)",
    y = "Site 2 (s2)",
    fill = "Mean P(X1>u, X2>u | X3>u)"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    panel.grid = element_blank()
  )

filename_plot <- paste0(
  foldername_plot,
  "cond_probs_avg_sim.png"
)
ggsave(
  filename = filename_plot,
  plot = last_plot(),
  width = 12,
  height = 8,
  dpi = 300
)

ggplot(
  df_cond_probs_avg_obs,
  aes(x = s1, y = s2, fill = prob)
) +
  geom_tile(color = NA) +
  facet_wrap(~ s3) +
  scale_fill_gradient(
    low = "white",
    high = "steelblue",
    limits = c(0, 1),
    na.value = "grey90"
  ) +
  coord_equal() +
  labs(
    x = "Site 1 (s1)",
    y = "Site 2 (s2)",
    fill = "Mean P(X1>u, X2>u | X3>u)"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    panel.grid = element_blank()
  )



filename_plot <- paste0(
  foldername_plot,
  "cond_probs_avg_obs.png"
)
ggsave(
  filename = filename_plot,
  plot = last_plot(),
  width = 12,
  height = 8,
  dpi = 300
)




# moyenne pondérée sur une liste de P (array p x p x p) en pondérant par denom(s3)
weighted_avg_cond_probs <- function(P_list, denom_list) {
  # P_list: liste de arrays [s1,s2,s3]
  # denom_list: liste de vecteurs denom par épisode (longueur p), denom[k]=#(Xs3>u)

  p <- dim(P_list[[1]])[1]
  sites <- dimnames(P_list[[1]])$s1

  num <- array(0, dim = c(p,p,p), dimnames = list(s1=sites,s2=sites,s3=sites))
  den <- array(0, dim = c(1,1,p), dimnames = list(NULL,NULL,sites))

  for (e in seq_along(P_list)) {
    P <- P_list[[e]]
    d <- denom_list[[e]] # longueur p

    for (k in seq_len(p)) {
      if (!is.finite(d[k]) || d[k] <= 0) next
      # P[,,k] = crossprod(Ek)/d[k] donc crossprod(Ek) = P* d[k]
      num[,,k] <- num[,,k] + P[,,k] * d[k]
      den[1,1,k] <- den[1,1,k] + d[k]
    }
  }

  out <- array(NA_real_, dim = c(p,p,p), dimnames = list(s1=sites,s2=sites,s3=sites))
  for (k in seq_len(p)) {
    if (den[1,1,k] > 0) out[,,k] <- num[,,k] / den[1,1,k]
  }
  out
}

cond_probs_all_sites_with_denom <- function(X, u) {
  sites <- colnames(X)
  p <- length(sites)

  E <- (X > u)
  # si NA -> FALSE (plus propre)
  E[is.na(E)] <- FALSE

  denom <- colSums(E)  # #(Xs>u) par site

  P <- array(NA_real_, dim = c(p, p, p),
             dimnames = list(s1 = sites, s2 = sites, s3 = sites))

  for (k in seq_len(p)) {
    nk <- denom[k]
    if (nk == 0) next
    Ek <- E[E[, k], , drop = FALSE]
    P[, , k] <- crossprod(Ek) / nk
  }

  list(P = P, denom = denom)
}


library(ggplot2)

cond_probs_to_df <- function(P, group = NA_character_, source = NA_character_) {
  dn <- dimnames(P)
  df <- expand.grid(s1 = dn$s1, s2 = dn$s2, s3 = dn$s3,
                    KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
  df$prob <- as.vector(P)
  if (!is.na(group)) df$group <- group
  if (!is.na(source)) df$source <- source
  df
}

make_heatmap <- function(df, title = NULL) {
  ggplot(df, aes(x = s1, y = s2, fill = prob)) +
    geom_tile(color = NA) +
    facet_wrap(~ s3) +
    scale_fill_gradient(low = "white", high = "steelblue",
                        limits = c(0,1), na.value = "grey90") +
    coord_equal() +
    labs(x = "Site 1 (s1)", y = "Site 2 (s2)", fill = "P(Xs1>u, Xs2>u | Xs3>u)", title = title) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
          panel.grid = element_blank())
}

results_by_adv <- vector("list", length(all_group_names))
names(results_by_adv) <- all_group_names

if (!dir.exists(foldername_plot)) {
  dir.create(foldername_plot, recursive = TRUE)
}
for (gname in all_group_names) {

  indices_group <- which(selected_points$speed_class == gname)
  if (length(indices_group) == 0) next

  list_episodes_group <- list_episodes[indices_group]
  s0_list_group <- s0_list[indices_group]
  u_list_group  <- u_list[indices_group]
  adv_group     <- adv_matrix[indices_group, , drop = FALSE]

  obs_tmp <- lapply(seq_along(list_episodes_group), function(i) {
    cond_probs_all_sites_with_denom(list_episodes_group[[i]], u_list_group[i])
  })
  P_obs_list <- lapply(obs_tmp, `[[`, "P")
  d_obs_list <- lapply(obs_tmp, `[[`, "denom")

  P_obs_avg <- weighted_avg_cond_probs(P_obs_list, d_obs_list)

  Nsim <- 100
  sims_group <- vector("list", Nsim)
  u_sim <- numeric(Nsim)

  for (i in seq_len(Nsim)) {
    idx <- sample(seq_along(list_episodes_group), 1)
    u_i <- u_list_group[idx]
    u_sim[i] <- u_i

    sims_group[[i]] <- simulate_many_episodes(
      N = 1,
      u_emp = u_i,
      params_vario = params_kmh,
      params_margins = params_margins,
      coords = grid_omsev,
      times = times * 5 / 60,
      adv = adv_group[idx, ],
      t0 = 0,
      s0 = s0_list_group[idx]
    )[[1]]
  }

  sim_tmp <- lapply(seq_along(sims_group), function(i) {
    cond_probs_all_sites_with_denom(sims_group[[i]]$X, u_sim[i])
  })
  P_sim_list <- lapply(sim_tmp, `[[`, "P")
  d_sim_list <- lapply(sim_tmp, `[[`, "denom")

  P_sim_avg <- weighted_avg_cond_probs(P_sim_list, d_sim_list)

  df_obs <- cond_probs_to_df(P_obs_avg, group = gname, source = "obs")
  df_sim <- cond_probs_to_df(P_sim_avg, group = gname, source = "sim")

  p_obs <- make_heatmap(df_obs, title = paste0("Obs — ", gname))
  p_sim <- make_heatmap(df_sim, title = paste0("Sim — ", gname))

  ggsave(file.path(foldername_plot, paste0("cond_probs_obs_", gname, ".png")),
         p_obs, width = 12, height = 8, dpi = 300)
  ggsave(file.path(foldername_plot, paste0("cond_probs_sim_", gname, ".png")),
         p_sim, width = 12, height = 8, dpi = 300)

  results_by_adv[[gname]] <- list(P_obs_avg = P_obs_avg, P_sim_avg = P_sim_avg)
}


df_both <- rbind(df_obs, df_sim)

ggplot(df_both, aes(x = s1, y = s2, fill = prob)) +
  geom_tile() +
  facet_grid(s3 ~ source) +
  coord_equal() +
  scale_fill_gradient(low="white", high="steelblue", limits=c(0,1), na.value="grey90") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle=90, vjust=0.5, hjust=1),
        panel.grid = element_blank())



library(ggplot2)

# ---- Choix (à adapter)
s1_pick <- "iem"
s2_pick <- "mse"
s3_pick <- "archie"

# grille de p et axe
p_grid <- seq(0.85, 0.995, by = 0.005)
x_grid <- -log(1 - p_grid)

# data.frame résultat
df <- data.frame(
  x = rep(x_grid, times = 2),
  p = rep(p_grid, times = 2),
  prob = NA_real_,
  source = rep(c("obs", "sim"), each = length(p_grid))
)

num_obs <- numeric(length(p_grid))
den_obs <- numeric(length(p_grid))

for (e in seq_along(list_episodes_group)) {
  X <- list_episodes_group[[e]]

  # seuils z(p) basés sur la colonne s3
  z_grid <- as.numeric(quantile(X[, s3_pick], probs = p_grid, na.rm = TRUE, type = 8))

  for (t in seq_along(p_grid)) {
    z <- z_grid[t]

    ind_s3 <- X[, s3_pick] > z
    ind_s3[is.na(ind_s3)] <- FALSE
    nk <- sum(ind_s3)
    if (nk == 0) next

    den_obs[t] <- den_obs[t] + nk

    joint <- (X[ind_s3, s1_pick] > z) & (X[ind_s3, s2_pick] > z)
    joint[is.na(joint)] <- FALSE

    num_obs[t] <- num_obs[t] + sum(joint)
  }
}

prob_obs <- ifelse(den_obs > 0, num_obs / den_obs, NA_real_)
df$prob[df$source == "obs"] <- prob_obs

# =========================
# 2) SIM : pool sur sims_group
# =========================
num_sim <- numeric(length(p_grid))
den_sim <- numeric(length(p_grid))

for (e in seq_along(sims_group)) {
  X <- sims_group[[e]]$X

  z_grid <- as.numeric(quantile(X[, s3_pick], probs = p_grid, na.rm = TRUE, type = 8))

  for (t in seq_along(p_grid)) {
    z <- z_grid[t]

    ind_s3 <- X[, s3_pick] > z
    ind_s3[is.na(ind_s3)] <- FALSE
    nk <- sum(ind_s3)
    if (nk == 0) next

    den_sim[t] <- den_sim[t] + nk

    joint <- (X[ind_s3, s1_pick] > z) & (X[ind_s3, s2_pick] > z)
    joint[is.na(joint)] <- FALSE

    num_sim[t] <- num_sim[t] + sum(joint)
  }
}

prob_sim <- ifelse(den_sim > 0, num_sim / den_sim, NA_real_)
df$prob[df$source == "sim"] <- prob_sim

# ---- sécurité anti-df vide
df_plot <- subset(df, is.finite(prob))

if (nrow(df_plot) == 0) {
  stop("df_plot est vide: aucun dépassement pour ces sites/niveaux p. Essaye de baisser p_grid ou changer s3_pick.")
}

# =========================
# 3) Plot (SANS facet)
# =========================
p <- ggplot(df_plot, aes(x = x, y = prob, linetype = source)) +
  geom_line(linewidth = 0.9) +
  geom_point(size = 1.5) +
  labs(
    x = "-log(1-p)",
    y = paste0("P(", s1_pick, ">z & ", s2_pick, ">z | ", s3_pick, ">z)"),
    title = paste0("Conditional joint exceedance — (", s1_pick, ", ", s2_pick, " | ", s3_pick, ")")
  ) +
  theme_minimal() +
  theme(panel.grid = element_blank())

print(p)

ggsave(
  filename = file.path(foldername_plot, paste0("cond_joint_curve_", s1_pick, "_", s2_pick, "_given_", s3_pick, ".png")),
  plot = p,
  width = 10, height = 6, dpi = 300
)




library(ggplot2)

s1_pick <- "cnrs"
s2_pick <- "mse"
s3_pick <- "archie"

# grille de p (si trop haut -> aucune valeur, baisse un peu)
p_grid <- seq(0.80, 0.99, by = 0.01)
x_grid <- -log(1 - p_grid)

# tous les groupes d’advection
all_group_names <- unique(selected_points$speed_class)

for (gname in all_group_names) {

  indices_group <- which(selected_points$speed_class == gname)
  if (length(indices_group) == 0) next

  list_episodes_group <- list_episodes[indices_group]
  s0_list_group <- s0_list[indices_group]
  u_list_group  <- u_list[indices_group]
  adv_group     <- adv_matrix[indices_group, , drop = FALSE]

  # -----------------------
  # OBS (pool sur épisodes)
  # -----------------------
  num_obs <- numeric(length(p_grid))
  den_obs <- numeric(length(p_grid))

  for (e in seq_along(list_episodes_group)) {
    X <- list_episodes_group[[e]]

    # si jamais les noms changent:
    if (!(s1_pick %in% colnames(X) && s2_pick %in% colnames(X) && s3_pick %in% colnames(X))) next

    z_grid <- as.numeric(quantile(X[, s3_pick], probs = p_grid, na.rm = TRUE, type = 8))

    for (t in seq_along(p_grid)) {
      z <- z_grid[t]

      ind_s3 <- X[, s3_pick] > z
      ind_s3[is.na(ind_s3)] <- FALSE
      nk <- sum(ind_s3)
      if (nk == 0) next

      den_obs[t] <- den_obs[t] + nk

      joint <- (X[ind_s3, s1_pick] > z) & (X[ind_s3, s2_pick] > z)
      joint[is.na(joint)] <- FALSE
      num_obs[t] <- num_obs[t] + sum(joint)
    }
  }

  prob_obs <- ifelse(den_obs > 0, num_obs / den_obs, NA_real_)

  # -----------------------
  # SIM (re-simuler Nsim)
  # -----------------------
  Nsim <- 100
  sims_group <- vector("list", Nsim)

  for (i in seq_len(Nsim)) {
    idx <- sample(seq_along(list_episodes_group), 1)

    sims_group[[i]] <- simulate_many_episodes(
      N = 1,
      u_emp = u_list_group[idx],
      params_vario = params_kmh,
      params_margins = params_margins,
      coords = grid_omsev,
      times = times * 5 / 60,
      adv = adv_group[idx, ],
      t0 = 0,
      s0 = s0_list_group[idx]
    )[[1]]
  }

  num_sim <- numeric(length(p_grid))
  den_sim <- numeric(length(p_grid))

  for (e in seq_along(sims_group)) {
    X <- sims_group[[e]]$X
    if (!(s1_pick %in% colnames(X) && s2_pick %in% colnames(X) && s3_pick %in% colnames(X))) next

    z_grid <- as.numeric(quantile(X[, s3_pick], probs = p_grid, na.rm = TRUE, type = 8))

    for (t in seq_along(p_grid)) {
      z <- z_grid[t]

      ind_s3 <- X[, s3_pick] > z
      ind_s3[is.na(ind_s3)] <- FALSE
      nk <- sum(ind_s3)
      if (nk == 0) next

      den_sim[t] <- den_sim[t] + nk

      joint <- (X[ind_s3, s1_pick] > z) & (X[ind_s3, s2_pick] > z)
      joint[is.na(joint)] <- FALSE
      num_sim[t] <- num_sim[t] + sum(joint)
    }
  }

  prob_sim <- ifelse(den_sim > 0, num_sim / den_sim, NA_real_)

  # -----------------------
  # Plot
  # -----------------------
  df <- rbind(
    data.frame(x = x_grid, p = p_grid, prob = prob_obs, source = "obs"),
    data.frame(x = x_grid, p = p_grid, prob = prob_sim, source = "sim")
  )
  df <- subset(df, is.finite(prob))

  # si c'est vide pour ce groupe, on saute
  if (nrow(df) == 0) next

  p <- ggplot(df, aes(x = x, y = prob, linetype = source)) +
    geom_line(linewidth = 0.9) +
    geom_point(size = 1.6) +
    labs(
      x = "-log(1-p)",
      y = paste0("P(", s1_pick, ">u & ", s2_pick, ">u | ", s3_pick, ">u)"),
      title = paste0("Advection group = ", gname, "  —  (", s1_pick, ", ", s2_pick, " | ", s3_pick, ")")
    ) +
    theme_minimal() +
    theme(panel.grid = element_blank())

  print(p)

  ggsave(
    filename = file.path(foldername_plot,
                         paste0("curve_", gname, "_", s1_pick, "_", s2_pick, "_given_", s3_pick, ".png")),
    plot = p, width = 10, height = 6, dpi = 300
  )
}













library(ggplot2)

# ---- choix : site conditionnant
s3_pick <- "archie"  # à adapter

# ---- grille p
p_grid <- seq(0.80, 0.99, by = 0.01)
x_grid <- -log(1 - p_grid)

# ---- tous les groupes d’advection
all_group_names <- unique(selected_points$speed_class)

for (gname in all_group_names) {

  indices_group <- which(selected_points$speed_class == gname)
  if (length(indices_group) == 0) next

  list_episodes_group <- list_episodes[indices_group]
  s0_list_group <- s0_list[indices_group]
  u_list_group  <- u_list[indices_group]
  adv_group     <- adv_matrix[indices_group, , drop = FALSE]

  # sites
  sites <- colnames(list_episodes_group[[1]])
  if (!(s3_pick %in% sites)) next
  Psites <- length(sites)
  k <- which(sites == s3_pick)

  # on veut toutes les paires (i<j) hors s3
  idx_other <- setdiff(seq_len(Psites), k)
  pairs <- t(combn(idx_other, 2))  # matrice nPairs x 2
  nPairs <- nrow(pairs)

  # =========================
  # OBS: pool sur épisodes
  # =========================
  meanpair_obs <- rep(NA_real_, length(p_grid))

  for (t in seq_along(p_grid)) {
    num_sum <- 0
    den_sum <- 0

    for (e in seq_along(list_episodes_group)) {
      X <- list_episodes_group[[e]]
      if (!(s3_pick %in% colnames(X))) next

      z <- as.numeric(quantile(X[, s3_pick], probs = p_grid[t], na.rm = TRUE, type = 8))

      ind_s3 <- X[, s3_pick] > z
      ind_s3[is.na(ind_s3)] <- FALSE
      nk <- sum(ind_s3)
      if (nk == 0) next

      # sous-échantillon conditionnel
      Xc <- X[ind_s3, , drop = FALSE]

      # exceedances bool (NA -> FALSE)
      E <- (Xc > z)
      E[is.na(E)] <- FALSE

      # counts de co-dépassements pour toutes les paires
      C <- crossprod(E)  # Psites x Psites, comptes sur nk jours

      # moyenne sur les paires (i<j) hors s3
      cvals <- C[pairs]  # vector de longueur nPairs (indexation matricielle)
      num_sum <- num_sum + sum(cvals)
      den_sum <- den_sum + nk * nPairs
    }

    meanpair_obs[t] <- if (den_sum > 0) num_sum / den_sum else NA_real_
  }

  Nsim <- 100
  sims_group <- vector("list", Nsim)
  s0_sim <- character(Nsim)
  s0_possible <- sites[-which(sites == s3_pick)]
  # get idx possible for s0 in s0_list_group
  s0_indices <- which(s0_list_group %in% s0_possible)
  for (i in seq_len(Nsim)) {
    idx <- sample(s0_indices, 1)

    # stocker s0 choisi
    s0_sim[i] <- s0_list_group[idx]

    sims_group[[i]] <- simulate_many_episodes(
      N = 1,
      u_emp = u_list_group[idx],
      params_vario = params_kmh,
      params_margins = params_margins,
      coords = grid_omsev,
      times = times * 5 / 60,
      adv = adv_group[idx, ],
      t0 = 0,
      s0 = s0_list_group[idx]
    )[[1]]
  }


  meanpair_sim <- rep(NA_real_, length(p_grid))

  for (t in seq_along(p_grid)) {
    num_sum <- 0
    den_sum <- 0

    for (e in seq_along(sims_group)) {
      X <- sims_group[[e]]$X
      if (!(s3_pick %in% colnames(X))) next

      z <- as.numeric(quantile(X[, s3_pick], probs = p_grid[t], na.rm = TRUE, type = 8))

      ind_s3 <- X[, s3_pick] > z
      ind_s3[is.na(ind_s3)] <- FALSE
      nk <- sum(ind_s3)
      if (nk == 0) next

      Xc <- X[ind_s3, , drop = FALSE]
      E <- (Xc > z)
      E[is.na(E)] <- FALSE

      C <- crossprod(E)
      cvals <- C[pairs]
      num_sum <- num_sum + sum(cvals)
      den_sum <- den_sum + nk * nPairs
    }

    meanpair_sim[t] <- if (den_sum > 0) num_sum / den_sum else NA_real_
  }

  # =========================
  # Plot par groupe
  # =========================
  df <- rbind(
    data.frame(x = x_grid, p = p_grid, prob = meanpair_obs, source = "obs"),
    data.frame(x = x_grid, p = p_grid, prob = meanpair_sim, source = "sim")
  )
  df <- subset(df, is.finite(prob))
  if (nrow(df) == 0) next

  p <- ggplot(df, aes(x = x, y = prob, linetype = source)) +
    geom_line(linewidth = 0.9) +
    geom_point(size = 1.6) +
    labs(
      x = "-log(1-p)",
      y = paste0("Mean over pairs of P(Xs1>z & Xs2>z | ", s3_pick, ">z)"),
      title = paste0("Advection group = ", gname, " — conditioning site ", s3_pick)
    ) +
    theme_minimal() +
    theme(panel.grid = element_blank())

  print(p)

  ggsave(
    filename = file.path(foldername_plot, paste0("mean_pairs_curve_", gname, "_given_", s3_pick, ".png")),
    plot = p, width = 10, height = 6, dpi = 300
  )
}








library(ggplot2)

s1_pick <- "iem"
s2_pick <- "mse"
s3_pick <- "cnrs"

p_grid <- seq(0.80, 0.99, by = 0.01)
x_grid <- -log(1 - p_grid)

# bootstrap
B <- 300
L <- 10

num_obs <- numeric(length(p_grid))
den_obs <- numeric(length(p_grid))

for (e in seq_along(list_episodes_group)) {
  X <- list_episodes_group[[e]]
  z_grid <- as.numeric(quantile(X[, s3_pick], probs = p_grid, na.rm = TRUE, type = 8))

  for (t in seq_along(p_grid)) {
    z <- z_grid[t]

    ind_s3 <- X[, s3_pick] > z
    ind_s3[is.na(ind_s3)] <- FALSE
    nk <- sum(ind_s3)
    if (nk == 0) next

    den_obs[t] <- den_obs[t] + nk

    joint <- (X[ind_s3, s1_pick] > z) & (X[ind_s3, s2_pick] > z)
    joint[is.na(joint)] <- FALSE
    num_obs[t] <- num_obs[t] + sum(joint)
  }
}

prob_obs <- ifelse(den_obs > 0, num_obs / den_obs, NA)

keep <- (s0_sim != s3_pick) & !is.na(s0_sim)
sims_use <- sims_group[keep]
# sims_use <- sims_group

num_sim <- numeric(length(p_grid))
den_sim <- numeric(length(p_grid))

for (e in seq_along(sims_use)) {
  X <- sims_use[[e]]$X
  if (!(s1_pick %in% colnames(X) && s2_pick %in% colnames(X) && s3_pick %in% colnames(X))) next

  z_grid <- as.numeric(quantile(X[, s3_pick], probs = p_grid, na.rm = TRUE, type = 8))

  for (t in seq_along(p_grid)) {
    z <- z_grid[t]
    ind_s3 <- X[, s3_pick] > z
    ind_s3[is.na(ind_s3)] <- FALSE
    nk <- sum(ind_s3)
    if (nk == 0) next

    den_sim[t] <- den_sim[t] + nk

    joint <- (X[ind_s3, s1_pick] > z) & (X[ind_s3, s2_pick] > z)
    joint[is.na(joint)] <- FALSE
    num_sim[t] <- num_sim[t] + sum(joint)
  }
}

prob_sim <- ifelse(den_sim > 0, num_sim / den_sim, NA)


boot_mat <- matrix(NA, nrow = B, ncol = length(p_grid))

for (b in seq_len(B)) {

  num_b <- numeric(length(p_grid))
  den_b <- numeric(length(p_grid))

  for (e in seq_along(sims_use)) {
    X0 <- sims_use[[e]]$X
    n <- nrow(X0)
    if (n < L) next
    if (!(s1_pick %in% colnames(X0) && s2_pick %in% colnames(X0) && s3_pick %in% colnames(X0))) next

    starts <- 1:(n - L + 1)
    idx <- integer(0)
    while (length(idx) < n) {
      s <- sample(starts, 1)
      idx <- c(idx, s:(s + L - 1))
    }
    idx <- idx[1:n]
    X <- X0[idx, , drop = FALSE]

    z_grid <- as.numeric(quantile(X[, s3_pick], probs = p_grid, na.rm = TRUE, type = 8))

    for (t in seq_along(p_grid)) {
      z <- z_grid[t]
      ind_s3 <- X[, s3_pick] > z
      ind_s3[is.na(ind_s3)] <- FALSE
      nk <- sum(ind_s3)
      if (nk == 0) next

      den_b[t] <- den_b[t] + nk

      joint <- (X[ind_s3, s1_pick] > z) & (X[ind_s3, s2_pick] > z)
      joint[is.na(joint)] <- FALSE
      num_b[t] <- num_b[t] + sum(joint)
    }
  }

  boot_mat[b, ] <- ifelse(den_b > 0, num_b / den_b, NA)
}

ci_lo <- apply(boot_mat, 2, quantile, probs = 0.025, na.rm = TRUE)
ci_hi <- apply(boot_mat, 2, quantile, probs = 0.975, na.rm = TRUE)

df_obs <- data.frame(x = x_grid, p = p_grid, prob = prob_obs, source = "obs")
df_sim <- data.frame(x = x_grid, p = p_grid, prob = prob_sim, source = "sim")
df_ci  <- data.frame(x = x_grid, p = p_grid, lo = ci_lo, hi = ci_hi)

df_obs <- subset(df_obs, is.finite(prob))
df_sim <- subset(df_sim, is.finite(prob))
df_ci  <- subset(df_ci, is.finite(lo) & is.finite(hi))

p <- ggplot() +
  geom_ribbon(data = df_ci, aes(x = x, ymin = lo, ymax = hi), alpha = 0.2, fill = "red") +
  geom_line(data = df_sim, aes(x = x, y = prob, linetype = source), linewidth = 0.9, color = "red") +
  geom_point(data = df_sim, aes(x = x, y = prob), size = 1.6, color = "red") +
  geom_line(data = df_obs, aes(x = x, y = prob, linetype = source), linewidth = 0.9) +
  geom_point(data = df_obs, aes(x = x, y = prob), size = 1.6) +
  scale_linetype_manual(values = c(obs = "solid", sim = "dashed")) +
  labs(
    x = "-log(1-p)",
    y = paste0("P(", s1_pick, ">z & ", s2_pick, ">z | ", s3_pick, ">z)"),
    title = "Conditional joint exceedance (sim band = block bootstrap 95%)",
    linetype = "source"
  ) +
  theme_minimal() +
  theme(panel.grid = element_blank())

print(p)
