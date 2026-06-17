
################################################################################
# Simulations over all episodes
################################################################################
adv_matrix <- as.matrix(adv_df_transfo[, c("vx_t", "vy_t")])
summary(adv_matrix)
speed <- sqrt(adv_matrix[, 1]^2 + adv_matrix[, 2]^2)
summary(speed)
speed[speed>150]
grid_omsev <- grid_coords_m

# convert advection from km/h to m/5min
adv_matrix <- adv_matrix * (1000 * 5 / 60)

# remove "invalid"
id_not_invalid <- which(adv_class$group != "invalid")
adv_matrix_filtered <- adv_matrix[id_not_invalid, , drop = FALSE]
list_episodes_filtered <- list_episodes[id_not_invalid]
s0_list_filtered <- s0_list[id_not_invalid]
u_list_filtered <- u_list[id_not_invalid]

M <- 100 # repetitions per episode
Nsim <- length(list_episodes_filtered)  # number of episodes to simulate
sims_by_ep <- vector("list", Nsim)
u_sim <- numeric(Nsim)
timesteps <- 0:11 #* 5 / 60 # every 5 minutes for 1 hour
set.seed(123)
for (j in seq_len(Nsim)) {
  ep_idx <- j
  s0_j  <- s0_list_filtered[ep_idx]
  u_j   <- u_list_filtered[ep_idx]
  adv_j <- adv_matrix_filtered[ep_idx, ]
  u_sim[j] <- u_j
  sims_by_ep[[j]] <- vector("list", M)

  for (m in seq_len(M)) {
  sims_by_ep[[j]][[m]] <- simulate_many_episodes(
      N = 1,
      u_emp = u_j,
      params_vario = params_m5min,
      params_margins = params_margins,
      coords = grid_omsev,
      times = timesteps,
      adv = adv_j,
      t0 = 0,
      s0 = s0_j
    )[[1]]
  }
}

apply_obs_na_mask <- function(Xsim, Xobs, sites) {
  Xsim2 <- Xsim
  mask_na <- is.na(as.matrix(Xobs)[, sites, drop=FALSE])
  Xsim2[mask_na] <- NA
  Xsim2
}

sites <- colnames(list_episodes_filtered[[1]])
sites <- colnames(list_episodes_filtered[[1]])
Xobs_all <- do.call(rbind, lapply(list_episodes_filtered, function(X) as.matrix(X)[, sites, drop=FALSE]))
u_by_s0_from_list <- tapply(u_list, s0_list, median, na.rm = TRUE)

u_by_s0_from_list

u_global <- u_by_s0_from_list
u_global <- as.vector(u_by_s0_from_list[sites])


sites <- colnames(list_episodes_filtered[[1]])

get_X_matrix <- function(obj) {
  if (is.null(obj)) return(NULL)
  if (is.matrix(obj) || is.data.frame(obj)) return(as.matrix(obj))
  if (is.list(obj) && !is.null(obj$X)) return(as.matrix(obj$X))
  NULL
}

sim_obj_list <- unlist(sims_by_ep, recursive = FALSE)
Xsim_all <- do.call(rbind, lapply(sim_obj_list, function(s) {
  X <- get_X_matrix(s)
  if (is.null(X)) return(NULL)
  as.matrix(X)[, sites, drop=FALSE]
}))


cond_probs_episode <- function(X, u, sites = colnames(X)) {
  X <- as.matrix(X)[, sites, drop = FALSE]
  p <- length(sites)

  if (length(u) == 1) u <- rep(u, p)
  if (is.null(names(u))) names(u) <- sites
  u <- u[sites]

  P <- array(NA_real_, dim = c(p,p,p), dimnames=list(s1=sites,s2=sites,s3=sites))
  denom <- rep(0, p); names(denom) <- sites

  for (k in seq_len(p)) {
    ok_k <- !is.na(X[,k])
    Ek_cond <- ok_k & (X[,k] > u[k])
    nk <- sum(Ek_cond)
    denom[k] <- nk
    if (nk == 0) next

    Xk <- X[Ek_cond, , drop=FALSE]
    E <- (Xk > matrix(u, nrow=nrow(Xk), ncol=p, byrow=TRUE))
    E[is.na(E)] <- FALSE

    P[, , k] <- crossprod(E) / nk 

  }

  list(P=P, denom=denom)
}



weighted_avg_cond_probs <- function(res_list) {
  p <- dim(res_list[[1]]$P)[1]
  sites <- dimnames(res_list[[1]]$P)$s1

  num <- array(0, dim = c(p, p, p), dimnames = list(s1=sites, s2=sites, s3=sites))
  den <- rep(0, p); names(den) <- sites

  for (r in res_list) {
    P <- r$P
    d <- r$denom
    for (k in seq_len(p)) {
      if (!is.finite(d[k]) || d[k] <= 0) next
      num[, , k] <- num[, , k] + P[, , k] * d[k]
      den[k] <- den[k] + d[k]
    }
  }

  out <- array(NA_real_, dim = c(p, p, p), dimnames = list(s1=sites, s2=sites, s3=sites))
  for (k in seq_len(p)) {
    if (den[k] > 0) out[, , k] <- num[, , k] / den[k]
  }
  out
}


selected_points_filtered<- selected_points[id_not_invalid, , drop = FALSE]
grp <- selected_points_filtered$speed_class
stopifnot(length(grp) == length(list_episodes_filtered))

P_obs_by_grp <- list()
P_sim_by_grp <- list()

for (g in levels(grp)) {
  idx <- which(grp == g)
  if (length(idx) == 0) next

  obs_res_g <- lapply(idx, function(j) {
    cond_probs_episode(list_episodes_filtered[[j]], u_list_filtered[j], sites)
  })
  P_obs_by_grp[[g]] <- weighted_avg_cond_probs(obs_res_g)

  sim_res_g <- do.call(c, lapply(idx, function(j) {
    u <- u_list_filtered[j]
    lapply(seq_along(sims_by_ep[[j]]), function(m) {
      cond_probs_episode(sims_by_ep[[j]][[m]]$X, u, sites)
    })
  }))
  P_sim_by_grp[[g]] <- weighted_avg_cond_probs(sim_res_g)
}


library(ggplot2)

s1_pick <- "iem"
s2_pick <- "mse"
s3_pick <- "cnrs"

p_grid <- seq(0.80, 0.99, by = 0.01)
x_grid <- -log(1 - p_grid)

get_X_matrix <- function(obj) {
  if (is.null(obj)) return(NULL)

  if (is.matrix(obj) || is.data.frame(obj)) return(as.matrix(obj))

  if (is.list(obj) && !is.null(obj$X)) {
    if (is.matrix(obj$X) || is.data.frame(obj$X)) return(as.matrix(obj$X))
  }

  NULL
}

cond_curve <- function(list_obj, s1, s2, s3, p_grid) {
  num <- den <- numeric(length(p_grid))

  for (i in seq_along(list_obj)) {
    X <- get_X_matrix(list_obj[[i]])
    if (is.null(X)) next
    if (!all(c(s1, s2, s3) %in% colnames(X))) next

    s3v <- X[, s3]
    if (all(is.na(s3v))) next

    z <- as.numeric(quantile(s3v, probs = p_grid, na.rm = TRUE, names = FALSE))
    if (any(!is.finite(z))) next

    for (k in seq_along(p_grid)) {
      ind <- s3v > z[k]
      ind[is.na(ind)] <- FALSE
      nk <- sum(ind)
      if (nk == 0) next

      den[k] <- den[k] + nk

      s1v <- X[ind, s1]
      s2v <- X[ind, s2]
      num[k] <- num[k] + sum((s1v > z[k]) & (s2v > z[k]), na.rm = TRUE)
    }
  }

  ifelse(den > 0, num / den, NA_real_)
}

sim_obj_list <- unlist(sims_by_ep, recursive = FALSE)


prob_obs <- cond_curve(list_episodes_filtered, s1_pick, s2_pick, s3_pick, p_grid)
prob_sim <- cond_curve(sim_obj_list, s1_pick, s2_pick, s3_pick, p_grid)

df_curve <- data.frame(
  x      = rep(x_grid, 2),
  prob   = c(prob_obs, prob_sim),
  source = rep(c("Observed", "Simulated"), each = length(x_grid))
)

ggplot(df_curve, aes(x = x, y = prob, linetype = source)) +
  geom_line(linewidth = 1, na.rm = TRUE) +
  geom_point(size = 1.5, na.rm = TRUE) +
  labs(
    title = paste0(
      "Conditional curve: P(", s1_pick, ">q, ", s2_pick, ">q | ", s3_pick, ">q)"
    ),
    x = "-log(1-p)",
    y = paste0("P(", s1_pick, ", ", s2_pick, " | ", s3_pick, ")"),
    linetype = "Source"
  ) +
  theme_minimal() +
  theme(panel.grid = element_blank())

head(df_curve)




u_by_s0_from_list <- tapply(u_list, s0_list, median, na.rm = TRUE)

u_by_s0_from_list

u_global <- u_by_s0_from_list


cond_probs_episode_q <- function(X, u_vec, sites = colnames(X)) {
  X <- as.matrix(X)[, sites, drop = FALSE]
  p <- length(sites)

  E <- sweep(X, 2, u_vec[sites], `>`)
  E[is.na(E)] <- FALSE

  denom <- colSums(E)

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

weighted_avg_cond_probs <- function(res_list) {
  p <- dim(res_list[[1]]$P)[1]
  sites <- dimnames(res_list[[1]]$P)$s1

  num <- array(0, dim = c(p, p, p), dimnames = list(s1=sites, s2=sites, s3=sites))
  den <- rep(0, p); names(den) <- sites

  for (r in res_list) {
    P <- r$P
    d <- r$denom
    for (k in seq_len(p)) {
      if (!is.finite(d[k]) || d[k] <= 0) next
      num[, , k] <- num[, , k] + P[, , k] * d[k]
      den[k] <- den[k] + d[k]
    }
  }

  out <- array(NA_real_, dim = c(p, p, p), dimnames = list(s1=sites, s2=sites, s3=sites))
  for (k in seq_len(p)) {
    if (den[k] > 0) out[, , k] <- num[, , k] / den[k]
  }
  out
}



obs_res <- lapply(seq_along(list_episodes_filtered), function(j) {
  X <- list_episodes_filtered[[j]]
  cond_probs_episode_q(X, u_global, sites = sites)
})
P_obs_avg <- weighted_avg_cond_probs(obs_res)

reps <- lapply(seq_along(sims_by_ep), function(j) {
  u <- u_global
  lapply(seq_along(sims_by_ep[[j]]), function(m) {
    Xsim <- get_X_matrix(sims_by_ep[[j]][[m]])
    cond_probs_episode_q(Xsim, u, sites = sites)
  })
})

P_list <- lapply(reps, `[[`, "P")
den_list <- lapply(reps, `[[`, "denom")

p <- dim(P_list[[1]])[1]
Pm_num <- array(0, dim = dim(P_list[[1]]), dimnames = dimnames(P_list[[1]]))
Pm_den <- rep(0, p); names(Pm_den) <- dimnames(P_list[[1]])$s3

for (m in seq_along(P_list)) {
  Pm_ <- P_list[[m]]
  d_  <- den_list[[m]]

  for (k in seq_len(p)) {
    if (!is.finite(d_[k]) || d_[k] <= 0) next
    tmp <- Pm_[,,k]
    tmp[is.na(tmp)] <- 0
    Pm_num[,,k] <- Pm_num[,,k] + tmp * d_[k]
    Pm_den[k] <- Pm_den[k] + d_[k]
  }
}

Pm <- array(NA_real_, dim = dim(Pm_num), dimnames = dimnames(Pm_num))
for (k in seq_len(p)) {
  if (Pm_den[k] > 0) Pm[,,k] <- Pm_num[,,k] / Pm_den[k]
}

denom_m <- Reduce(`+`, den_list) / length(den_list)

cond_probs_episode_q_na <- function(X, u_vec, sites = colnames(X)) {
  X <- as.matrix(X)[, sites, drop = FALSE]
  p <- length(sites)

  if (length(u_vec) == 1) u_vec <- rep(u_vec, p)
  if (is.null(names(u_vec))) names(u_vec) <- sites
  u <- u_vec[sites]

  P <- array(NA_real_, dim = c(p, p, p),
             dimnames = list(s1 = sites, s2 = sites, s3 = sites))
  denom <- rep(0, p); names(denom) <- sites

  for (k in seq_len(p)) {
    xk <- X[, k]
    ok_k <- !is.na(xk)

    ind_k <- ok_k & (xk > u[k])
    nk <- sum(ind_k)
    denom[k] <- nk
    if (nk == 0) next

    Xk <- X[ind_k, , drop = FALSE]

    E <- sweep(Xk, 2, u, `>`)

    for (i in seq_len(p)) {
      ei <- E[, i]
      ok_i <- !is.na(ei)
      for (j in seq_len(p)) {
        ej <- E[, j]
        ok_j <- !is.na(ej)
        ok_ij <- ok_i & ok_j
        if (!any(ok_ij)) next
        P[i, j, k] <- sum(ei[ok_ij] & ej[ok_ij]) / nk
      }
    }
  }

  list(P = P, denom = denom)
}


sim_res_by_ep <- lapply(seq_along(sims_by_ep), function(j) {
  Xobs_j <- list_episodes_filtered[[j]]

  reps <- lapply(seq_along(sims_by_ep[[j]]), function(m) {
    Xsim <- get_X_matrix(sims_by_ep[[j]][[m]])
    Xsim <- Xsim[, sites, drop=FALSE]
    Xsim_masked <- apply_obs_na_mask(Xsim, Xobs_j, sites)

    cond_probs_episode_q_na(Xsim_masked, u_global, sites = sites)
  })

  P_list <- lapply(reps, `[[`, "P")
  den_list <- lapply(reps, `[[`, "denom")

  p <- dim(P_list[[1]])[1]
  num <- array(0, dim=dim(P_list[[1]]), dimnames=dimnames(P_list[[1]]))
  den <- rep(0, p); names(den) <- dimnames(P_list[[1]])$s3

  for (mm in seq_along(P_list)) {
    Pm <- P_list[[mm]]; d <- den_list[[mm]]
    for (k in seq_len(p)) if (d[k] > 0) {
      tmp <- Pm[,,k]; tmp[is.na(tmp)] <- 0
      num[,,k] <- num[,,k] + tmp * d[k]
      den[k] <- den[k] + d[k]
    }
  }
  Pbar <- array(NA_real_, dim=dim(num), dimnames=dimnames(num))
  for (k in seq_len(p)) if (den[k] > 0) Pbar[,,k] <- num[,,k] / den[k]

  list(P = Pbar, denom = Reduce(`+`, den_list)/length(den_list))
})

P_sim_avg <- weighted_avg_cond_probs(sim_res_by_ep)

s3_pick <- "poly"
P_obs_s3 <- P_obs_avg[, , s3_pick]
P_sim_s3 <- P_sim_avg[, , s3_pick]

cond_probs_to_df2d <- function(P2d, source) {
  df <- melt(P2d, varnames = c("s1", "s2"), value.name = "prob")
  df$source <- source
  df
}

df_plot <- rbind(
  cond_probs_to_df2d(P_obs_s3, "Observed"),
  cond_probs_to_df2d(P_sim_s3, "Simulated")
)
head(df_plot)

ggplot(df_plot, aes(x = s1, y = s2, fill = value)) +
  geom_tile() +
  facet_wrap(~ source, nrow = 1) +
  coord_equal() +
  scale_fill_gradient(low = "white", high = "#357470",
                      limits = c(0, 1), na.value = "grey90") +
  labs(
    fill = "Probability"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size =15),
    axis.text.y = element_text(size = 15),
    panel.grid = element_blank(),
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 16),
    strip.text = element_text(size = 16, face = "bold"),
    axis.title = element_text(size = 16)
  )

foldername <- paste0(im_folder, "/swg/omsev/2025/cond_probs/")
# if (!dir.exists(foldername)) {
#   dir.create(foldername, recursive = TRUE)
# }
ggsave(paste0(foldername, "cond_probs_", s3_pick, ".png"), width = 15, height = 7)






