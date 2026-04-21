Nsim = 100
# get sims_group for one adv group
sims_sign <- readRDS(file = paste0(
  data_folder, "/omsev/simulated_episodes/sims_episodes_advgroup_",
  "significant", "_n", Nsim, ".rds"
) )

sims_weak <- readRDS(file = paste0(
  data_folder, "/omsev/simulated_episodes/sims_episodes_advgroup_",
  "weak", "_n", Nsim, ".rds"
) )

sims_still <- readRDS(file = paste0(
  data_folder, "/omsev/simulated_episodes/sims_episodes_advgroup_",
  "still", "_n", Nsim, ".rds"
) )

# combine all with sims_sign, sims_weak, sims_still
sims_group <- c(sims_sign, sims_weak, sims_still)
# add adv_group info
group_adv <- rep(c("significant", "weak", "still"), each = Nsim)
for( i in seq_along(sims_group)) {
  sims_group[[i]]$adv_group <- group_adv[i]
}


X <- sims_group[[1]]$X  # matrix of size n x p
site <- "crbm"
u <- 1  # threshold
# Computes P(X_s>u & X_s'>u | X_s*>u) = P(X_s>u & X_s'>u & X_s*>u) / P(X_s*>u)
cond_probs_all_sites <- function(X, u, site) {
  sites <- colnames(X)
  p <- length(sites)
  n <- nrow(X)

  E_site <- sum(X[, site] > u) # P(X_s*>u)
  cond_probs <- matrix(0, nrow = p, ncol = p,
                       dimnames = list(sites, sites))
  if(E_site == 0) {
    return(cond_probs)
  } else {
    for(s in seq_len(p)) {
      for(s_prime in seq_len(p)) {
        E_ssp_u <- sum(X[, site] > u & X[, s] > u & X[, s_prime] > u)
        cond_probs[s, s_prime] <- E_ssp_u / E_site
      }
    }
  }
  return(cond_probs)
}

X
X <- sims_group[[100]]$X
u <- sims_group[[100]]$u_emp
site <- "chu7"
cond_probs_sim <- cond_probs_all_sites(X, u, site)
episode_obs <- list_episodes[[2]]
cond_probs_obs <- cond_probs_all_sites(episode_obs, u, site)


sites <- colnames(X)
sim_X <- sims_group[[2]]$X
sim_s0 <- sims_group[[2]]$s0
sim_u <- sims_group[[2]]$u_emp
# for all sites s, s' in sites and all s* except s0 compute P(X_s>u, X_s'>u | X_s*>u)
# create a list of matrices with names the sites excluding s0
cond_probs_sim_excl_s0 <- list()

for(s_star in sites[sites != sim_s0]) {
    cond_probs_sim_excl_s0[[s_star]] <- cond_probs_all_sites(sim_X, sim_u, s_star)
}

library(ggplot2)
library(dplyr)
library(tidyr)

df_plot <- data.frame()
for(s in sites) {
  for(s_star in sites) {
    if(s != sim_s0 & s_star != sim_s0) {
      prob_sstar <- mean(sim_X[, s_star] > sim_u)
      cond_prob <- cond_probs_sim_excl_s0[[s_star]][s, s_star]
      df_plot <- rbind(
        df_plot,
        data.frame(
          site = s,
          site_prime = s_star,
          prob_sstar = prob_sstar,
          neg_log_prob_sstar = -log(prob_sstar),
          cond_prob = cond_prob
        )
      )
    }
  }
}

df_plot <- df_plot %>%
  filter(!is.infinite(neg_log_prob_sstar))

ggplot(df_plot, aes(x = neg_log_prob_sstar, y = cond_prob)) +
  geom_point(alpha = 0.5) +
  geom_smooth(method = "loess") +
  labs(
    title = paste0("Conditional exceedance probabilities given site ", sim_s0,
                   " exceeds threshold ", sim_u),
    x = "-log(P(X_s*>u))",
    y = "P(X_s>u & X_s'>u | X_s*>u)"
  ) +
  theme_minimal()



cond_prob_bivar <- function(X, u, s_star) {
  sites <- colnames(X)
  I <- (X > u) * 1L
  k <- which(sites == s_star)
  idx <- which(I[, k] == 1L)
  p <- ncol(I)

  out <- matrix(0, p, p, dimnames = list(sites, sites))
  if (length(idx) == 0) return(out)

  J <- I[idx, , drop = FALSE]
  out <- crossprod(J) / nrow(J)
  out
}

cond_prob_bivar_all_sstar <- function(X, u, exclude = NULL) {
  sites <- colnames(X)
  if (!is.null(exclude)) sites <- setdiff(sites, exclude)

  res <- setNames(vector("list", length(sites)), sites)
  for (s_star in sites) {
    res[[s_star]] <- cond_prob_bivar(X, u, s_star)
  }
  res
}

library(dplyr)
library(tidyr)

to_long_bivar <- function(list_mats, X, u) {
  sites <- colnames(X)

  bind_rows(lapply(names(list_mats), function(s_star) {
    prob_sstar <- mean(X[, s_star] > u)
    as.data.frame(list_mats[[s_star]]) |>
      mutate(s = row.names(.)) |>
      pivot_longer(-s, names_to = "t", values_to = "cond_prob") |>
      mutate(
        s_star = s_star,
        prob_sstar = prob_sstar,
        neg_log_prob_sstar = -log(prob_sstar)
      )
  }))
}


sites <- colnames(sim_X)

mats <- cond_prob_bivar_all_sstar(sim_X, sim_u, exclude = sim_s0)

df_bivar <- to_long_bivar(mats, sim_X, sim_u) %>%
  filter(is.finite(neg_log_prob_sstar)) %>%
  filter(s != s_star, t != s_star) %>%
  identity()

library(ggplot2)

ggplot(df_bivar, aes(x = neg_log_prob_sstar, y = cond_prob)) +
  geom_point(alpha = 0.3) +
  geom_smooth(method = "loess", se = FALSE) +
  facet_wrap(~ s_star) +
  labs(
    title = "P(X_s>u & X_t>u | X_{s*}>u) vs -log P(X_{s*}>u)",
    x = "-log(P(X_{s*}>u))",
    y = "P(X_s>u & X_t>u | X_{s*}>u)"
  ) +
  theme_minimal()


all_obs <- bind_rows(lapply(seq_along(list_episodes), function(e) {
  Xobs <- list_episodes[[e]]
  u_obs <- u_list[[e]]
  mats_e <- cond_prob_bivar_all_sstar(Xobs, u_obs)
  to_long_bivar(mats_e, Xobs, u_obs) %>% mutate(episode = e)
}))

# for one s_star
s_star <- "cefe"

df_obs_sstar <- all_obs %>%
  filter(s_star == s_star) %>%
  filter(s != s_star, t != s_star) %>%
  filter(is.finite(neg_log_prob_sstar)) %>%
  identity()

# mean over episodes
df_obs_sstar_mean <- df_obs_sstar %>%
  group_by(neg_log_prob_sstar) %>%
  summarise(
    cond_prob = mean(cond_prob, na.rm = TRUE),
    .groups = "drop"
  )
x11()
plot(df_obs_sstar_mean$neg_log_prob_sstar, df_obs_sstar_mean$cond_prob)
