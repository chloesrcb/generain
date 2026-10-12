library(dplyr)
library(ggplot2)
library(mev)
library(POT)
library(patchwork)

best_cens <- 0.23

q_levels <- c(
  0.95,
  0.96,
  0.97,
  0.98,
  0.99
)

################################################################################
# Pool all rain gauges
################################################################################

y <- unlist(rain, use.names = FALSE)

y <- y[
  !is.na(y) &
  is.finite(y) &
  y > 0
]

cat("Number of positive observations:", length(y), "\n")
cat("Range:", range(y), "\n")


################################################################################
# EGPD - pooled fit
################################################################################

inits <- init_values(y, 0)

fit_egpd <- fit.extgp(
  y,
  model = 1,
  method = "mle",
  init = c(1, inits[1], inits[2]),
  censoring = c(best_cens, Inf),
  plots = FALSE,
  confint = FALSE,
  ncpus = 1
)

par_egpd <- fit_egpd$fit$mle

cat(
  "\nEGPD parameters:\n",
  "kappa =", par_egpd[1], "\n",
  "sigma =", par_egpd[2], "\n",
  "xi    =", par_egpd[3], "\n"
)

p <- ppoints(length(y))

qq_egpd_df <- data.frame(
  empirical = sort(y),
  fitted = qextgp(
    p,
    type = 1,
    kappa = par_egpd[1],
    sigma = par_egpd[2],
    xi = par_egpd[3]
  )
)


################################################################################
# Gamma truncated below threshold u
################################################################################

fit_gamma_trunc <- function(y, u) {

  yg <- y[y <= u]

  if(length(yg) < 20)
    return(NULL)

  nll <- function(theta) {

    shape <- exp(theta[1])
    rate  <- exp(theta[2])

    Fu <- pgamma(
      u,
      shape = shape,
      rate = rate
    )

    -sum(
      dgamma(
        yg,
        shape = shape,
        rate = rate,
        log = TRUE
      )
    ) +
      length(yg) * log(Fu)
  }

  init <- c(
    log(2),
    log(2 / mean(yg))
  )

  opt <- optim(
    init,
    nll
  )

  list(
    shape = exp(opt$par[1]),
    rate = exp(opt$par[2])
  )
}


################################################################################
# GPD above threshold u
################################################################################

fit_gpd <- function(y, u) {

  exc <- y[y > u] - u

  if(length(exc) < 20)
    return(NULL)

  fit <- tryCatch(
    mev::gp.fit(
      exc,
      threshold = 0
    ),
    error = function(e) NULL
  )

  if(is.null(fit))
    return(NULL)

  list(
    sigma = fit$estimate[1],
    xi = fit$estimate[2]
  )
}


################################################################################
# Quantile function Gamma + GPD
################################################################################

q_mix <- function(
    p,
    u,
    shape,
    rate,
    sigma,
    xi,
    pi
) {

  Fu <- pgamma(
    u,
    shape = shape,
    rate = rate
  )

  out <- numeric(length(p))

  lower <- p <= pi

  out[lower] <- qgamma(
    p[lower] / pi * Fu,
    shape = shape,
    rate = rate
  )

  upper <- !lower

  out[upper] <- u +
    qgpd(
      (p[upper] - pi) / (1 - pi),
      loc = 0,
      scale = sigma,
      shape = xi
    )

  out
}


################################################################################
# QQ dataframe Gamma + GPD
################################################################################

qq_mix_df <- function(y, q) {

  u <- as.numeric(
    quantile(
      y,
      probs = q,
      type = 8
    )
  )

  gfit <- fit_gamma_trunc(
    y,
    u
  )

  pfit <- fit_gpd(
    y,
    u
  )

  if(is.null(gfit) || is.null(pfit))
    return(NULL)

  pi <- mean(y <= u)

  p <- ppoints(length(y))

  data.frame(
    empirical = sort(y),

    fitted = q_mix(
      p,
      u,
      gfit$shape,
      gfit$rate,
      pfit$sigma,
      pfit$xi,
      pi
    ),

    model = paste0(
      "q=",
      q,
      " (u=",
      round(u, 3),
      ")"
    )
  )
}


################################################################################
# Fit all Gamma + GPD thresholds
################################################################################

mix_list <- lapply(
  q_levels,
  function(q) {
    qq_mix_df(
      y,
      q
    )
  }
)

mix_list <- mix_list[
  !sapply(
    mix_list,
    is.null
  )
]

df_mix <- bind_rows(
  mix_list
)


################################################################################
# Common plotting range
################################################################################

common_ylim <- range(
  c(
    qq_egpd_df$fitted,
    df_mix$fitted
  ),
  finite = TRUE
)


################################################################################
# Separate plots
################################################################################

p_egpd <- ggplot(
  qq_egpd_df,
  aes(
    empirical,
    fitted
  )
) +
  geom_point(
    colour = "#d95f02",
    size = 1.2,
    alpha = 0.7
  ) +
  geom_abline(
    slope = 1,
    intercept = 0,
    linetype = 2
  ) +
  theme_minimal() +
  labs(
    title = paste0(
      "EGPD - all rain gauges (cens = ",
      best_cens,
      ")"
    ),
    x = "Empirical quantiles",
    y = "Theoretical quantiles"
  ) +
  coord_cartesian(
    ylim = common_ylim
  )


p_mix <- ggplot(
  df_mix,
  aes(
    empirical,
    fitted,
    colour = model
  )
) +
  geom_point(
    size = 1.2,
    alpha = 0.7
  ) +
  geom_abline(
    slope = 1,
    intercept = 0,
    linetype = 2
  ) +
  theme_minimal() +
  labs(
    title = "Gamma (truncated) + GPD - all rain gauges",
    x = "Empirical quantiles",
    y = "Theoretical quantiles",
    colour = "Threshold"
  ) +
  coord_cartesian(
    ylim = common_ylim
  )


p_egpd / p_mix


################################################################################
# Combined comparison plot
################################################################################

df_mix_plot <- df_mix %>%
  mutate(
    type = "Gamma+GPD"
  )

df_egpd_plot <- qq_egpd_df %>%
  mutate(
    type = "EGPD"
  )


p_compare <- ggplot() +

  geom_point(
    data = df_mix_plot,
    aes(
      empirical,
      fitted,
      colour = model
    ),
    size = 1.4,
    alpha = 0.65
  ) +

  geom_point(
    data = df_egpd_plot,
    aes(
      empirical,
      fitted
    ),
    colour = "#f31414",
    shape = 17,
    size = 2,
    alpha = 0.8
  ) +

  geom_abline(
    slope = 1,
    intercept = 0,
    linetype = 2
  ) +

  theme_minimal() +

  labs(
   
    x = "Empirical quantiles",
    y = "Theoretical quantiles",
    colour = "Gamma+GPD threshold"
  ) +

  coord_cartesian(
    ylim = common_ylim
  ) +
  btf_theme

p_compare


################################################################################
# Save
################################################################################

save_path <- paste0(
  im_folder,
  "/EGPD/OMSEV/2019_2025/comparison_gamma_gpd/"
)

if(!dir.exists(save_path)) {
  dir.create(
    save_path,
    recursive = TRUE
  )
}


ggsave(
  p_egpd / p_mix,
  filename = paste0(
    save_path,
    "all_sites.png"
  ),
  width = 8,
  height = 12
)


ggsave(
  p_compare,
  filename = paste0(
    save_path,
    "all_sites_compare.png"
  ),
  width = 10,
  height = 8
)
