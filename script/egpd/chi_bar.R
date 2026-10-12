################################################################################
# CHI AND CHI-BAR --------------------------------------------------------------
################################################################################

library(dplyr)
library(ggplot2)

# ------------------------------------------------------------------------------
# Function to compute chi(u) and chi_bar(u)
#
# IMPORTANT:
# remove_zeros = TRUE means that dependence is studied CONDITIONAL on both
# observations being positive.
#
# lag = 0   : same time
# lag = 1   : X_t compared with Y_{t+1}
# etc.
# ------------------------------------------------------------------------------

compute_chi_chibar <- function(x,
                               y,
                               lag = 0,
                               probs = seq(0.90, 0.995, by = 0.002),
                               remove_zeros = TRUE) {

  # --------------------------------------------------------------------------
  # Align series according to lag
  # --------------------------------------------------------------------------

  if (lag > 0) {

    x <- x[1:(length(x) - lag)]
    y <- y[(lag + 1):length(y)]

  } else if (lag < 0) {

    lag_abs <- abs(lag)

    x <- x[(lag_abs + 1):length(x)]
    y <- y[1:(length(y) - lag_abs)]

  }

  # --------------------------------------------------------------------------
  # Remove missing values
  # --------------------------------------------------------------------------

  keep <- is.finite(x) & is.finite(y)

  x <- x[keep]
  y <- y[keep]

  # --------------------------------------------------------------------------
  # Remove zeros
  #
  # IMPORTANT: remove rows where either value is zero.
  # Do NOT remove zeros independently from x and y.
  # --------------------------------------------------------------------------

  if (remove_zeros) {

    keep <- x > 0 & y > 0

    x <- x[keep]
    y <- y[keep]

  }

  n <- length(x)

  if (n < 10) {
    stop("Not enough observations after filtering.")
  }

  # --------------------------------------------------------------------------
  # Empirical probability integral transform
  #
  # rank/(n+1) avoids exact values 0 and 1
  # --------------------------------------------------------------------------

  Ux <- rank(x, ties.method = "average") / (n + 1)
  Uy <- rank(y, ties.method = "average") / (n + 1)

  # --------------------------------------------------------------------------
  # Compute chi and chi-bar
  # --------------------------------------------------------------------------

  out <- lapply(probs, function(u) {

    exc_x <- Ux > u
    exc_y <- Uy > u

    joint <- exc_x & exc_y

    p_joint <- mean(joint)

    # empirical chi
    #
    # Since Ux and Uy have uniform margins:
    # P(Ux > u) ~= 1-u
    #
    chi <- p_joint / (1 - u)

    # chi-bar
    #
    # undefined if there is no joint exceedance
    #
    if (p_joint > 0 && p_joint < 1) {

      chibar <- 2 * log(1 - u) / log(p_joint) - 1

    } else {

      chibar <- NA_real_

    }

    data.frame(
      u = u,
      chi = chi,
      chibar = chibar,
      n = n,
      n_exc_x = sum(exc_x),
      n_exc_y = sum(exc_y),
      n_joint = sum(joint),
      p_joint = p_joint
    )

  })

  bind_rows(out)
}







sites <- c("cnrs", "um", "poly", "mse", "cefe", "hydro", "um35", "brives", "archiw", "iem", "chu1")

probs <- seq(0.7, 0.999, by = 0.002)

# all unique pairs
pairs <- combn(sites, 2, simplify = FALSE)

res_all_chibar <- lapply(pairs, function(pair) {

  s1 <- pair[1]
  s2 <- pair[2]

  res <- compute_chi_chibar(
    x = rain[[s1]],
    y = rain[[s2]],
    lag = 0,
    probs = probs,
    remove_zeros = TRUE
  )

  res$site1 <- s1
  res$site2 <- s2
  res$pair <- paste(s1, s2, sep = " - ")

  res
})

res_all_chibar <- bind_rows(res_all_chibar)

head(res_all_chibar)



################################################################################
# SUMMARY ACROSS PAIRS ---------------------------------------------------------
################################################################################

summary_chibar <- res_all_chibar %>%
  group_by(u) %>%
  summarise(
    median_chibar = median(chibar, na.rm = TRUE),
    q25 = quantile(chibar, 0.25, na.rm = TRUE),
    q75 = quantile(chibar, 0.75, na.rm = TRUE),
    min_chibar = min(chibar, na.rm = TRUE),
    max_chibar = max(chibar, na.rm = TRUE),
    .groups = "drop"
  )



################################################################################
# PLOT CHI - ALL PAIRS ---------------------------------------------------------
################################################################################

ggplot(
  res_all_chibar,
  aes(
    x = u,
    y = chi,
    group = pair,
    colour = pair
  )
) +
  geom_hline(
    yintercept = 0,
    linetype = "dotted",
    linewidth = 0.5
  ) +
  geom_line(
    linewidth = 0.8,
    alpha = 0.8
  ) +
  coord_cartesian(
    xlim = c(0.95, 0.999),
    ylim = c(0, 1)
  ) +
  labs(
    x = "Quantile q",
    y = expression(hat(chi)(q)),
    colour = "Pair"
  ) +
  btf_theme +
  theme(
    legend.position = "none"
  )

ggsave(
  filename = paste0(foldername, "/chi_all_pairs.pdf"),
  width = 8,
  height = 7
)





ggplot(
  res_all_chibar,
  aes(
    x = u,
    y = chibar,
    group = pair,
    colour = pair
  )
) +
  geom_hline(
    yintercept = 0,
    linetype = "dotted",
    linewidth = 0.5
  ) +
  geom_line(
    linewidth = 0.8,
    alpha = 0.8
  ) +
  coord_cartesian(
    xlim = c(0.95, 0.999),
    ylim = c(0, 1)
  ) +
  labs(
    x = "Quantile q",
    y = expression(bar(chi)(q)),
    colour = "Pair"
  ) +
  btf_theme +
  theme(
    legend.position = "none"
  ) +
  # add y lines for y=1
  geom_hline(yintercept = 1, linetype = "dashed", linewidth = 0.5)


ggsave(
  filename = paste0(foldername, "/chibar_all_pairs.pdf"),
  width = 8,
  height = 5
)
