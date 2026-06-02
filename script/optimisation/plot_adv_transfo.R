# do transformations with eta1 and eta2 on magnitude
eta1 = 3.8
eta2 = 2.2

adv_transformed <- adv_df
adv_transformed$adv_speed <- sqrt(adv_df$vx_final^2 + adv_df$vy_final^2)
adv_transformed$adv_direction <- atan2(adv_df$vy_final, adv_df$vx_final) * (180 / pi)
adv_transformed$adv_speed <- eta1 * adv_transformed$adv_speed^eta2

# plot adv_transformed speed vs original speed
ggplot(adv_transformed, aes(x = sqrt(vx_final^2 + vy_final^2), y = adv_speed)) +
  geom_point() +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
  labs(x = "Original advection speed (km/h)", y = "Transformed advection speed", title = "Advection speed transformation") +
  theme_bw()
# plot transformed advection vectors
library(dplyr)
library(ggplot2)

adv_plot <- bind_rows(
  adv_df %>%
    transmute(
      dx = vx_final,
      dy = vy_final,
      speed = sqrt(dx^2 + dy^2),
      type = "Empirical"
    ),
  adv_transformed %>%
    transmute(
      dx = adv_speed * cos(adv_direction * pi / 180),
      dy = adv_speed * sin(adv_direction * pi / 180),
      speed = adv_speed,
      type = "Transformed"
    )
)

adv_rose <- adv_plot %>%
  mutate(
    direction_math = atan2(dy, dx) * 180 / pi,
    direction_math = ifelse(direction_math < 0, direction_math + 360, direction_math),

    # convention rose : 0 = N, 90 = E, 180 = S, 270 = W
    direction_rose = (90 - direction_math) %% 360,

    dir_bin = cut(
      direction_rose,
      breaks = seq(0, 360, by = 22.5),
      include.lowest = TRUE
    ),
    speed_bin = cut(
      speed,
      breaks = c(0, 2, 5, 10, 20, 50, 100, Inf),
      labels = c("0–2", "2–5", "5–10", "10–20", "20–50", "50–100", ">100")
    )
  )

ggplot(adv_rose, aes(x = direction_rose, fill = speed_bin)) +
  geom_histogram(
    binwidth = 22.5,
    boundary = 0,
    color = "white",
    linewidth = 0.25,
    alpha = 0.7
  ) +
  coord_polar(start = 0, direction = 1) +
  facet_wrap(~type) +
  scale_x_continuous(
    limits = c(0, 360),
    breaks = c(0, 90, 180, 270),
    labels = c("N", "E", "S", "W")
  )

ggplot(adv_plot, aes(x = 0, y = 0, xend = dx, yend = dy, color = type)) +
  geom_segment(
    arrow = arrow(length = unit(0.16, "cm")),
    linewidth = 0.5
  ) +
  geom_hline(yintercept = 0, color = "grey70") +
  geom_vline(xintercept = 0, color = "grey70") +
  coord_equal(xlim = c(-10, 10), ylim = c(-10, 10)) +
  scale_alpha_continuous(range = c(0.15, 0.8), guide = "none") +
  labs(
    x = expression(V[x]~"(km/h)"),
    y = expression(V[y]~"(km/h)"),
    color = NULL,
    title = "Advection vectors",
    subtitle = "Empirical vs transformed advection"
  ) +
  theme_bw() +
  theme(
    legend.position = "top",
    plot.title = element_text(face = "bold")
  )



adv_rose <- adv_rose[!is.na(adv_rose$dir_bin) & !is.na(adv_rose$speed_bin), ]



ggplot(adv_rose, aes(x = direction_rose, fill = speed_bin)) +
  geom_histogram(
    binwidth = 22.5,
    boundary = 0,
    color = "white",
    linewidth = 0.25,
    alpha = 0.7
  ) +
  coord_polar(start = 0, direction = 1) +
  facet_wrap(~type) +
  scale_x_continuous(
    limits = c(0, 360),
    breaks = c(0, 90, 180, 270),
    labels = c("N", "E", "S", "W")
  ) +
  labs(
    x = NULL,
    y = NULL,
    fill = "Speed (km/h)"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(face = "bold", size = 12),
    axis.text.y = element_blank(),
    axis.ticks = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(color = "grey85", linewidth = 0.3),
    legend.position = "right",
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 9)
  )

# save rose plot
foldername <- paste0(im_folder, "advection/withtransform/")
if (!dir.exists(foldername)) {
  dir.create(foldername, recursive = TRUE)
}
filename_rose <- paste0(foldername, "omsev_adv_rose_q", q * 100, "_delta", delta, "_dmin", min_spatial_dist,
"eta1", eta1, "_eta2", eta2, ".png")


ggsave(filename_rose, width = 8, height = 4, dpi = 300)
