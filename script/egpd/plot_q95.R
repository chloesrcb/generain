library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)

sites <- colnames(rain)

q <- 0.95

u95 <- apply(rain[, sites], 2, function(x) {
  x_pos <- x[!is.na(x) & x > 0]
  quantile(x_pos, probs = 0.95, na.rm = TRUE)
})

E <- sweep(as.matrix(rain[, sites]), 2, u95, `>`)
E[is.na(E)] <- FALSE

location_gauges <- location_gauges %>%
  filter(Station %in% sites) %>%
  arrange(match(Station, sites))

sites <- location_gauges$Station

E <- E[, sites, drop = FALSE]

df_spatial <- expand.grid(
  site1 = sites_order,
  site2 = sites_order,
  stringsAsFactors = FALSE
) %>%
  mutate(
    i1 = match(site1, sites_order),
    i2 = match(site2, sites_order)
  ) %>%
  filter(i1 <= i2) %>%  
  rowwise() %>%
  mutate(
    n_exceed = sum(E[, site1] & E[, site2], na.rm = TRUE)
  ) %>%
  ungroup()

df_spatial$site1 <- factor(df_spatial$site1, levels = sites_order)
df_spatial$site2 <- factor(df_spatial$site2, levels = rev(sites_order))
min_pos <- min(df_spatial$n_exceed[df_spatial$n_exceed > 0], na.rm = TRUE)

p_spatial <- ggplot(df_spatial, aes(x = site1, y = site2)) +
  geom_point(
    data = subset(df_spatial, n_exceed > 0),
    aes(size = n_exceed),
    alpha = 0.8,
    color = "#357470"
  ) +
  geom_text(
    data = subset(df_spatial, n_exceed == 0),
    aes(label = "0"),
    size = 5,
    color = "grey30"
  ) +
  geom_text(
    data = subset(df_spatial, n_exceed == min_pos),
    aes(label = n_exceed),
    size = 4,
    vjust = -0.4
  ) +
  scale_size_continuous(name = "Number of exceedances") +
  coord_equal() +
  labs(
    x = expression(paste("Site ", s[1])),
    y = expression(paste("Site ", s[2]))
  ) +
  btf_theme +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "bottom",
    legend.key.size = unit(0.5, "cm"),
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 12)
  )

p_spatial


# save the plot
foldername <- paste0(im_folder, "/swg/omsev/2025/")
if (!dir.exists(foldername)) dir.create(foldername, recursive = TRUE)
ggsave(
  paste0(foldername, "joint_exceedances_spatial_q", q * 100, ".png"),
  plot = p_spatial,
  width = 8,
  height = 8
)


max_lag <- 10

df_temporal <- lapply(sites, function(s) {
  e <- E[, s]

  data.frame(
    site = s,
    lag = 0:max_lag,
    n_exceed = sapply(0:max_lag, function(h) {
      if (h == 0) {
        sum(e, na.rm = TRUE)
      } else {
        sum(e[1:(length(e) - h)] & e[(1 + h):length(e)], na.rm = TRUE)
      }
    })
  )
}) %>%
  bind_rows()

p_temporal <- ggplot(df_temporal, aes(x = lag, y = n_exceed, color = site)) +
  geom_line(linewidth = 0.7) +
  geom_point(size = 1.3) +
  scale_x_continuous(breaks = 0:max_lag) +
  labs(
    x = "Temporal lag (5 min intervals)",
    y = "Number of exceedances",
    color = "Site"
  ) +
  # add line for number of exceedances = 20 
  geom_hline(yintercept = 20, linetype = "dashed", color = "#2c2b2b", linewidth = 0.7) +
  annotate(
    "text",
    x = 0,
    y = 20 + 4,
    label = "20",
    color = "#0e0b0b",
    size = 4,
    hjust = 1
  ) +
  xlim(0, max_lag) +
  theme(
    axis.text = element_text(size = 9),
    axis.title = element_text(size = 11),
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 12),
  ) +
  btf_theme


ggsave(
  paste0(foldername, "joint_exceedances_temporal_q", q * 100, ".png"),
  plot = p_temporal,
  width = 10,
  height = 8
)
