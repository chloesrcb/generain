# remove objects from the environment
rm(list = ls())

# clear console
cat("\014")

# Load libraries and set theme
source("./script/load_libraries.R")
functions_folder <- "./R"
files <- list.files(functions_folder, full.names = TRUE)
# load all functions in files
invisible(lapply(files, function(f) source(f, echo = FALSE)))
library(latex2exp)
library(boot)
library(spdep)
library(dplyr)
library(tidyr)
library(kSamples)
library(twosamples)
library(ggrepel)

################################################################################
# DATA -------------------------------------------------------------------------
################################################################################
data_folder <- "../phd_extremes/data/"
filename_omsev <- paste0(data_folder,
                         "omsev/omsev_5min/rain_mtp_5min_2019_2025.csv")
rain_raw <- read.csv(filename_omsev)
rain <- rain_raw
head(rain)
tail(rain)

filename_loc <- paste0(data_folder,
                           "omsev/loc_rain_gauges.csv")
# get location of each rain gauge
location_gauges <- read.csv(filename_loc)

# compare rain.all5 and rain dataframes on same period
rain$dates <- as.POSIXct(rain$dates, format = "%Y-%m-%d %H:%M:%S", tz = "UTC")
rownames(rain) <- rain$dates
# remove rain$dates
rain$dates <- NULL

################################################################################
rain_chu1_pos <- rain$cnrs[rain$cnrs > 0]
min_precision <- min(rain_chu1_pos, na.rm = TRUE)
sort(unique(rain_chu1_pos))

results <- list()
sites_name <- colnames(rain)

# Get possible censoring values according to precision
min_censuring <- round(min_precision, 2)
# k_prec_censuring <- seq(min_censuring, min_censuring * 4, by = min_censuring)
k_prec_censuring <- seq(min_censuring, 2*min_censuring, by = 0.01)
# p <- min_censuring
# k_prec_censuring <- round(c(p, 1.5 * p, 2 * p, 2.5 * p, 3 * p), 4)
save_path <- paste0(im_folder,"/EGPD/OMSEV/2019_2025/continuous/")

if (!dir.exists(save_path)) {
  dir.create(save_path, recursive = TRUE)
}
# Fit EGPD for each site
for (site_name in sites_name) {
  cat("Processing:", site_name, "\n")
  y <- as.data.frame(na.omit(rain[[site_name]]))
  site_result <- tryCatch({
    process_site(y = y, site_name =  site_name,
                 save_path = save_path, R = 1000, forced_cens = k_prec_censuring)
  }, error = function(e) {
    warning(paste("Failed for site:", site_name))
    NULL
  })
  if (!is.null(site_result)) {
    results[[site_name]] <- site_result
  }
}


df_all_results <- do.call(rbind, results)
print(df_all_results)
# save results
write.csv(df_all_results, file = paste0(im_folder,"/EGPD/OMSEV/2019_2025/egpd_results.csv"), 
            row.names = FALSE)

df_all_results <- read.csv(paste0(im_folder,"/EGPD/OMSEV/2019_2025/egpd_results.csv"))


df_plot <- df_all_results %>% arrange(xi) 

foldername <- paste0(im_folder,"/EGPD/OMSEV/2019_2025/")
ggplot(df_plot, aes(x = reorder(Site, kappa), y = xi)) +
  geom_point(size = 3, color = "#0072B2") +
  geom_errorbar(aes(ymin = xi_low, ymax = xi_high), width = 0.2, color = "#0072B2") +
  coord_flip() +
  theme_minimal() +
  xlab("Site") + ylab("Parameter") +
  btf_theme 
ggsave(paste0(save_path, "xi_by_station.png"), width = 8, height = 6, dpi = 400)


ggplot(df_plot, aes(x = reorder(Site, kappa), y = sigma)) +
  geom_point(size = 3, color = "#009E73") +
  geom_errorbar(aes(ymin = sigma_low, ymax = sigma_high), width = 0.2, color = "#009E73") +
  coord_flip() +
  theme_minimal() +
  xlab("Site") + ylab("Parameter") +
  btf_theme 
ggsave(paste0(save_path, "sigma_by_station.png"), width = 8, height = 6, dpi = 400)


ggplot(df_plot, aes(x = reorder(Site, kappa), y = kappa)) +
  geom_point(size = 3, color = "#D55E00") +
  geom_errorbar(aes(ymin = kappa_low, ymax = kappa_high), width = 0.2, color = "#D55E00") +
  coord_flip() +
  theme_minimal() +
    xlab("Site") + ylab("Parameter") +
  btf_theme 
ggsave(paste0(save_path, "kappa_by_station.png"), width = 8, height = 6, dpi = 400)

# Do boxplots of left-censored values and RMSE

df_cens_rmse <- df_all_results %>%
  select(Site, BestCens, RMSE) %>%
  pivot_longer(cols = c(BestCens, RMSE), names_to = "Metric", values_to = "Value")

labels <- c(
  BestCens = "Best left-censoring",
  RMSE = "RMSE"
)
ggplot(df_cens_rmse, aes(x = Metric, y = Value)) +
  geom_boxplot(fill=btfgreen, alpha = 0.7) +
  theme_minimal() +
  xlab("") + ylab("Value") +
  scale_x_discrete(labels = labels) +
  theme(legend.position = "none") +
  btf_theme   

ggsave(paste0(save_path, "boxplot_cens_rmse.pdf"), width = 8, height = 6, dpi = 400)

# Do boxplots -------------------------------------------------------
df_melted <- df_all_results %>%
  select(Site, xi, sigma, kappa) %>%
  pivot_longer(cols = c(xi, sigma, kappa), names_to = "Parameter", values_to = "Value")
latex_names <- c(
    xi = TeX("$\\widehat{\\xi}$"),
    sigma = TeX("$\\widehat{\\sigma}$"),
    kappa = TeX("$\\widehat{\\kappa}$")
  )
ggplot(df_melted, aes(x = Parameter, y = Value)) +
  geom_boxplot(fill=btfgreen, alpha = 0.7) +
  theme_minimal() +
  xlab("Parameter") + ylab("Value") +
  scale_x_discrete(labels = latex_names) +
  theme(legend.position = "none") +
  ylim(0, 1.5) +
  btf_theme 
  
ggsave(paste0(save_path, "boxplot_parameters.pdf"), width = 8, height = 6, dpi = 400)
################################################################################
# Statistical tests -----------------------------------------------------------
################################################################################

library(spdep)
coords <- cbind(location_gauges$Longitude, location_gauges$Latitude)
# around 5 years
sites_name_5y <- c("cefe", "cnrs", "crbm", "iem", "mse",
                   "poly", "um", "archiw")

# around 3 years
sites_name_3y <- c("archie", "chu1", "chu2", "chu3",
                   "chu4", "chu5", "chu6", "chu7",
                   "um35")

# around 2 years
sites_name_2y <- c("cines", "hydro", "brives", "archie")

# add duration column
df_all_results$duration <- dplyr::case_when(
  df_all_results$Site %in% sites_name_5y ~ "5 ans",
  df_all_results$Site %in% sites_name_3y ~ "3 ans",
  df_all_results$Site %in% sites_name_2y ~ "2 ans",
  TRUE ~ "Inconnue"
)


# do kruskal wallis test -------------------------------------------------
kruskal.test(xi ~ duration, data = df_all_results)
kruskal.test(sigma ~ duration, data = df_all_results)
kruskal.test(kappa ~ duration, data = df_all_results)

# Adjust parameters by duration means
df_results <- df_all_results %>%
  mutate(
    xi_adj = xi - mean(xi, na.rm = TRUE),
    sigma_adj = sigma - mean(sigma, na.rm = TRUE),
    kappa_adj = kappa - mean(kappa, na.rm = TRUE)
  ) %>% ungroup()

location_gauges <- location_gauges %>%
  select(Station, Longitude, Latitude)
colnames(location_gauges) <- c("Site", "Longitude", "Latitude")


# Do Moran test -----------------------------------------------------------
df_geo <- left_join(df_results, location_gauges, by = "Site")
coords_mat <- as.matrix(df_geo[, c("Longitude", "Latitude")])
nb <- knn2nb(knearneigh(coords_mat, k = 4))   # k = 4 nearest neighbors
lw <- nb2listw(nb, style = "W")
moran.test(df_geo$kappa_adj, lw, alternative = "two.sided")
moran.test(df_geo$xi_adj, lw, alternative = "two.sided")
moran.test(df_geo$sigma_adj, lw, alternative = "two.sided")

ggplot(df_geo) +
  geom_point(aes(Longitude, Latitude, color = xi_adj), size = 4) +
  scale_color_viridis_c() +
  coord_fixed() 
  
z <- scale(df_geo$sigma_adj)[, 1]
lag_z <- lag.listw(lw, z)

df_moran <- data.frame(
  z = z,
  lag_z = lag_z,
  Site = df_geo$Site
)

ggplot(df_moran, aes(z, lag_z, label = Site)) +
  geom_point(size = 3) +
  geom_smooth(method = "lm", se = FALSE) +
  geom_hline(yintercept = 0, linetype = 2) +
  geom_vline(xintercept = 0, linetype = 2) +
  labs(
    x = "Standardized sigma",
    y = "Spatial lag of sigma",
    title = "Moran scatterplot"
  ) +
  theme_minimal()





ggplot(df_moran, aes(z, lag_z, label = Site)) +
  geom_point(size = 3) +
  geom_text_repel(size = 3) +
  geom_smooth(method = "lm", se = FALSE) +
  geom_hline(yintercept = 0, linetype = 2) +
  geom_vline(xintercept = 0, linetype = 2) +
  labs(
    x = "Standardized sigma",
    y = "Spatial lag of sigma"
  ) +
  theme_bw()


moran_sigma <- moran.test(df_geo$sigma_adj, lw, alternative = "two.sided")
moran_kappa <- moran.test(df_geo$kappa_adj, lw, alternative = "two.sided")
moran_xi    <- moran.test(df_geo$xi_adj, lw, alternative = "two.sided")


facet_labels <- c(
  sigma_adj = sprintf("sigma\nI = %.2f, p = %.3f",
                      moran_sigma$estimate[1],
                      moran_sigma$p.value),
  kappa_adj = sprintf("kappa\nI = %.2f, p = %.3f",
                      moran_kappa$estimate[1],
                      moran_kappa$p.value),
  xi_adj = sprintf("xi\nI = %.2f, p = %.3f",
                   moran_xi$estimate[1],
                   moran_xi$p.value)
)


# Plot the spatial distribution of the parameters with Moran's I values

df_map <- df_geo %>%
  select(Site, Longitude, Latitude, sigma_adj, kappa_adj, xi_adj) %>%
  pivot_longer(
    cols = c(sigma_adj, kappa_adj, xi_adj),
    names_to = "Parameter",
    values_to = "value"
  )

coords <- as.matrix(df_geo[, c("Longitude", "Latitude")])

segments <- do.call(
  rbind,
  lapply(seq_along(nb), function(i) {
    if(length(nb[[i]]) == 0) return(NULL)

    data.frame(
      x    = coords[i,1],
      y    = coords[i,2],
      xend = coords[nb[[i]],1],
      yend = coords[nb[[i]],2]
    )
  })
)



ggplot() +
  geom_segment(
    data = segments,
    aes(x = x, y = y, xend = xend, yend = yend),
    color = "grey50",
    alpha = 0.4,
    linewidth = 0.4
  ) +
  geom_point(
    data = df_geo,
    aes(Longitude, Latitude, fill = xi_adj),
    shape = 21,
    size = 5,
    color = "black"
  ) +
  geom_text_repel(
    data = df_geo,
    aes(Longitude, Latitude, label = Site),
    size = 3,
    max.overlaps = Inf
  ) +
  scale_fill_gradient2(
    name = expression(xi~adjacency),
    low = "darkred",
    mid = "white",
    high = btfgreen,
    midpoint = 0
  ) +
  coord_equal() +
  labs(
    title = sprintf(
      "Moran's I = %.2f, p = %.3f",
      moran_xi$estimate[1],
      moran_xi$p.value
    ),
    x = NULL,
    y = NULL
  ) +
  theme_bw() +
  theme(
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank(),
    panel.grid = element_blank()
  ) +
  theme(panel.border = element_blank())

# save plot
ggsave(paste0(save_path, "moran_xi_map.pdf"), width = 8, height = 6, dpi = 400)


# same for sigma_adj and kappa_adj

ggplot() +
  geom_segment(
    data = segments,
    aes(x = x, y = y, xend = xend, yend = yend),
    color = "grey50",
    alpha = 0.4,
    linewidth = 0.4
  ) +
  geom_point(
    data = df_geo,
    aes(Longitude, Latitude, fill = sigma_adj),
    shape = 21,
    size = 5,
    color = "black"
  ) +
  geom_text_repel(
    data = df_geo,
    aes(Longitude, Latitude, label = Site),
    size = 3,
    max.overlaps = Inf
  ) +
  scale_fill_gradient2(
    name = expression(sigma~adjacency),
    low = "darkred",
    mid = "white",
    high = btfgreen,
    midpoint = 0
  ) +
  coord_equal() +
  labs(
    title = sprintf(
      "Moran's I = %.2f, p = %.3f",
      moran_sigma$estimate[1],
      moran_sigma$p.value
    ),
    x = NULL,
    y = NULL
  ) +
  theme_bw() +
  theme(
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank(),
    panel.grid = element_blank()
  ) +
  theme(panel.border = element_blank())

ggsave(paste0(save_path, "moran_sigma_map.pdf"), width = 8, height = 6, dpi = 400)


ggplot() +
  geom_segment(
    data = segments,
    aes(x = x, y = y, xend = xend, yend = yend),
    color = "grey50",
    alpha = 0.4,
    linewidth = 0.4
  ) +
  geom_point(
    data = df_geo,
    aes(Longitude, Latitude, fill = kappa_adj),
    shape = 21,
    size = 5,
    color = "black"
  ) +
  geom_text_repel(
    data = df_geo,
    aes(Longitude, Latitude, label = Site),
    size = 3,
    max.overlaps = Inf
  ) +
  scale_fill_gradient2(
    name = expression(kappa~adjacency),
    low = "darkred",
    mid = "white",
    high = btfgreen,
    midpoint = 0
  ) +
  coord_equal() +
  labs(
    title = sprintf(
      "Moran's I = %.2f, p = %.3f",
      moran_kappa$estimate[1],
      moran_kappa$p.value
    ),
    x = NULL,
    y = NULL
  ) +
  theme_bw() +
  theme(
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank(),
    panel.grid = element_blank()
  ) +
  theme(panel.border = element_blank())

ggsave(paste0(save_path, "moran_kappa_map.pdf"), width = 8, height = 6, dpi = 400)



# Anderson-Darling test ---------------------------------------------------

rain_df <- rain_raw %>%
  dplyr::select(dates, all_of(sites)) %>%
  mutate(dates = as.Date(dates))

rain_sync <- rain_df 

rain_long_sync <- rain_sync %>%
  pivot_longer(-dates, names_to = "Site", values_to = "Rain") %>%
  filter(Rain > 0)

sites_sync <- unique(rain_long_sync$Site)

u <- quantile(rain_long_sync$Rain, 0.9)


cvm_results <- combn(sites_sync, 2, function(pair) {
  x1 <- rain_long_sync$Rain[rain_long_sync$Site == pair[1]]
  x2 <- rain_long_sync$Rain[rain_long_sync$Site == pair[2]]
  x1_tail <- x1[x1 > u]
  x2_tail <- x2[x2 > u]
  cvmtest <- twosamples::cvm_test(x1_tail, x2_tail, nboots = 10000)
  adtest <- twosamples::ad_test(x1_tail, x2_tail, nboots = 20000)

  data.frame(
    Site1 = pair[1],
    Site2 = pair[2],
    cvm_statistic = cvmtest[1],
    cvm_p.value = cvmtest[2],
    ad_statistic = adtest[1],
    ad_p.value = adtest[2]
  )
}, simplify = FALSE) |> bind_rows()


mean_cvm_ok <- mean(cvm_results$cvm_p.value > 0.05)
mean_cvm_ok
mean_ad_ok <- mean(cvm_results$ad_p.value > 0.01)

mean_ad_ok


# Do a fit over all sites together ------------------------------------------------
rain_no3gauges <- rain 

k_prec_censuring_all <- seq(0.22, 0.25, by = 0.01)
y_all <- as.data.frame(na.omit(as.vector(as.matrix(rain_no3gauges))))
site_result_all <- tryCatch({
  process_site(y = y_all, site_name = "All_sites", save_path = save_path, R = 1000, forced_cens = k_prec_censuring)
}, error = function(e) {
  warning("Failed for all sites")
  NULL
})

# Do Moran test on exceedances and proportion of 0 ---------------------------------
u <- quantile(rain[rain > 0], 0.95, na.rm = TRUE)
u
df_exceed <- rain_long_sync %>%
  filter(Rain > 0) %>%
  mutate(exceed = as.integer(Rain > u))

moran.test(df_exceed$exceed, lw, alternative = "two.sided")



df_site_exceed <- rain_long_sync %>%
  group_by(Site) %>%
  summarise(
    p_exceed95 = mean(Rain > u, na.rm = TRUE),
    .groups = "drop"
  )

df_geo_exceed <- df_geo %>%
  select(Site, Longitude, Latitude) %>%
  left_join(df_site_exceed, by = "Site")

moran.test(df_geo_exceed$p_exceed95, lw, alternative = "two.sided")

rain_long_sync <- rain_sync %>%
  pivot_longer(-dates, names_to = "Site", values_to = "Rain")


df_site_zero <- rain_long_sync %>%
  group_by(Site) %>%
  summarise(
    p0 = mean(Rain == 0, na.rm = TRUE),
    .groups = "drop"
  )

df_geo_zero <- df_geo %>%
  select(Site, Longitude, Latitude) %>%
  left_join(df_site_zero, by = "Site")

moran.test(df_geo_zero$p0, lw, alternative = "two.sided")

# Plot the spatial distribution of the parameters with Moran's I values
# for excesses and p0

df_map <- df_geo_exceed %>%
  select(Site, Longitude, Latitude, p_exceed95) %>%
  pivot_longer(
    cols = c(p_exceed95),
    names_to = "Parameter",
    values_to = "value"
  )

ggplot() +
  geom_segment(
    data = segments,
    aes(x = x, y = y, xend = xend, yend = yend),
    color = "grey50",
    alpha = 0.4,
    linewidth = 0.4
  ) +
  geom_point(
    data = df_geo_exceed,
    aes(Longitude, Latitude, fill = p_exceed95),
    shape = 21,
    size = 5,
    color = "black"
  ) +
  geom_text_repel(
    data = df_geo_exceed,
    aes(Longitude, Latitude, label = Site),
    size = 3,
    max.overlaps = Inf
  ) +
  scale_fill_gradient2(
    name = expression(p[exceed~95]),
    low = "darkred",
    mid = "white",
    high = btfgreen,
    midpoint = mean(df_geo_exceed$p_exceed95),
    limits = c(0.02, 0.07)
  ) +
  coord_equal() +
  labs(
    title = sprintf(
      "Moran's I: %.2f, p-value: %.3f",
      moran.test(df_geo_exceed$p_exceed95, lw)$estimate[1],
      moran.test(df_geo_exceed$p_exceed95, lw)$p.value
    ),
    x = NULL,
    y = NULL
  ) +
  theme_minimal() +
  theme(
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank(),
    panel.grid = element_blank()
  ) 

ggsave(paste0(save_path, "moran_p95_map.pdf"), width = 8, height = 6, dpi = 400)


df_map_zero <- df_geo_zero %>%
  select(Site, Longitude, Latitude, p0) %>%
  pivot_longer(
    cols = c(p0),
    names_to = "Parameter",
    values_to = "value"
  )

ggplot() +
  geom_segment(
    data = segments,
    aes(x = x, y = y, xend = xend, yend = yend),
    color = "grey50",
    alpha = 0.4,
    linewidth = 0.4
  ) +
  geom_point(
    data = df_geo_zero,
    aes(Longitude, Latitude, fill = p0),
    shape = 21,
    size = 5,
    color = "black"
  ) +
  geom_text_repel(
    data = df_geo_zero,
    aes(Longitude, Latitude, label = Site),
    size = 3,
    max.overlaps = Inf
  ) +
  scale_fill_gradient2(
    name = expression(p[0]),
    low = "darkred",
    mid = "white",
    high = btfgreen,
    midpoint = mean(df_geo_zero$p0),
    limits = c(0.98, 0.995)
  ) +
  coord_equal() +
  labs(
    title = sprintf(
      "Moran's I: %.2f, p-value: %.3f",
      moran.test(df_geo_zero$p0, lw)$estimate[1],
      moran.test(df_geo_zero$p0, lw)$p.value
    ),
    x = NULL,
    y = NULL
  ) +
  theme_minimal() +
  theme(
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank(),
    panel.grid = element_blank()
  )

# save plots
ggsave(paste0(save_path, "moran_p0_map.pdf"), width = 8, height = 6, dpi = 400)
