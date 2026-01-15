# remove objects from the environment
rm(list = ls())

# clear console
cat("\014")

# Load libraries and set theme
source("./script/load_libraries.R")

# Get all files in the folder "R"
functions_folder <- "./R"
files <- list.files(functions_folder, full.names = TRUE)
# # load all functions in files
invisible(lapply(files, function(f) source(f, echo = FALSE)))
library(latex2exp)

# get rain data from omsev
# filename_omsev <- paste0(data_folder,
#                          "omsev/omsev_5min/rain_mtp_5min_2019_2024.RData")

# load(filename_omsev)


# get rain data
filename_rain <- paste0(data_folder, "omsev/omsev_5min/rain_mtp_5min_2019_2024_cleaned.csv")
rain_omsev <- read.csv(filename_rain)
head(rain_omsev)

# get location of each rain gauge
filename_loc <- paste0(data_folder, "omsev/loc_rain_gauges.csv")
location_gauges <- read.csv(filename_loc)
location_gauges$Station <- c("iem", "mse", "poly", "um", "cefe", "cnrs",
                             "crbm", "archiw", "archie", "um35", "chu1",
                             "chu2", "chu3", "chu4", "chu5", "chu6", "chu7",
                             "cines", "brives", "hydro")

# save rain as csv
# filename_rain <- paste0(data_folder, "omsev/omsev_5min/rain_mtp_5min_2019_2024_cleaned.csv")
# # write.csv(rain, file = filename_rain, row.names = FALSE)

# rain_omsev <- read.csv(filename_rain)
# head(rain_omsev)


# rain_test <- rain_omsev
# rain_test$nb_sites_non_NA <- apply(rain_test[ , -1], 1, function(x) sum(!is.na(x)))
# date_2_sites <- rain_test$dates[which(rain_test$nb_sites_non_NA >= 2)[1]]
# date_3_sites <- rain_test$dates[which(rain_test$nb_sites_non_NA >= 3)[1]]
# date_4_sites <- rain_test$dates[which(rain_test$nb_sites_non_NA >= 4)[1]]
# date_5_sites <- rain_test$dates[which(rain_test$nb_sites_non_NA >= 5)[1]]
# # date_6_sites <- rain_test$dates[which(rain_test$nb_sites_non_NA >= 6)[1]]
# # date_7_sites <- rain_test$dates[which(rain_test$nb_sites_non_NA >= 7)[1]]


# # # begin in 2020-01-01
# rain_omsev <- rain_omsev[rain_omsev$dates >= date_5_sites, ]
# head(rain_omsev)





# put dates as rownames
rownames(rain_omsev) <- rain_omsev$dates
rain <- rain_omsev[-1] # remove dates column
# rain$mse
# Get distances matrix
dist_mat <- get_dist_mat(location_gauges)
df_dist <- reshape_distances(dist_mat)

# remove cines, hydro, brives
colnames(rain)
rain <- rain[, !(colnames(rain) %in% c("cines", "hydro", "brives"))]

location_gauges <- location_gauges[location_gauges$Station != "cines" &
                                   location_gauges$Station != "hydro" &
                                   location_gauges$Station != "brives", ]
dist_mat <- get_dist_mat(location_gauges)

df_dist <- reshape_distances(dist_mat)


library(ggplot2)

site_names <- location_gauges$Station

# change order of X and Y based on location_gauges
df_dist <- df_dist %>%
  mutate(
    X = factor(X, levels = rev(site_names)),
    Y = factor(Y, levels = site_names),
    diagonal = X == Y
  )

ggplot(df_dist, aes(x = X, y = Y, size = value)) +
  geom_point(data = df_dist %>% filter(!diagonal), alpha = 0.7, color = btfgreen) +
  scale_size_continuous(range = c(1, 10), breaks = c(100, 250, 500, 750, 1000, 1500)) +
  labs(title = "",
       x = "Site",
       y = "Site",
       size = "Distance") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.text.y = element_text(angle = 0, hjust = 1)
  )
# save plot
filename <- paste(im_folder, "rain/OMSEV/distance_matrix.png", sep = "")
ggsave(filename, width = 10, height = 5, units = "in", dpi = 300)

# get only those with distance > 500 m
thresholds_dist <- c(500, 600, 750, 1000, 1100, 1200, 1250, 1300, 1400, 1500, 2300)

count_pairs <- sapply(thresholds_dist, function(th) {
  nrow(df_dist %>% filter(value > th))
})

df_counts <- data.frame(
  threshold_m = thresholds_dist,
  nb_pairs = count_pairs
)
thresholds_dist <- seq(500, 1500, by = 50)

df_counts <- data.frame(
  threshold_m = thresholds_dist,
  nb_pairs = sapply(thresholds_dist, function(th) nrow(df_dist %>% filter(value >= th)))
)

ggplot(df_counts, aes(x = threshold_m, y = nb_pairs)) +
  geom_line(color = btfgreen, size = 1.2) +
  geom_point(color = "#668167", size = 3) +
  scale_x_continuous(breaks = thresholds_dist) +
  labs(
    title = "",
    x = "Distance threshold (m)",
    y = "Number of possible neighbor pairs"
  ) +
  theme_minimal() +
  # add a line at y = 50
  geom_hline(yintercept = 30, linetype = "dashed", color = "red") +
  annotate("text", x = 500, y = 31, label = "30 pairs", color = "red", hjust = 0)

# save plot
filename <- paste(im_folder, "rain/OMSEV/cumulative_distance_pairs_redline30.png", sep = "")
ggsave(filename, width = 10, height = 5, units = "in", dpi = 300)

################################################################################
# QUANTILE ---------------------------------------------------------------------
################################################################################

sites_names <- colnames(rain)
# site_combinations <- combn(sites_names, 2, simplify = FALSE)

# for (site in sites_names) {
#   filename <- paste0(im_folder, "mrlplot/omsev/", site, ".png")
#   png(filename, width = 10, height = 5, units = "in", res = 300)
#   par(mfrow = c(1, 1))
#   rain_pair <- rain[, c(site)]
#   # remove na
#   rain_pair <- rain_pair[complete.cases(rain_pair), ]
#   rain_pair_no0 <- rain_pair[rowSums(rain_pair) > 0, ]
#   mrlplot(rain_pair)

#   # get corresponding quantile for threshold u

#   dev.off()
# }

# library(ggplot2)


# sites_names <- colnames(rain)
# site_combinations <- combn(sites_names, 2, simplify = FALSE)
# q <- 0.955
# excess_counts_spa <- data.frame(site1 = character(), site2 = character(), 
#                       n_excesses = integer(), stringsAsFactors = FALSE)
# for (pair in site_combinations) {
#   site1 <- pair[1]
#   site2 <- pair[2]
#   filename <- paste0(im_folder, "chiplot/omsev/", site1, "_", site2, ".png")
#   png(filename, width = 10, height = 5, units = "in", res = 300)

#   par(mfrow = c(1, 1))
#   rain_pair <- rain[, c(site1, site2)]
#   # remove na
#   rain_pair <- rain_pair[complete.cases(rain_pair), ]
#   colnames(rain_pair)
#   rain_pair_no0 <- rain_pair[rowSums(rain_pair) > 0, ]
#   chiplot(rain_pair_no0, xlim = c(0.9, 1), ylim1 = c(0, 1), which = 1,
#           qlim = c(0.9, 0.995), main1 = "Without zeros")
#   abline(v = q, col = "red", lty = 2)
#   dev.off()
#   n <- nrow(rain_pair_no0)
#   data_unif <- cbind(rank(rain_pair_no0[, 1]) / (n + 1),
#                     rank(rain_pair_no0[, 2]) / (n + 1))

#   count_excesses <- sum(data_unif[, 1] > q & data_unif[, 2] > q)
#   excess_counts_spa <- rbind(excess_counts_spa, data.frame(site1 = site1,
#                                                    site2 = site2,
#                                                    n_excesses = count_excesses,
#                                                    stringsAsFactors = FALSE))
 
# }

# excess_save <- excess_counts_spa
# # get latex table for excess counts
# excess_counts_spa <- excess_counts_spa %>%
#   mutate(site1 = factor(site1, levels = sites_names),
#          site2 = factor(site2, levels = sites_names)) %>%
#   arrange(site1, site2)


# ggplot(excess_counts_spa %>% filter(n_excesses > 0),
#        aes(x = site2, y = site1, size = n_excesses)) +
#   geom_point(alpha = 0.7, color = btfgreen) +
#   geom_text(aes(label = ifelse(n_excesses < 40, n_excesses, "")),
#             size = 3, vjust = -1) +
#   scale_size_continuous(range = c(1, 10)) +
#   theme_minimal() +
#   labs(title = "", size = "Number of excesses") +
#   theme(axis.text.x = element_text(angle = 45, hjust = 1))
# # save plot
# filename <- paste(im_folder, "WLSE/omsev/spatial/excess_counts_spatial_", 
#                     q, ".png", sep = "")
# ggsave(filename, width = 20, height = 15, units = "cm")




# templag <- 0:10
# q <- 0.94
# excess_counts_temp <- data.frame(site1 = character(), tau = integer(), 
#                       n_excesses = integer(), stringsAsFactors = FALSE)
# for (site in sites_names) {
#   for (tau in templag) {
#     site1 <- site
#     site2 <- paste0(site, "_lag", tau)
#     # create a lagged version of the site
#     rain_lagged <- cbind(rain[[site]][1:(nrow(rain) - tau)],
#                          rain[[site]][(tau + 1):nrow(rain)])
#     colnames(rain_lagged) <- c(site1, site2)
#     rain_pair <- as.data.frame(rain_lagged)
#     # remove na
#     rain_pair <- rain_pair[complete.cases(rain_pair), ]
#     rain_pair_no0 <- rain_pair[rowSums(rain_pair) > 0, ]
#     filename <- paste0(im_folder, "chiplot/omsev/", site1, "_", site2, ".png")
#     png(filename, width = 10, height = 5, units = "in", res = 300)

#     par(mfrow = c(1, 1))
#     rain_nolag <- rain[[site1]][1:(length(rain[[site1]]) - tau)]
#     rain_lag <- rain[[site1]][(tau + 1):length(rain[[site1]])]
#     rain_pair <- cbind(rain_nolag, rain_lag)
#     # remove na
#     rain_pair <- rain_pair[complete.cases(rain_pair), ]
#     colnames(rain_pair)
#     rain_pair_no0 <- rain_pair[rowSums(rain_pair) > 0, ]
#     chiplot(rain_pair_no0, xlim = c(0.9, 1), ylim1 = c(0, 1), which = 1,
#             qlim = c(0.9, 0.98), main1 = "Without zeros")
#     abline(v = q, col = "red", lty = 2)
#     dev.off()
#     n <- nrow(rain_pair_no0)
#     data_unif <- cbind(rank(rain_pair_no0[, 1]) / (n + 1),
#                       rank(rain_pair_no0[, 2]) / (n + 1))

#     count_excesses <- sum(data_unif[, 1] > q & data_unif[, 2] > q)
#     excess_counts_temp <- rbind(excess_counts_temp, data.frame(site = site,
#                                                    tau = tau,
#                                                    n_excesses = count_excesses,
#                                                    stringsAsFactors = FALSE))
#   }
 
# }

# library(ggplot2)
# library(scales)

# ymin <- min(excess_counts_temp$n_excesses, na.rm = TRUE)

# ggplot(excess_counts_temp, aes(x = tau, y = n_excesses, color = site)) +
#   geom_line() +
#   geom_point() +
#   theme_minimal() +
#   scale_y_continuous(
#     limits = c(0, NA),
#     breaks = pretty(c(0, excess_counts_temp$n_excesses), n = 15)
#   ) +
#   scale_x_continuous(
#     breaks = templag
#   ) +
#   labs(
#     title = "",
#     x = "Temporal lag",
#     y = "Number of joint excesses",
#     color = "Site"
#   )

################################################################################
# WLSE results -----------------------------------------------------------------
################################################################################
# foldername <- paste0(data_folder, "omsev/WLSE/")
# df_spa_result <- read.csv(paste0(foldername, "wlse_spa_summary.csv"))
# df_temp_result <- read.csv(paste0(foldername, "wlse_temp_summary.csv"))
# df_result_all <- merge(df_spa_result, df_temp_result, by = NULL)

# # select one row
# df_result <- df_result_all[df_result_all$q_spa == 0.97 &
#                              df_result_all$q_temp == 0.95, ]

# beta1 <- df_result$beta1
# beta2 <- df_result$beta2
# alpha1 <- df_result$alpha1
# alpha2 <- df_result$alpha2

# # estimates
# wlse_omsev <- c(beta1, beta2, alpha1, alpha2) # m / 5 min

################################################################################
# VARIOGRAM --------------------------------------------------------------------
################################################################################
q <- 0.95 # quantile
set_st_excess <- get_spatiotemp_excess(rain, quantile = q, remove_zeros = TRUE)

sites_coords <- location_gauges[, c("Longitude", "Latitude")]

rownames(sites_coords) <- location_gauges$Station
sites_coords_sf <- st_as_sf(sites_coords, coords = c("Longitude", "Latitude"),
                            crs = 4326)

sites_coords_sf <- st_transform(sites_coords_sf, crs = 2154)

coords_m <- st_coordinates(sites_coords_sf)

grid_coords_km <- sites_coords
grid_coords_m <- sites_coords
grid_coords_m$x_m <- (coords_m[, "X"] - min(coords_m[, "X"]))
grid_coords_m$y_m <- (coords_m[, "Y"] - min(coords_m[, "Y"]))
grid_coords_km$x_km <- (coords_m[, "X"] - min(coords_m[, "X"])) / 1000
grid_coords_km$y_km <- (coords_m[, "Y"] - min(coords_m[, "Y"]))  / 1000

# get distance matrix
grid_coords_m <- grid_coords_m[, c("x_m", "y_m")]
grid_coords_km <- grid_coords_km[, c("x_km", "y_km")]
colnames(grid_coords_m) <- c("Longitude", "Latitude")
colnames(grid_coords_km) <- c("Longitude", "Latitude")
dist_mat <- get_dist_mat(grid_coords_m, latlon = FALSE)

# Spatial chi
df_dist <- reshape_distances(dist_mat)
df_dist_km <- df_dist
df_dist_km$value <- df_dist$value / 1000

stopifnot(nrow(grid_coords_m) == ncol(rain))
rownames(grid_coords_m) <- colnames(rain)



library(data.table)
library(dplyr)
library(ggplot2)

dmins <- seq(200, 2500, by = 100)

episode_pair_stats <- function(sel, sites_coords, dmin_m,
                               delta_steps, step_minutes = 5,
                               latlon = FALSE, beta = 0) {

  if (nrow(sel) < 2) {
    return(list(
      n_episodes = nrow(sel),
      n_pairs = 0,
      n_pairs_dtlt_delta = 0,
      p_close_space_500m = NA_real_,
      p_close_space_1200m = NA_real_,
      p_close_space_500m_given_dtlt_delta = NA_real_,
      p_conflict_dmin_delta = NA_real_,
      spatial_d10 = NA_real_,
      temporal_d10 = NA_real_
    ))
  }

  dist_matrix <- get_dist_mat(sites_coords, latlon = latlon)
  if (latlon) dist_matrix <- dist_matrix / 1000  # attention unités si latlon

  s <- sel$s0
  t <- sel$t0

  idx <- combn(seq_along(s), 2)
  s1 <- s[idx[1,]]; s2 <- s[idx[2,]]
  t1 <- t[idx[1,]]; t2 <- t[idx[2,]]

  spatial_d <- dist_matrix[cbind(s1, s2)]

  # temporal distance en "steps" et en minutes
  temporal_steps <- abs(t1 - t2)
  temporal_min <- temporal_steps

  # même fenêtre que ta règle de sélection
  delta_eff_steps <- delta_steps + 2 * beta
  is_dt_close <- temporal_steps < delta_eff_steps

  # indicateurs globaux (toutes paires)
  p_close_500_all <- mean(spatial_d < 500, na.rm = TRUE)
  p_close_1200_all <- mean(spatial_d < 1200, na.rm = TRUE)

  # indicateur conditionnel aux paires temporellement proches
  n_pairs_dt <- sum(is_dt_close, na.rm = TRUE)
  p_close_500_given_dt <- if (n_pairs_dt == 0) NA_real_ else
    mean((spatial_d < 500)[is_dt_close], na.rm = TRUE)

  # proportion de paires "conflit" au sens de ta règle (si on ne déclusterait pas)
  p_conflict <- mean((spatial_d < dmin_m) & is_dt_close, na.rm = TRUE)

  list(
    n_episodes = nrow(sel),
    n_pairs = length(spatial_d),
    n_pairs_dtlt_delta = n_pairs_dt,
    p_close_space_500m  = p_close_500_all,
    p_close_space_1200m = p_close_1200_all,
    p_close_space_500m_given_dtlt_delta = p_close_500_given_dt,
    p_conflict_dmin_delta = p_conflict,
    spatial_d10  = as.numeric(quantile(spatial_d, 0.10, na.rm = TRUE)),
    temporal_d10 = as.numeric(quantile(temporal_min, 0.10, na.rm = TRUE))
  )
}


episode_size <- 12
set_st_excess <- get_spatiotemp_excess(rain, quantile = q, remove_zeros = TRUE)
res <- lapply(dmins, function(dm) {
  sel <- get_s0t0_pairs(
    sites_coords = grid_coords_m,      # en mètres + rownames OK
    data = rain,
    min_spatial_dist = dm,
    episode_size = episode_size,
    set_st_excess = set_st_excess,
    n_max_episodes = 10000,
    latlon = FALSE,
    beta = 0
  )

  st <- episode_pair_stats(
    sel = sel,
    sites_coords = grid_coords_m,
    dmin_m = dm,
    delta_steps = episode_size,
    step_minutes = 5,
    latlon = FALSE,
    beta = 0
  )

  data.frame(
    dmin_m = dm,
    n_episodes = st$n_episodes,
    n_pairs_dtlt_delta = st$n_pairs_dtlt_delta,
    p_close_space_500m_all = st$p_close_space_500m,
    p_close_space_500m_given_dtlt_delta = st$p_close_space_500m_given_dtlt_delta,
    p_conflict_dmin_delta = st$p_conflict_dmin_delta,
    spatial_d10_m = st$spatial_d10,
    temporal_d10_min = st$temporal_d10
  )
})

df_tradeoff <- bind_rows(res)

pA <- ggplot(df_tradeoff, aes(x = dmin_m, y = n_episodes)) +
  geom_line(size = 1.1, color = btfgreen) +
  geom_point(size = 2, color = btfgreen) +
  geom_vline(xintercept = 1200, linetype = "dashed", color = "red") +
  geom_hline(yintercept = 300, linetype = "dashed", color = "#5a4b4b") +
  theme_minimal() +
  labs(
    x = expression(d[min]~"(m)"),
    y = "Number of selected episodes"
  )

foldername <- paste0(im_folder, "optim/omsev/choice_config/")
if (!dir.exists(foldername)) {
  dir.create(foldername, recursive = TRUE)
}
filename <- paste0(foldername, "tradeoff_dmin_episodes.png")
ggsave(filename, plot = pA, width = 7, height = 5, units = "in", dpi = 300)


library(dplyr)
library(ggplot2)

dmin_fixed <- 1600
delta_grid <- seq(6, 48, by = 2)   # en "steps" de 5 minutes

set_st_excess <- get_spatiotemp_excess(rain, quantile = q, remove_zeros = TRUE)

res_delta <- lapply(delta_grid, function(delta_steps) {

  sel <- get_s0t0_pairs(
    sites_coords = grid_coords_m,
    data = rain,
    min_spatial_dist = dmin_fixed,
    episode_size = delta_steps,
    set_st_excess = set_st_excess,
    n_max_episodes = 10000,
    latlon = FALSE,
    beta = 0
  )

  st <- episode_pair_stats(
    sel = sel,
    sites_coords = grid_coords_m,
    dmin_m = dmin_fixed,
    delta_steps = delta_steps,
    step_minutes = 5,
    latlon = FALSE,
    beta = 0
  )

  data.frame(
    delta_steps = delta_steps,
    delta_minutes = 5 * delta_steps,
    n_episodes = st$n_episodes,
    n_pairs_dtlt_delta = st$n_pairs_dtlt_delta,
    p_close_space_500m_given_dtlt_delta = st$p_close_space_500m_given_dtlt_delta,
    p_conflict_dmin_delta = st$p_conflict_dmin_delta,
    spatial_d10_m = st$spatial_d10,
    temporal_d10_min = st$temporal_d10
  )
})

df_tradeoff_delta <- bind_rows(res_delta)

p_deltaA <- ggplot(df_tradeoff_delta, aes(x = delta_minutes, y = n_episodes)) +
  geom_line(size = 1.1, color = btfgreen) +
  geom_point(size = 2, color = btfgreen) +
  geom_vline(xintercept = 60, linetype = "dashed", color = "red") +  # 12*5=60
  theme_minimal() +
  labs(
    x = expression(delta~"(minutes)"),
    y = "Number of selected episodes"
  )

print(p_deltaA)

ggsave(paste0(foldername, "tradeoff_delta_episodes_dmin1600.png"),
       plot = p_deltaA, width = 7, height = 5, units = "in", dpi = 300)

################################################################################

# in rain remove when all data are NA
q <- 0.95 # quantile
set_st_excess <- get_spatiotemp_excess(rain, quantile = q, remove_zeros = TRUE)

# verify that the excess is above the threshold
# get list of sites and times
list_s <- set_st_excess$list_s
list_t <- set_st_excess$list_t
list_u <- set_st_excess$list_u

min_u <- min(unlist(list_u))
max_u <- max(unlist(list_u))

# Spatio-temporal neighborhood parameters
min_spatial_dist <- 1600 # m
delta <- 12 # in * 5 min
episode_size <- delta # size of the episode
sites_coords <- location_gauges[, c("Longitude", "Latitude")]
tail(rain)
s0t0_set <- get_s0t0_pairs(grid_coords_m, rain,
                            min_spatial_dist = min_spatial_dist,
                            episode_size = episode_size,
                            set_st_excess = set_st_excess,
                            n_max_episodes = 10000,
                            latlon = FALSE)

selected_points <- s0t0_set
selected_points[12,]
length(selected_points$s0) # number of selected points

selected_points <- selected_points %>%
  mutate(t0_date = as.POSIXct(t0_date, format = "%Y-%m-%d %H:%M:%S", tz = "UTC"))


# Assuming s0t0_set is a data.table or data.frame
library(data.table)

# Make sure rain is in matrix form
# site names in columns, time in rows
# column names must match `s0t0_set$s0`
stopifnot(all(s0t0_set$s0 %in% colnames(rain)))

# For each (s0, t0, u_s0), check if rain[t0, s0] > u_s0
excess_check_s0t0 <- s0t0_set[, {
  rain_val <- rain[t0, s0]
  is_excess <- rain_val > u_s0
  list(rain_value = rain_val, is_excess = is_excess)
}, by = .(s0, t0, u_s0)]

# Check how many are not true exceedances
invalid_exceedances <- excess_check_s0t0[is_excess == FALSE]

# Summary
cat("Total s0t0 pairs:", nrow(s0t0_set), "\n")
cat("Invalid exceedances (rain ≤ threshold):", nrow(invalid_exceedances), "\n")

# print all u_s0 values
cat("Unique u_s0 values:", length(unique(s0t0_set$u_s0)), "\n")
cat("Minimum u_s0 value:", min(s0t0_set$u_s0), "\n")
cat("Maximum u_s0 value:", max(s0t0_set$u_s0), "\n")

# check that for all s0, t0 we have an excess above corresponding threshold
for (i in 1:length(selected_points$s0)) {
  s0 <- selected_points$s0[i]
  t0 <- selected_points$t0[i]
  u_s0 <- selected_points$u_s0[i]
  # check that the excess is above the threshold
  if (rain[t0, s0] <= u_s0) {
    stop(paste("Excess is not above threshold for s0 =", s0, "and t0 =", t0))
  }
}

# Histogram of dates t0 of selected points
df_hist_t0 <- data.frame(t0_date = as.POSIXct(selected_points$t0_date, format = "%Y-%m-%d %H:%M:%S", tz = "UTC"))
df_hist_t0$month <- format(df_hist_t0$t0_date, "%Y-%m")
selected_months <- unique(df_hist_t0$month)[seq(1, length(unique(df_hist_t0$month)), by = 3)]  # tous les 3 mois

ggplot(df_hist_t0, aes(x = month)) +
  geom_bar(fill = btfgreen, alpha = 0.5, color = "#5f5d5d") +
  btf_theme +
  xlab("Month of episode start") +
  ylab("Count") +
  scale_x_discrete(breaks = selected_months) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

filename <- paste(im_folder, "optim/omsev/episodes_histogram/t0_histogram_q",
                  q * 100, "_delta", delta, "_dmin", min_spatial_dist,
                  ".png", sep = "")

ggsave(filename, width = 20, height = 15, units = "cm")


# Threshold histogram
df_threshold <- data.frame(u_s0 = selected_points$u_s0)
unique(df_threshold$u_s0)
breaks <- seq(floor(min(df_threshold$u_s0)), ceiling(max(df_threshold$u_s0)), by = 0.1)

ggplot(df_threshold, aes(x = u_s0)) +
  geom_histogram(breaks = breaks, fill = btfgreen, alpha = 0.5, color = "#5f5d5d") +
  btf_theme +
  xlab(TeX(paste0("Threshold for quantile $q = ", q, "$"))) +
  ylab("Count")
filename <- paste(im_folder, "optim/omsev/threshold_histogram_q",
                  q * 100, "_delta", delta, "_dmin", min_spatial_dist,
                  ".png", sep = "")
ggsave(filename, width = 20, height = 15, units = "cm")

n_episodes <- length(selected_points$s0)
t0_list <- selected_points$t0
s0_list <- selected_points$s0

list_episodes_points <- get_extreme_episodes(selected_points, rain,
                                     episode_size = episode_size, unif = FALSE)

list_episodes <- list_episodes_points$episodes
episode <- list_episodes[[1]]
head(episode)
head(selected_points$t0_date)
t0_1 <- selected_points$t0_date[1]
# check excess for the episode
s0_1 <- selected_points$s0[1]
u_s0_1 <- selected_points$u_s0[1]
episode[1, s0_1] > u_s0_1 # should be

rain[as.POSIXct(rownames(rain), format = "%Y-%m-%d %H:%M:%S", tz = "UTC") == t0_1, s0_1] > u_s0_1 # should be true

selected_points$t0_date <- as.POSIXct(selected_points$t0_date, format = "%Y-%m-%d %H:%M:%S", tz = "UTC")


# save 5min dates and hour
library(lubridate)

# Round to the next hour
selected_points$t0_date_rounded <- ceiling_date(selected_points$t0_date, "hour")

datetimes <- unique(selected_points$t0_date)


datetimes_hour <- unique(selected_points$t0_date_rounded)


# save datetime list to csv
datetime_filename <- paste(data_folder, "/omsev/t0_episodes_q", q * 100,
                           "_delta", delta, "_dmin", min_spatial_dist,
                           ".csv", sep = "")
write.csv(data.frame(t0_date = datetimes_hour), datetime_filename, row.names = FALSE)


# save datetime list to csv
datetime_filename <- paste(data_folder, "/omsev/t0_5min_episodes_q", q * 100,
                           "_delta", delta, "_dmin", min_spatial_dist,
                           ".csv", sep = "")
write.csv(data.frame(t0_date = datetimes), datetime_filename, row.names = FALSE)

# get comephore advection for each episode
adv_omsev_filename <- paste(data_folder, "/omsev/adv_estim/bary_omsev/episode_advection_q",
                          q * 100, "_delta", delta, "_dmin", min_spatial_dist,
                          ".csv", sep = "")
adv_omsev_df <- read.csv(adv_omsev_filename, sep = ",")
head(adv_omsev_df)
adv_filename <- paste(data_folder, "/omsev/adv_estim/combined_comephore_omsev/episode_advection_q",
                      q * 100, "_delta", delta, "_dmin", min_spatial_dist,
                      ".csv", sep = "")
adv_df <- read.csv(adv_filename, sep = ",")
head(adv_df)

# number of 0 advection
nrow(adv_df[adv_df$mean_dx_kmh == 0 & adv_df$mean_dy_kmh == 0, ])
nrow(adv_df)

adv_df <- adv_df %>% mutate(t0 = ymd_hms(t0_omsev))

library(dplyr)
library(lubridate)
library(fuzzyjoin)

# Sélection des lignes où advection = 0
adv_zero <- adv_df %>%
  filter(vx_final == 0 & vy_final == 0)

# Jointure approx sur ±1h (choisit le plus proche)
adv_zero_match <- difference_inner_join(
  adv_zero, adv_omsev_df,
  by = "t0",
  max_dist = dminutes(60),
  distance_col = "time_diff"
) %>%
  group_by(t0.x) %>%
  slice_min(time_diff) %>%
  ungroup()

# Remettre les valeurs omsev dans adv_df
adv_df <- adv_df %>%
  left_join(
    adv_zero_match %>%
      select(t0 = t0.x,
             dx_omsev = vx_final,
             dy_omsev = vy_final),
    by = "t0"
  ) %>%
  mutate(mean_dx_kmh = ifelse(mean_dx_kmh == 0 & !is.na(dx_omsev), dx_omsev, mean_dx_kmh),
         mean_dy_kmh = ifelse(mean_dy_kmh == 0 & !is.na(dy_omsev), dy_omsev, mean_dy_kmh)) %>%
  select(-dx_omsev, -dy_omsev)


nrow(adv_df[adv_df$mean_dx_kmh == 0 & adv_df$mean_dy_kmh == 0, ])


head(adv_df)



# get only matching episodes from selected_points
matching_indices <- match(selected_points$t0_date_rounded, adv_df$t0)
# remove NA indices
matching_indices <- matching_indices[!is.na(matching_indices)]
# get only matching rows
adv_df <- adv_df[matching_indices, ]
rownames(adv_df) <- NULL  # reset row names

selected_episodes <- selected_points
selected_episodes$adv_x <- rep(NA, nrow(selected_episodes))
selected_episodes$adv_y <- rep(NA, nrow(selected_episodes))

# get adv values for each episode according to the t0_date
for (i in 1:nrow(selected_episodes)) {
  t0_date <- selected_episodes$t0_date_rounded[i]
  adv_row <- adv_df[adv_df$t0 == t0_date, ]
  if (nrow(adv_row) == 0) {
    print(i)
    print(paste("No advection data found for t0_date =", t0_date))
  } else {
    # if there are multiple rows, take the first one
    adv_row <- adv_row[1, ]
    adv_x <- adv_row$mean_dx_kmh
    adv_y <- adv_row$mean_dy_kmh
    selected_episodes$adv_x[i] <- adv_x
    selected_episodes$adv_y[i] <- adv_y
  }
}

head(selected_episodes)
tail(selected_episodes)

V_episodes <- data.frame(
  v_x = selected_episodes$adv_x,
  v_y = selected_episodes$adv_y
)

colnames(V_episodes) <- c("vx", "vy")

tau_vect <- 0:10
thresholds_by_site <- apply(rain, 2, function(col) {
  col <- col[!is.na(col) & col > 0]
  if (length(col) < 30) return(NA)
  quantile(col, probs = q)
})

# plot wind df vectors
library(ggplot2)
library(grid)  # pour arrow()

wind_df_plot <- ggplot(selected_episodes, aes(x = 0, y = 0)) +
  geom_segment(aes(xend = adv_x, yend = adv_y), 
               arrow = arrow(length = unit(0.2, "cm")),
               color = btfgreen, size = 0.5) +
  geom_point(aes(x = adv_x, y = adv_y), color = btfgreen, size = 1) +
  btf_theme +
  coord_equal() +
  xlab("Advection in x (km/h)") +
  ylab("Advection in y (km/h)") +
  ggtitle("") +
  theme(plot.title = element_text(hjust = 0.5)) 

wind_df_plot

# save wind df plot
filename <- paste(im_folder, "optim/omsev/advection_com_emp_plot_q",
                  q * 100, "_delta", delta, "_dmin", min_spatial_dist,
                  ".png", sep = "")
ggsave(filename, plot = wind_df_plot, width = 20, height = 15, units = "cm",
       dpi = 600, device = "png")

# number of 0 advection
nrow(selected_episodes[selected_episodes$adv_x == 0 & selected_episodes$adv_y == 0, ])
nrow(selected_episodes)

selected_episodes_nona <- selected_episodes[!is.na(selected_episodes$adv_x) & !is.na(selected_episodes$adv_y), ]

# remove those with 0 advection
# selected_points_adv <- selected_points_nona[
  # selected_points_nona$adv_x != 0 | selected_points_nona$adv_y != 0, ]
s0_list <- selected_episodes_nona$s0
list_episodes_points <- get_extreme_episodes(selected_episodes_nona, rain,
                              episode_size = episode_size, unif = FALSE)

list_episodes <- list_episodes_points$episodes
episode <- list_episodes[[1]]
u0_list <- selected_episodes_nona$u_s0
tmax <- max(tau_vect)
df_coords <- as.data.frame(grid_coords_km)
# Compute the lags and excesses for each conditional point
list_results <- mclapply(1:length(s0_list), function(i) {
  s0 <- s0_list[i]
  col_s0 <- which(colnames(rain) == s0)
  s0_coords <- df_coords[col_s0, ]
  # t0 <- t0_list[i]
  episode <- list_episodes[[i]]
  ind_t0_ep <- 0 # index of t0 in the episode
  lags <- get_conditional_lag_vectors(df_coords, s0_coords, ind_t0_ep,
                                tau_vect, latlon = FALSE)
  # hnorm is in meters
  # lags$hnorm <- lags$hnorm

  # excesses <- empirical_excesses_rpar(episode, q, lags, t0 = ind_t0_ep)
  u <- u0_list[i]
  excesses <- empirical_excesses_rpar(episode, threshold = u,
                                  df_lags = lags, t0 = ind_t0_ep)
  # tau is in 5 minutes
  lags$tau <- lags$tau * 5 / 60 # convert to hours
  list(lags = lags, excesses = excesses)
}, mc.cores = detectCores() - 1)

list_lags <- lapply(list_results, `[[`, "lags")
list_excesses <- lapply(list_results, `[[`, "excesses")

s0 <- s0_list[1]
col_s0 <- which(colnames(rain) == s0)
s0_coords <- df_coords[col_s0, ]
excesses <- list_excesses[[2]]
sum(excesses$kij)
df_lags <- list_lags[[1]]
tail(df_lags)

# plot sum of kij
kij_values <- sapply(list_excesses, function(excess) sum(excess$kij))
hist(kij_values, breaks = 30, main = "Histogram of sum of kij",
     xlab = "Sum of kij", col = btfgreen, border = "#5f5d5d")



filename_com_res <- paste(data_folder, "/comephore/optim_results/lalpha/free_eta/combined_optim_results.csv", sep = "")
com_results <- read.csv(filename_com_res)

# get estimates from comephore
params_com <- com_results[com_results$q == q*100 &
                           com_results$delta == 30 &
                           com_results$dmin == 5, ]

# params_com <- params_com[2,]
init_params_com <- c(params_com$beta1, params_com$beta2,
                     params_com$alpha1, params_com$alpha2,
                     params_com$eta1, params_com$eta2)
# check for na in adv and wind
V_episodes <- V_episodes[!is.na(V_episodes$vx) & !is.na(V_episodes$vy), ]

# plot V_episodes
# head(V_episodes)
# plot(V_episodes$vx, V_episodes$vy, xlab = "vx (km/h)", ylab = "vy (km/h)",
#      main = "Advection vectors from episodes", col = btfgreen, pch = 16)
# abline(h = 0, v = 0, col = "gray", lty = 2)
# With etas
hmax <- max(dist_mat) / 1000 # convert to km
df_lags <- list_lags[[1]]
head(df_lags)

# remove 0 advection episodes from all lists
# nonzero_indices <- which(V_episodes$vx != 0 | V_episodes$vy != 0)
# V_episodes_no0 <- V_episodes[nonzero_indices, ]
# list_lags_no0 <- list_lags[nonzero_indices]
# list_excesses_no0 <- list_excesses[nonzero_indices]
# list_episodes_no0 <- list_episodes[nonzero_

# remove adv < 0.1 km/h
adv_magnitudes <- sqrt(V_episodes$vx^2 + V_episodes$vy^2)
#plot histogram of adv magnitudes
hist(adv_magnitudes, breaks = 30, main = "Histogram of advection magnitudes",
     xlab = "Advection magnitude (km/h)", col = btfgreen, border = "#5f5d5d")
nonzero_indices <- which(adv_magnitudes >= 0.1)
V_episodes_no0 <- V_episodes[nonzero_indices, ]
list_lags_no0 <- list_lags[nonzero_indices]
list_excesses_no0 <- list_excesses[nonzero_indices]
list_episodes_no0 <- list_episodes[nonzero_indices]


convert_params <- function(beta1, beta2, alpha1, alpha2, 
                           eta1, eta2, 
                           c_x = 1, c_t = 1) {
  beta1_new <- beta1 / (c_x^alpha1)
  beta2_new <- beta2 / (c_t^alpha2)
  
  eta1_new <- eta1 * (c_t / c_x)^2  # depending on the model
  
  list(beta1 = beta1_new, 
       beta2 = beta2_new, 
       eta1 = eta1_new, 
       eta2 = eta2)
}



init_params_com <- c(params_com$beta1, params_com$beta2,
                     params_com$alpha1, params_com$alpha2,
                     params_com$eta1, params_com$eta2)
# 0 adv
V_episodes_adv0 <- data.frame(vx = rep(0, nrow(V_episodes)), vy = rep(0, nrow(V_episodes)))
# change_init <- c(0.2, 1, 1, 1, 8, 2)

init_params <- c(0.38, 0.85, 0.035, 0.69, 0.1, 1)
estimates_adv0 <- c(0.50886955, 4.87047484, 0.02650187, 0.70188455, 0, 1)

result <- optim(par = init_params_com[1:4], fn = neg_ll_composite_fixed_eta,
        list_lags = list_lags, list_episodes = list_episodes,
        list_excesses = list_excesses, hmax = hmax,
        wind_df = V_episodes,
        latlon = FALSE,
        fixed_eta1 = init_params_com[5],
        fixed_eta2 = init_params_com[6],
        distance = "lalpha",
        method = "L-BFGS-B",
        lower = c(1e-08, 1e-08, 1e-08, 1e-08, -1e-08, 1e-08),
        upper = c(10, 10, 1.999, 1.999, 10, 10),
        control = list(maxit = 10000,
                      trace = 1),
        hessian = F)

result





# Different initial parameter sets
init_list <- list(
  init_params_com,              # comephore estimates
  c(0.5, 5, 0.2, 0.7, 2, 4),   # close to first solution
  c(1, 1, 0.5, 0.5, 1, 1),     # more balanced start
  c(0.1, 0.1, 0.1, 0.1, 5, 5), # skewed towards small betas/alphas, large etas
  c(2, 8, 0.3, 1.2, 3, 0.5)    # different spread
)

# Run optim for each set of initials
results <- lapply(init_list, function(init_params_com) {
  optim(par = init_params_com, fn = neg_ll_composite,
        list_lags = list_lags, list_episodes = list_episodes,
        list_excesses = list_excesses, hmax = NA,
        wind_df = V_episodes,
        latlon = FALSE,
        distance = "lalpha",
        method = "L-BFGS-B",
        lower = c(1e-08, 1e-08, 1e-08, 1e-08, 1e-08, 1e-08),
        upper = c(10, 10, 1.999, 1.999, 10, 10),
        control = list(maxit = 10000, trace = 0),
        hessian = FALSE)
})

# Summarize results in a table: initial values and estimates
results_table <- data.frame(
  init_beta1 = sapply(init_list, function(x) x[1]),
  init_beta2 = sapply(init_list, function(x) x[2]),
  init_alpha1 = sapply(init_list, function(x) x[3]),
  init_alpha2 = sapply(init_list, function(x) x[4]),
  init_eta1 = sapply(init_list, function(x) x[5]),
  init_eta2 = sapply(init_list, function(x) x[6]),
  est_beta1 = sapply(results, function(res) res$par[1]),
  est_beta2 = sapply(results, function(res) res$par[2]),
  est_alpha1 = sapply(results, function(res) res$par[3]),
  est_alpha2 = sapply(results, function(res) res$par[4]),
  est_eta1 = sapply(results, function(res) res$par[5]),
  est_eta2 = sapply(results, function(res) res$par[6]),
  neg_loglik = sapply(results, function(res) res$value)
)
print(results_table)

params <- c("beta1","beta2","alpha1","alpha2","eta1","eta2")
plot_data <- data.frame(init = unlist(results_table[, paste0("init_", params)]),
                        est  = unlist(results_table[, paste0("est_", params)]),
                        param = rep(params, each = nrow(results_table)))

# Box plot with params labels in latex
library(ggplot2)
library(latex2exp)
# Define parameter labels in LaTeX
param_labels <- c(
  beta1 = TeX("$\\beta_1$"),
  beta2 = TeX("$\\beta_2$"),
  alpha1 = TeX("$\\alpha_1$"),
  alpha2 = TeX("$\\alpha_2$"),
  eta1 = TeX("$\\eta_1$"),
  eta2 = TeX("$\\eta_2$")
)
ggplot(plot_data, aes(x = param, y = est)) +
  geom_boxplot(alpha = 0.5, fill = btfgreen) +
  scale_x_discrete(labels = param_labels) +
  btf_theme +
  labs(x = "Parameter (km/h)", y = "Estimated Values",
       title = "")


foldername <- paste(im_folder, "optim/omsev/estimations_robustness/boxplot_estimations_free_eta_q",
                  q * 100, "_delta", delta, "_dmin", min_spatial_dist,
                  ".png", sep = "")
ggsave(foldername, width = 20, height = 15, units = "cm")


init_params_com <- c(params_com$beta1, params_com$beta2,
                     params_com$alpha1, params_com$alpha2)

# change_init <- c(0.2, 1, 1, 1, 8, 2)
result <- optim(par = init_params_com[1:4], fn = neg_ll_composite_fixed_eta,
        list_lags = list_lags, list_episodes = list_episodes,
        list_excesses = list_excesses, hmax = hmax,
        wind_df = V_episodes,
        latlon = FALSE,
        fixed_eta1 = params_com$eta1,
        fixed_eta2 = params_com$eta2,
        distance = "lalpha",
        method = "L-BFGS-B",
        lower = c(1e-08, 1e-08, 1e-08, 1e-08),
        upper = c(10, 10, 1.999, 1.999),
        control = list(maxit = 10000,
                      trace = 1),
        hessian = F)

result


init_list_fixed <- list(
  init_params_com,
  c(0.5, 5, 0.2, 0.7),
  c(1, 1, 0.5, 0.5),
  c(0.1, 0.1, 0.1, 0.1),
  c(2, 8, 0.3, 1.2)
)

results_fixed <- lapply(init_list_fixed, function(init_params_com) {
  optim(par = init_params_com[1:4], fn = neg_ll_composite_fixed_eta,
        list_lags = list_lags, list_episodes = list_episodes,
        list_excesses = list_excesses, hmax = hmax,
        wind_df = V_episodes,
        latlon = FALSE,
        fixed_eta1 = params_com$eta1,
        fixed_eta2 = params_com$eta2,
        distance = "lalpha",
        method = "L-BFGS-B",
        lower = c(1e-08, 1e-08, 1e-08, 1e-08),
        upper = c(10, 10, 1.999, 1.999),
        control = list(maxit = 10000, trace = 0),
        hessian = FALSE)
})

for (i in seq_along(results_fixed)) {
  cat("\n--- Fixed-eta Init set", i, "---\n")
  print(results_fixed[[i]]$par)
  cat("Neg log-lik:", results_fixed[[i]]$value, "\n")
}

# plot results_fixed
results_table_fixed <- data.frame(
  init_beta1 = sapply(init_list_fixed, function(x) x[1]),
  init_beta2 = sapply(init_list_fixed, function(x) x[2]),
  init_alpha1 = sapply(init_list_fixed, function(x) x[3]),
  init_alpha2 = sapply(init_list_fixed, function(x) x[4]),
  est_beta1 = sapply(results_fixed, function(res) res$par[1]),
  est_beta2 = sapply(results_fixed, function(res) res$par[2]),
  est_alpha1 = sapply(results_fixed, function(res) res$par[3]),
  est_alpha2 = sapply(results_fixed, function(res) res$par[4]),
  neg_loglik = sapply(results_fixed, function(res) res$value)
)
print(results_table_fixed)
params <- c("beta1","beta2","alpha1","alpha2")
plot_data_fixed <- data.frame(init = unlist(results_table_fixed[, paste0("init_", params)]),
                        est  = unlist(results_table_fixed[, paste0("est_", params)]),
                        param = rep(params, each = nrow(results_table_fixed)))

# Box plot
param_labels <- c(
  beta1 = TeX("$\\beta_1$"),
  beta2 = TeX("$\\beta_2$"),
  alpha1 = TeX("$\\alpha_1$"),
  alpha2 = TeX("$\\alpha_2$")
)
ggplot(plot_data_fixed, aes(x = param, y = est)) +
  geom_boxplot(alpha = 0.5) +
  geom_jitter(width = 0.1) +
  scale_x_discrete(labels = param_labels) +
  btf_theme +
  labs(x = "Parameters (km/h)", y = "Estimated Values",
       title = "")

foldername <- paste(im_folder, "optim/omsev/estimations_robustness/boxplot_estimations_fixed_eta_q",
                  q * 100, "_delta", delta, "_dmin", min_spatial_dist,
                  ".png", sep = "")
ggsave(foldername, width = 20, height = 15, units = "cm")


# Conversion factors
c_x_km <- 1      # for km/h
c_t_h <- 1
c_x_m <- 1000    # for m/5min
c_t_5min <- 12   # 1 hour = 12 * 5min
c_v_m5min <- 1000 / 12  # from km/h to m/5min
# 0.5677 & 4.7658 & 0.1519 & 0.7126 & 2.0907 & 3.7794 & 19174.24 \\
beta1_hat <- result$par[1]
beta2_hat <- result$par[2]
alpha1_hat <- result$par[3]
alpha2_hat <- result$par[4]
eta1_hat <- result$par[5]
eta2_hat <- result$par[6]
# Parameters in m/5min
params_m5min <- convert_params(beta1_hat, beta2_hat, alpha1_hat, alpha2_hat,
                               eta1 = eta1_hat, eta2 = eta2_hat,
                               c_x = c_x_m, c_t = c_t_5min)
final_res <- c(params_m5min$beta1, params_m5min$beta2,
                alpha1_hat, alpha2_hat)
print(final_res)
print(result$value)
################################################################################
# OTHER MODEL with omega
################################################################################

omega_com <- 1
# With omega
beta1_com <- init_params_com[1]
beta2_com <- init_params_com[2]
alpha1_com <- init_params_com[3]
alpha2_com <- init_params_com[4]
eta1_com <- init_params_com[5]
eta2_com <- init_params_com[6]
init_param <- c(beta1_com, beta2_com, alpha1_com, alpha2_com, omega_com, eta2_com)

hmax <- max(dist_mat) / 1000 # convert to km
df_lags <- list_lags[[1]]
head(df_lags)
# change_init <- c(0.2, 1, 1, 1, 8, 2)
result <- optim(par = init_param, fn = neg_ll_composite_omega,
        list_lags = list_lags, list_episodes = list_episodes,
        list_excesses = list_excesses, hmax = hmax,
        V_episode = V_episodes, W_episode = W_episodes,
        latlon = TRUE,
        fixed_omega = NA,
        fixed_eta = 1,
        distance = "lalpha",
        method = "L-BFGS-B",
        lower = c(1e-08, 1e-08, 1e-08, 1e-08, 1e-08, 1e-08),
        upper = c(10, 10, 1.999, 1.999, 1, 10),
        control = list(maxit = 10000,
                      trace = 1),
        hessian = F)

result$par

# nothing fixed
# [1] 0.5461316 4.7832909 0.0970711 0.7212678 0.9004272 3.3980109

# eta2= 2
# [1] 0.5520984 4.7973882 0.1083775 0.7234912 0.7719444 2.0600000

init_param <- c(beta1_com, beta2_com, alpha1_com, alpha2_com, omega_com, eta1_com, eta2_com)

result <- optim(par = init_param, fn = neg_ll_composite_new,
        list_lags = list_lags, list_episodes = list_episodes,
        list_excesses = list_excesses, hmax = hmax,
        V_episode = V_episodes, W_episode = W_episodes,
        latlon = TRUE,
        fixed_omega = NA,
        fixed_eta1 = init_param[6],
        fixed_eta2 = init_param[7],
        distance = "lalpha",
        method = "L-BFGS-B",
        lower = c(1e-08, 1e-08, 1e-08, 1e-08, 1e-08, 1e-08, 1e-08),
        upper = c(10, 10, 1.999, 1.999, 1, 10, 10),
        control = list(maxit = 10000,
                      trace = 1),
        hessian = F)

# count number of zero adv and zero wind
sum(V_episodes$vx == 0 & V_episodes$vy == 0)
sum(W_episodes$wx == 0 & W_episodes$wy == 0)

V_episodes_etas <- V_episodes
V_episodes_etas$vx <- eta1_com * sign(V_episodes$vx) * abs(V_episodes$vx)^eta2_com
V_episodes_etas$vy <- eta1_com * sign(V_episodes$vy) * abs(V_episodes$vy)^eta2_com

# put NA for outliers in V_episodes_etas
V_episodes_etas$vx[abs(V_episodes_etas$vx) > 50] <- NA
V_episodes_etas$vy[abs(V_episodes_etas$vy) > 50] <- NA

na_indices <- is.na(V_episodes_etas$vx) | is.na(V_episodes_etas$vy)
V_episodes_etas <- V_episodes_etas[!na_indices, ]
list_lags_etas <- list_lags[!na_indices]
list_excesses_etas <- list_excesses[!na_indices]
list_episodes_etas <- list_episodes[!na_indices]


result <- optim(par = init_params_com[1:4], fn = neg_ll_composite_fixed_eta,
        list_lags = list_lags, list_episodes = list_episodes,
        list_excesses = list_excesses, hmax = 7,
        wind_df = V_episodes,
        latlon = TRUE,
        fixed_eta1 = init_params_com[5],
        fixed_eta2 = init_params_com[6],
        distance = "lalpha",
        method = "L-BFGS-B",
        lower = c(1e-08, 1e-08, 1e-08, 1e-08, 1e-08, 1e-08),
        upper = c(10, 10, 1.999, 1.999,  10, 10),
        control = list(maxit = 10000,
                      trace = 1),
        hessian = F)

result


################################################################################

# Fictive advection vectors for representation
# Fictive advection vectors for representation
# Create a grid of directions (angles) and speeds > 1 km/h
# Focus on northward directions (angles between -90° and +90° from north)
angles_deg <- seq(-90, 90, by = 30)  # angles from north (-90° = west, 90° = east)
speeds_kmh <- c(1, 2, 5, 7)  # speeds > 1 km/h

fictive_adv <- expand.grid(angle = angles_deg, speed = speeds_kmh)

# Convert to Cartesian coordinates
# North is positive y-axis, East is positive x-axis
fictive_adv$adv_x <- fictive_adv$speed * sin(fictive_adv$angle * pi / 180)
fictive_adv$adv_y <- fictive_adv$speed * cos(fictive_adv$angle * pi / 180)

# Keep only adv_x and adv_y columns
fictive_adv <- fictive_adv[, c("adv_x", "adv_y")]


# group strong.N = 47 episodes
# $par
#     beta1     beta2    alpha1    alpha2 
# 1.8281024 4.1616951 0.4804306 0.6912886 

# $value
# [1] 1839.098
# $beta1
#      beta1 
# 0.06617724 

# $beta2
#     beta2 
# 0.7468692 
beta1_hat <- 1.83
beta2_hat <- 4.16
alpha1_hat <- 0.48
alpha2_hat <- 0.69
eta1_hat <- result$par[5]
eta2_hat <- result$par[6]
params_m5min <- convert_params(beta1_hat, beta2_hat, alpha1_hat, alpha2_hat,
                               eta1 = eta1_hat, eta2 = eta2_hat,
                               c_x = c_x_m, c_t = c_t_5min)

# Estimated parameters from optimization
theta_hat_kmh <- result$par[1:4]  # beta1, beta2, alpha1, alpha2
theta_hat_m5min <- c(params_m5min$beta1, params_m5min$beta2,
                     alpha1_hat, alpha2_hat)
theta_hat <- theta_hat_kmh
eta1_hat <- 0.6665347
eta2_hat <- 1.75743


# Lag distances (km)
h_vals <- seq(0, 15, by = 0.05)

# Time lags (converted to hours)
tau_vals <- c(0, 1, 3, 5, 7, 10) * 5 / 60

# Standard spatial directions (unit vectors)
directions_named <- list(
  EW = c(1, 0),
  NS = c(0, 1),
  Diagonal = c(1, 1) / sqrt(2)
)

# Output folder for saving plots (update path as needed)
foldername <- paste0(im_folder, "optim/omsev/variogram/q", q * 100,
                     "_delta", delta, "_dmin", min_spatial_dist, "_fictive_adv/")
if (!dir.exists(foldername)) dir.create(foldername, recursive = TRUE)
foldername_hnorm <- paste0(foldername, "hnorm/")
if (!dir.exists(foldername_hnorm)) dir.create(foldername_hnorm, recursive = TRUE)
foldername_hnormV <- paste0(foldername, "hnormV/")
if (!dir.exists(foldername_hnormV)) dir.create(foldername_hnormV, recursive = TRUE)


# Loop over all fictive advection vectors
for (i in seq_len(nrow(fictive_adv))) {
  fictive_v <- c(fictive_adv$adv_x[i], fictive_adv$adv_y[i])
  
  # Normalize advection direction vector
  norm_v <- sqrt(sum(fictive_v^2))
  if (norm_v > 0) {
    dir_adv <- fictive_v / norm_v
  } else {
    dir_adv <- c(1, 0)  # default if zero vector
  }
  
  # Add the advection direction to standard directions
  directions_all <- c(directions_named, list(Advection = dir_adv))
  
  # Loop over all directions
  # for (dname in names(directions_all)) {
  
    direction <- NULL
    
    # Compute gamma grid with corrected distance for each direction and advection vector
    df_gamma <- compute_gamma_grid(h_vals, tau_vals, direction,
                                   theta_hat,
                                   eta1 = eta1_hat, eta2 = eta2_hat,
                                   fictive_v = fictive_v)
    
    # Plot gamma as function of corrected distance (dist_corr)
    subtitle_txt <- paste0("Advection: V = (",
                           round(fictive_v[1], 3), ", ", round(fictive_v[2], 3), ") km/h")
    
    # p <- ggplot(df_gamma, aes(x = h_normV, y = gamma, color = factor(tau_min))) +
    #     geom_line(size = 1.2) +
    #     labs(
    #       subtitle = subtitle_txt,
    #       x = expression("||h|| (km)"),
    #       y = expression(gamma(h, tau)),
    #       color = expression(tau ~ "(min)")
    #     ) +
    #     theme_minimal()


    
    # filename <- paste0(foldername_hnormV, "variogram_fictive", i, "_dir_", dname, ".pdf")
    # ggsave(filename, plot = p, width = 20, height = 15, units = "cm", dpi = 600)


    p <- ggplot(df_gamma, aes(x = h_norm, y = gamma, color = factor(tau_min))) +
        geom_line(size = 1.2) +
        labs(
          subtitle = subtitle_txt,
          x = expression("||h|| (km)"),
          y = expression(gamma(h, tau)),
          color = expression(tau ~ "(min)")
        ) +
        theme_minimal()


    
    filename <- paste0(foldername_hnorm, "variogram_fictive", i, ".pdf")
    ggsave(filename, plot = p, width = 20, height = 15, units = "cm", dpi = 600)
  # }
}



# plot without direction
df_gamma <- compute_gamma_grid_no_direction(h_vals, tau_vals,
                                       theta_hat,
                                       eta1 = eta1_hat, eta2 = eta2_hat,
                                       fictive_v = c(3, 4))