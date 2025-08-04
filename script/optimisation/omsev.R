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
filename_omsev <- paste0(data_folder,
                         "omsev/omsev_5min/rain_mtp_5min_2019_2024.RData")

load(filename_omsev)

# get location of each rain gauge
filename_loc <- paste0(data_folder, "omsev/loc_rain_gauges.csv")
location_gauges <- read.csv(filename_loc)
location_gauges$Station <- c("iem", "mse", "poly", "um", "cefe", "cnrs",
                             "crbm", "archiw", "archie", "um35", "chu1",
                             "chu2", "chu3", "chu4", "chu5", "chu6", "chu7",
                             "cines", "brives", "hydro")



# rain.all5 is the data frame name for 5 min data
# get only stations records and dates
rain <- as.data.frame(rain.all5[, c(1, 6:(ncol(rain.all5) - 1))])

# remove rain before september 2019 and after january 2024
# rain$dates <- as.POSIXct(rain.all5$dates, tz = "Europe/Paris")
rain$dates <- with_tz(rain$dates, tzone = "UTC")
rain <- rain[rain$dates >= "2019-01-01" & rain$dates <= "2024-02-01", ]
rownames(rain) <- rain$dates
head(rain)
tail(rain)
# Remove non real zeros
filename_omsev1min <- paste0(data_folder,
                         "omsev/omsev_1min/rain_mtp_1min_2019_2024.RData")
load(filename_omsev1min)
rain1min <- Rain.all.treated
colnames(rain1min) <- c("dates", "archie", "archiw", "cefe", "chu1", "chu2",
                        "chu3", "chu4", "chu5", "chu6", "chu7", "cines", "cnrs",
                        "crbm", "hydro", "iem", "mse", "poly", "brives", "um",
                        "um35")

results <- data.frame(
  Station = character(),
  First_Date = as.Date(character()),
  Last_Date = as.Date(character()),
  stringsAsFactors = FALSE
)

rain1min$dates <- as.POSIXct(rain1min$dates, tz = "Europe/Paris")
rain1min$dates <- with_tz(rain1min$dates, tzone = "UTC")
# Loop through all columns except "Dates"
for (col in colnames(rain1min)[colnames(rain1min) != "dates"]) {
  values <- rain1min[[col]]
  non_na_indices <- which(!is.na(values))
  
  if (length(non_na_indices) > 0) {
    first_date <- rain1min$dates[non_na_indices[1]]
    last_date <- rain1min$dates[non_na_indices[length(non_na_indices)]]
    
    results <- rbind(results, data.frame(
      Station = col,
      First_Date = first_date,
      Last_Date = last_date
    ))
  }
}

results
library(dplyr)
library(tidyr)

# Pivot rain1min to long format
rain_long <- rain %>%
  pivot_longer(-dates, names_to = "Station", values_to = "Value")

# Join with results to get First_Date and Last_Date
rain_long <- rain_long %>%
  left_join(results, by = "Station")

# Replace false zeros (before first or after last date) with NA
rain_long <- rain_long %>%
  mutate(Value = if_else(
    Value == 0 & (dates < First_Date | dates > Last_Date),
    NA_real_,
    Value
  )) %>%
  select(dates, Station, Value)  # Optional: drop First/Last date columns

# Pivot back to wide format
rain_clean <- rain_long %>%
  pivot_wider(names_from = Station, values_from = Value)

head(rain_clean)

rain <- as.data.frame(rain_clean)
# put dates as rownames
rownames(rain) <- rain$dates
rain <- rain[-1] # remove dates column
# rain$mse
# Get distances matrix
dist_mat <- get_dist_mat(location_gauges)
df_dist <- reshape_distances(dist_mat)
# x11()
plot(df_dist$value)
max(df_dist$value)

library(dplyr)

# remove all rows with all NAs
rain <- rain %>%
  filter(!if_all(everything(), is.na))

# in rain remove when all data are NA
rain <- rain[rowSums(is.na(rain)) < ncol(rain), ]

# remove cines, hydro, brives
colnames(rain)
rain <- rain[, !(colnames(rain) %in% c("cines", "hydro", "brives"))]

location_gauges <- location_gauges[location_gauges$Station != "cines" &
                                   location_gauges$Station != "hydro" &
                                   location_gauges$Station != "brives", ]
dist_mat <- get_dist_mat(location_gauges)
df_dist <- reshape_distances(dist_mat)
plot(df_dist$value)
max(df_dist$value)

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
foldername <- paste0(data_folder, "omsev/WLSE/")
df_spa_result <- read.csv(paste0(foldername, "wlse_spa_summary.csv"))
df_temp_result <- read.csv(paste0(foldername, "wlse_temp_summary.csv"))
df_result_all <- merge(df_spa_result, df_temp_result, by = NULL)

# select one row
df_result <- df_result_all[df_result_all$q_spa == 0.95 &
                             df_result_all$q_temp == 0.9, ]

beta1 <- df_result$beta1
beta2 <- df_result$beta2
alpha1 <- df_result$alpha1
alpha2 <- df_result$alpha2

################################################################################
# VARIOGRAM --------------------------------------------------------------------
################################################################################

# estimates
param_omsev <- c(beta1, beta2, alpha1, alpha2) # m / min

# in rain remove when all data are NA<
rain <- rain[rowSums(is.na(rain)) < ncol(rain), ]
q <- 0.97 # quantile
set_st_excess <- get_spatiotemp_excess(rain, quantile = q, remove_zeros = TRUE)

# verify that the excess is above the threshold
# get list of sites and times
list_s <- set_st_excess$list_s
unique(unlist(list_s))
list_t <- set_st_excess$list_t

rain[list_t[[1]], list_s[[1]]]

list_u <- set_st_excess$list_u
list_u[[1]]
unique(unlist(list_u))


# Check all selected points to see if they exceed the threshold
excess_check <- sapply(seq_along(list_s), function(i) {
  value <- rain[list_t[[i]], list_s[[i]]]
  threshold <- list_u[[i]]
  value > threshold
})

# Print summary
all(excess_check)  # Should be TRUE if all exceedances are valid


# Spatio-temporal neighborhood parameters
min_spatial_dist <- 500 # m
delta <- 12 # in * 5 min
episode_size <- delta # size of the episode
sites_coords <- location_gauges[, c("Longitude", "Latitude")]
tail(rain)
s0t0_set <- get_s0t0_pairs(sites_coords, rain,
                            min_spatial_dist = min_spatial_dist,
                            episode_size = episode_size,
                            set_st_excess = set_st_excess,
                            n_max_episodes = 10000,
                            latlon = TRUE)

selected_points <- s0t0_set
selected_points[12,]

length(selected_points$s0) # number of selected points


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
print(n_episodes)
length(unique(selected_points$s0)) # can be same s0
length(unique(selected_points$t0)) # never same t0?
print(min(selected_points$u_s0)) # min threshold
t0_list <- selected_points$t0
s0_list <- selected_points$s0

list_episodes_points <- get_extreme_episodes(selected_points, rain,
                                     episode_size = episode_size, unif = FALSE,
                                     beta = 0)

list_episodes <- list_episodes_points$episodes
episode <- list_episodes[[1]]

selected_points$t0_date <- as.POSIXct(selected_points$t0_date, format = "%Y-%m-%d %H:%M:%S", tz = "UTC")

# selected_points$t0_date <- ifelse(
#   nchar(format(selected_points$t0_date, "%Y-%m-%d %H:%M:%S")) == 10,
#   paste0(format(selected_points$t0_date, "%Y-%m-%d"), " 00:00:00"),
#   format(selected_points$t0_date, "%Y-%m-%d %H:%M:%S")
# )

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


library(ggplot2)
library(dplyr)

# Create output folder if it doesn't exist
rain_plot_folder <- file.path(im_folder, "optim/omsev/rain_episodes/")
if (!dir.exists(rain_plot_folder)) dir.create(rain_plot_folder, recursive = TRUE)

# Loop over all datetimes and plot
for (i in seq_along(datetimes)) {
  # Define window (e.g., -2h to +10h around dt)
  dt <- as.POSIXct(datetimes[i], tz = "UTC")
  start_date <- dt - 2 * 3600
  end_date <- dt + 10 * 3600

  rain_subset <- rain[rownames(rain) >= start_date & rownames(rain) <= end_date, , drop = FALSE]
  rain_subset$datetime <- as.POSIXct(rownames(rain_subset), tz = "UTC")

  rain_subset_long <- rain_subset %>%
    pivot_longer(cols = -datetime, names_to = "Station", values_to = "Value") %>%
    filter(!is.na(Value))

  plot_title <- paste0("Rainfall around ", format(dt, "%Y-%m-%d %H:%M"))
  
  rain_plot <- ggplot(rain_subset_long, aes(x = datetime, y = Value, color = Station)) +
    geom_line(size = 0.5) +
    geom_vline(xintercept = as.numeric(dt), linetype = "dashed", color = "black") +  # dashed vertical line at dt
    labs(title = plot_title, x = "Time", y = "Rainfall (mm)") +
    theme_minimal() +
    theme(legend.position = "right")

  # Save plot
  plot_filename <- file.path(rain_plot_folder, paste0("rain_", format(dt, "%Y%m%d_%H%M"), ".png"))
  ggsave(plot_filename, plot = rain_plot, width = 12, height = 6, units = "in", dpi = 150)
}


# get comephore advection for each episode
adv_filename <- paste(data_folder, "/omsev/adv_estim/advection_results_q",
                      q * 100, "_delta", delta, "_dmin", min_spatial_dist,
                      ".csv", sep = "")
adv_df <- read.csv(adv_filename, sep = ",")
head(adv_df)

# convert adv in m/s to km/h
adv_df$mean_dx_kmph <- adv_df$mean_dx_mps * 3.6
adv_df$mean_dy_kmph <- adv_df$mean_dy_mps * 3.6
# Transform t0 to POSIXct
adv_df$t0 <- as.POSIXct(adv_df$t0, format = "%Y-%m-%d %H:%M:%S", tz = "UTC")

# get only matching episodes from selected_points
matching_indices <- match(selected_points$t0_date_rounded, adv_df$t0)
# remove NA indices
matching_indices <- matching_indices[!is.na(matching_indices)]
# get only matching rows
adv_df <- adv_df[matching_indices, ]
rownames(adv_df) <- NULL  # reset row names
head(adv_df)

nrow(adv_df) # number of episodes with advection estimates

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
    adv_x <- adv_row$mean_dx_kmph
    adv_y <- adv_row$mean_dy_kmph
    selected_episodes$adv_x[i] <- adv_x
    selected_episodes$adv_y[i] <- adv_y
  }
}

head(selected_episodes)
tail(selected_episodes)
# Remove episodes with NA advection values
# ind_NA_adv <- which(is.na(selected_episodes$adv_x) | is.na(selected_episodes$adv_y))
# ind_0_adv <- which(selected_episodes$adv_x == 0 & selected_episodes$adv_y == 0)
# selected_episodes_nona <- selected_episodes[-ind_0_adv, ]
selected_episodes_nona <- selected_episodes
wind_df <- selected_episodes_nona[, c("adv_x", "adv_y")]
colnames(wind_df) <- c("vx", "vy")
length(wind_df$vx) # should be the same as number of episodes

tau_vect <- 0:10
thresholds_by_site <- apply(rain, 2, function(col) {
  col <- col[!is.na(col) & col > 0]
  if (length(col) < 30) return(NA)
  quantile(col, probs = q)
})

# plot wind df vectors
library(ggplot2)
library(grid)  # pour arrow()

wind_df_plot <- ggplot(selected_episodes_nona, aes(x = 0, y = 0)) +
  geom_segment(aes(xend = adv_x, yend = adv_y), 
               arrow = arrow(length = unit(0.2, "cm")),
               color = btfgreen, size = 0.5) +
  geom_point(aes(x = adv_x, y = adv_y), color = btfgreen, size = 1) +
  btf_theme +
  coord_equal() +
  xlab("Advection in x (km/h)") +
  ylab("Advection in y (km/h)") +
  ggtitle("Advection vectors for selected episodes") +
  theme(plot.title = element_text(hjust = 0.5)) 

wind_df_plot

# save wind df plot
filename <- paste(im_folder, "optim/omsev/advection_com_emp_plot_q",
                  q * 100, "_delta", delta, "_dmin", min_spatial_dist,
                  ".png", sep = "")
ggsave(filename, plot = wind_df_plot, width = 20, height = 15, units = "cm",
       dpi = 600, device = "png")

selected_episodes_nona$speed <- sqrt(selected_episodes_nona$adv_x^2 + 
                                       selected_episodes_nona$adv_y^2)

selected_episodes_nona$angle_deg <- atan2(selected_episodes_nona$adv_x, selected_episodes_nona$adv_y) * 180 / pi
# Mettre l’angle dans [0, 360[
selected_episodes_nona$angle_deg <- (selected_episodes_nona$angle_deg + 360) %% 360


get_cardinal <- function(angle) {
  directions <- c("N", "NE", "E", "SE", "S", "SW", "W", "NW")
  idx <- round(angle / 45) %% 8 + 1
  return(directions[idx])
}

selected_episodes_nona$direction <- sapply(selected_episodes_nona$angle_deg, get_cardinal)


# remove those with 0 advection
# selected_points_adv <- selected_points_nona[
  # selected_points_nona$adv_x != 0 | selected_points_nona$adv_y != 0, ]
s0_list <- selected_episodes_nona$s0
list_episodes_points <- get_extreme_episodes(selected_episodes_nona, rain,
                              episode_size = episode_size, unif = FALSE)

list_episodes <- list_episodes_points$episodes

tmax <- max(tau_vect)
df_coords <- as.data.frame(sites_coords)
# Compute the lags and excesses for each conditional point
list_results <- mclapply(1:length(s0_list), function(i) {
  s0 <- s0_list[i]
  col_s0 <- which(colnames(rain) == s0)
  s0_coords <- df_coords[col_s0, ]
  # t0 <- t0_list[i]
  episode <- list_episodes[[i]]
  ind_t0_ep <- 0 # index of t0 in the episode
  lags <- get_conditional_lag_vectors(df_coords, s0_coords, ind_t0_ep,
                                tau_vect, latlon = TRUE)
  # hnorm is in meters
  lags$hnorm <- lags$hnorm / 1000 # convert to km

  # excesses <- empirical_excesses_rpar(episode, q, lags, t0 = ind_t0_ep)

  excesses <- empirical_excesses_rpar(episode, thresholds = thresholds_by_site,
                                      lags, t0 = ind_t0_ep)
  # tau is in 5 minutes
  lags$tau <- lags$tau * 5 / 60 # convert to hours
  list(lags = lags, excesses = excesses)
}, mc.cores = detectCores() - 1)

list_lags <- lapply(list_results, `[[`, "lags")
list_excesses <- lapply(list_results, `[[`, "excesses")

s0 <- s0_list[1]
col_s0 <- which(colnames(rain) == s0)
s0_coords <- df_coords[col_s0, ]
excesses <- list_excesses[[20]]
sum(excesses$kij)
df_lags <- list_lags[[1]]
tail(df_lags)

# wind_df test
# adv_df <- data.frame(vx = 1, vy = 1)

eta1_com <- 4
eta2_com <- 2

init_param <- c(beta1, beta2, alpha1, alpha2, eta1_com, eta2_com)
init_param <- c(beta1, beta2, alpha1, alpha2, 1, 1)

# init_param <- c(0.02, 0.5, alpha1, alpha2, 1, 1)

hmax <- max(dist_mat) / 1000 # convert to km
result <- optim(par = init_param, fn = neg_ll_composite,
        list_lags = list_lags, list_episodes = list_episodes,
        list_excesses = list_excesses, hmax = hmax,
        wind_df = wind_df,
        latlon = TRUE,
        directional = TRUE,
        fixed_eta1 = TRUE,
        fixed_eta2 = TRUE,
        method = "L-BFGS-B",
        lower = c(1e-08, 1e-08, 1e-08, 1e-08, 1e-08, 1e-08),
        upper = c(10, 10, 1.999, 1.999, 10, 10),
        control = list(maxit = 10000,
                      trace = 1,
                      parscale = c(1, 1, 1, 1, 1, 1)),
        hessian = T)

result

beta1_hat <- result$par[1]
beta2_hat <- result$par[2]
alpha1_hat <- result$par[3]
alpha2_hat <- result$par[4]
eta1_hat <- result$par[5]
eta2_hat <- result$par[6]
hessian <- result$hessian
# remove eta1 and eta2 from the hessian
hessian_noeta <- hessian[-(5:6), -(5:6)]
vcov_matrix <- solve(hessian_noeta)


library(numDeriv)

compute_gamma_ic_grid <- function(h_vals, tau_vals, direction, wind_i,
                                  theta_hat, vcov_matrix, eta1 = 1, eta2 = 1) {
  df_out <- data.frame()
  
  for (tau in tau_vals) {
    for (h in h_vals) {
      hx <- h * direction[1]
      hy <- h * direction[2]

      gamma_func <- function(theta) {
        beta1 <- theta[1]; beta2 <- theta[2]
        alpha1 <- theta[3]; alpha2 <- theta[4]
        if (length(theta) > 4) {
          eta1 <- theta[5]; eta2 <- theta[6]
        } 
        vx <- eta1 * abs(wind_i$vx)^eta2 * sign(wind_i$vx)
        vy <- eta1 * abs(wind_i$vy)^eta2 * sign(wind_i$vy)

        term1 <- beta1 * abs(hx - vx * tau)^alpha1
        term2 <- beta2 * abs(hy - vy * tau)^alpha1
        term3 <- beta2 * tau^alpha2
        return(term1 + term2 + term3)
      }

      grad_gamma <- grad(gamma_func, theta_hat)
      var_gamma <- t(grad_gamma) %*% vcov_matrix %*% grad_gamma
      se_gamma <- sqrt(var_gamma)
      gamma_hat <- gamma_func(theta_hat)

      df_out <- rbind(df_out, data.frame(
        h = h,
        tau = tau,
        gamma = gamma_hat,
        gamma_low = gamma_hat - 1.96 * se_gamma,
        gamma_high = gamma_hat + 1.96 * se_gamma
      ))
    }
  }

  return(df_out)
}

# Paramètres estimés
theta_hat <- result$par[1:4]  # beta1, beta2, alpha1, alpha2
library(ggplot2)
library(numDeriv)

# Named standard directions
directions_named <- list(
  EW = c(1, 0),
  NS = c(0, 1),
  Diagonal = c(1, 1) / sqrt(2)
)

# Selected episodes
episode_ids <- c(1, 3, 5, 10, 20, 100)

# Spatial and temporal lags
dmax <- max(dist_mat) / 1000  # convert m to km
dmax<- 1
h_vals <- seq(0, dmax, by = 0.05)
tau_vals <- c(0, 1, 3, 5, 10) * 5 / 60  # in hours (5-min intervals)



str_zero_wind <- ""
if (!any(wind_df$vx == 0 & wind_df$vy == 0)) {
  str_zero_wind <- "_no_0adv"
} 

# Output folder
foldername <- paste0(im_folder, "optim/omsev/variogram/q", q * 100,
                     "_delta", delta, "_dmin", min_spatial_dist, 
                     str_zero_wind, "/")

if (!dir.exists(foldername)) dir.create(foldername, recursive = TRUE)

# Main loop
for (i in episode_ids) {
  print(i)
  wind_i <- wind_df[i, ]
  vx_i <- wind_i$vx
  vy_i <- wind_i$vy
  V_vec <- c(vx_i, vy_i)

  norm_V <- sqrt(sum(V_vec^2))
  
  # Skip advection if wind is zero
  if (norm_V == 0) {
    # message("Episode ", i, ": zero wind → skipping advection direction.")
    directions_named_ext <- directions_named
  } else {
    direction_adv <- V_vec / norm_V
    directions_named_ext <- c(directions_named, list(Advection = direction_adv))
  }

  # Loop over all directions
  for (dname in names(directions_named_ext)) {
    direction <- directions_named_ext[[dname]]

    df_gamma <- compute_gamma_ic_grid(h_vals, tau_vals, direction,
                                      wind_i, theta_hat, vcov_matrix,
                                      eta1 = eta1_hat, eta2 = eta2_hat)
    df_gamma$tau_min <- df_gamma$tau * 60  # for plotting in minutes

    subtitle_txt <- paste0("Direction: ", dname,
                           ", Episode ", i, ", V = (",
                           round(vx_i * 1000 / 60, 4), ", ",
                           round(vy_i * 1000 / 60, 4), ") m/s")

    p <- ggplot(df_gamma, aes(x = h, y = gamma, color = factor(tau_min))) +
      geom_line(size = 1.2) +
      geom_ribbon(aes(ymin = gamma_low, ymax = gamma_high, fill = factor(tau_min)),
                  alpha = 0.2, color = NA) +
      labs(
        title = "",
        subtitle = subtitle_txt,
        x = "Distance (km)",
        y = expression(gamma(h, tau)),
        color = expression(tau ~ "(min)"),
        fill = expression(tau ~ "(min)")
      ) +
      theme_minimal()

    filename <- paste0(foldername, "variogram_direction_",
                       dname, "_episode_", i, ".pdf")
    ggsave(filename, plot = p, width = 20, height = 15,
           units = "cm", dpi = 600)
  }
}


################################################################################

# Fictive advection vectors for representation
fictive_adv <- data.frame(
  adv_x = c(0, 0.01, 0.1, -1, -1, 1, 10, -2, 1),
  adv_y = c(0, 0.02, 0.5, -5, 5, 1, 11, -1, -1)
)

# Estimated parameters from optimization
theta_hat <- result$par[1:4]  # beta1, beta2, alpha1, alpha2
eta1_hat <- result$par[5]
eta2_hat <- result$par[6]

# Lag distances (km)
h_vals <- seq(0, 3, by = 0.05)

# Time lags (converted to hours)
tau_vals <- c(0, 1, 3, 5, 7, 10) * 5 / 60

# Standard spatial directions (unit vectors)
directions_named <- list(
  EW = c(1, 0),
  NS = c(0, 1),
  Diagonal = c(1, 1) / sqrt(2)
)

# Function to compute variogram gamma over grid of h and tau
# Arguments:
# - h_vals: spatial distances
# - tau_vals: temporal lags
# - direction: spatial direction vector (unit vector)
# - theta_hat: parameter vector (beta1, beta2, alpha1, alpha2)
# - eta1, eta2: parameters controlling advection scaling
# - fictive_v: advection velocity vector
compute_gamma_grid <- function(h_vals, tau_vals, direction,
                               theta_hat, eta1 = 1, eta2 = 1, fictive_v = c(1, 1)) {
  df_out <- data.frame()
  
  # Calculate scaled advection velocity components
  vx <- eta1 * abs(fictive_v[1])^eta2 * sign(fictive_v[1])
  vy <- eta1 * abs(fictive_v[2])^eta2 * sign(fictive_v[2])
  
  for (tau in tau_vals) {
    for (h in h_vals) {
      # Compute the spatial lag vector along given direction
      hx <- h * direction[1]
      hy <- h * direction[2]
      
      # Calculate the corrected spatial lag considering advection
      # Distance between h vector and advected position tau*V
      # dist_corr <- sqrt((hx - vx * tau)^2 + (hy - vy * tau)^2)

      vx <- eta1 * abs(fictive_v[1])^eta2 * sign(fictive_v[1])
      vy <- eta1 * abs(fictive_v[2])^eta2 * sign(fictive_v[2])

      # Projection sur la direction
      h_proj <- hx * direction[1] + hy * direction[2]
      v_proj <- vx * direction[1] + vy * direction[2]

      dist_corr <- abs(h_proj - v_proj * tau)

      beta1 <- theta_hat[1]
      beta2 <- theta_hat[2]
      alpha1 <- theta_hat[3]
      alpha2 <- theta_hat[4]
      
      # Variogram model: sum of spatial and temporal components
      term1 <- beta1 * dist_corr^alpha1
      term2 <- beta2 * tau^alpha2
      
      gamma_hat <- term1 + term2
      
      # Save results with corrected distance and parameters
      df_out <- rbind(df_out, data.frame(
        h = h,
        tau = tau,
        tau_min = tau * 60,  # for plotting in minutes
        gamma = gamma_hat,
        dist_corr = dist_corr
      ))
    }
  }
  return(df_out)
}

# Output folder for saving plots (update path as needed)
foldername <- paste0(im_folder, "optim/omsev/variogram/q", q * 100,
                     "_delta", delta, "_dmin", min_spatial_dist, "_fictive_adv/")
if (!dir.exists(foldername)) dir.create(foldername, recursive = TRUE)

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
  for (dname in names(directions_all)) {
    direction <- directions_all[[dname]]
    
    # Compute gamma grid with corrected distance for each direction and advection vector
    df_gamma <- compute_gamma_grid(h_vals, tau_vals, direction,
                                   theta_hat,
                                   eta1 = eta1_hat, eta2 = eta2_hat,
                                   fictive_v = fictive_v)
    
    # Plot gamma as function of corrected distance (dist_corr)
    subtitle_txt <- paste0("Advection: V = (",
                           round(fictive_v[1], 3), ", ", round(fictive_v[2], 3), "), direction: ", dname)
    
    p <- ggplot(df_gamma, aes(x = dist_corr, y = gamma, color = factor(tau_min))) +
      geom_line(size = 1.2) +
      labs(
        subtitle = subtitle_txt,
        x = expression(abs(h - tau * V)),
        y = expression(gamma(h, tau)),
        color = expression(tau ~ "(min)")
      ) +
      theme_minimal()
    
    filename <- paste0(foldername, "variogram_fictive", i, "_dir_", dname, ".pdf")
    ggsave(filename, plot = p, width = 20, height = 15, units = "cm", dpi = 600)
  }
}



