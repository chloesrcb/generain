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

library(sf)
library(dplyr)
library(ggplot2)
library(gridExtra)

load("workspace.RData")
# library(generain)

# LOAD DATA ####################################################################
filename_com <- paste0(data_folder, "comephore/comephore_2008_2024_5km.csv")
# filename_com <- paste0(data_folder, "comephore/zoom_5km.csv")
comephore_raw <- read.csv(filename_com, sep = ",")
filename_loc <- paste0(data_folder, "comephore/loc_px_zoom_5km.csv")
loc_px <- read.csv(filename_loc, sep = ",")
colnames(comephore_raw)

# remove pixel in loc_px that are not in comephore_raw
loc_px <- loc_px[loc_px$pixel_name %in% colnames(comephore_raw)[-1], ]
nrow(loc_px) # number of pixels
# reindex
rownames(loc_px) <- NULL
# length(colnames(comephore_raw)) - 1
# length(loc_px$pixel_name)

df_comephore <- as.data.frame(comephore_raw)

# when date is just a date, put 00:00:00 as time
if (nchar(df_comephore$date[1]) == 10) {
  df_comephore$date <- paste(df_comephore$date, "00:00:00")
}

df_comephore$date <- as.POSIXct(df_comephore$date,
                                format = "%Y-%m-%d %H:%M:%S", tz = "UTC")
head(df_comephore)

length(unique(df_comephore$date))
# Voir les duplicatas
duplicated_dates <- df_comephore[duplicated(df_comephore$date), ]
head(duplicated_dates)

# Take only data after 2007
# colnames(df_comephore)[1] <- "date"
df_comephore <- df_comephore[df_comephore$date >= "2008-01-01", ]
tail(df_comephore)
# put date in index
rownames(df_comephore) <- df_comephore$date
comephore <- df_comephore[-1] # remove dates column


# # get wind data
filename_wind <- paste0(data_folder, "wind/data_gouv/wind_mtp.csv")
wind_mtp <- read.csv(filename_wind)

# if datetime is just a date put 00:00:00 as time
if (nchar(wind_mtp$datetime[1]) == 10) {
  wind_mtp$datetime <- paste(wind_mtp$datetime, "00:00:00")
}


# Convert datetime to POSIXct
wind_mtp$datetime <- as.POSIXct(wind_mtp$datetime,
                                format = "%Y-%m-%d %H:%M:%S", tz = "UTC")

# Get only data after 2008
wind_mtp <- wind_mtp[wind_mtp$datetime >= "2008-01-01", ]

# Apply function to the DD column
wind_mtp$cardDir <- sapply(wind_mtp$DD, convert_to_cardinal, nb_cardinal = 8)
wind_mtp$cardDir <- as.character(wind_mtp$cardDir)  # Ensure it's character
wind_mtp$cardDir[is.na(wind_mtp$DD)] <- NA
# Check if NA values are properly handled
summary(wind_mtp)
# plot(wind_mtp$FF)
# 1 nd = 0,514 m/s
head(wind_mtp$cardDir)


# get year 2024
# wind_mtp <- wind_mtp[wind_mtp$datetime >= "2024-01-01", ]

# FF are in m/s 1/10 so we convert to m/s
# wind_mtp$FF <- wind_mtp$FF 
wind_mtp$FF_kmh <- wind_mtp$FF * 3.6  # convert to km/h

# plot(wind_mtp$FF_kmh, type = "l", col = "blue",
#      main = "Wind speed in km/h", xlab = "Time", ylab = "Wind speed (km/h)")

# DISTANCE AND COORDS ##########################################################

# Get distances matrix
# dist_mat <- get_dist_mat(loc_px)

# Get number of sites
nsites <- nrow(loc_px)

# Get coords
sites_coords <- loc_px[, c("Longitude", "Latitude")]
rownames(sites_coords) <- loc_px$pixel_name
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




p1 <- ggplot(grid_coords_km, aes(x = Longitude, y = Latitude)) +
  geom_point(color = "blue", size = 2) +
  geom_text(aes(label = rownames(grid_coords_km)), hjust = -0.2, size = 1.5) +
  coord_fixed() +
  ggtitle("GPS coordinates WGS84") +
  theme_minimal()

p2 <- ggplot(grid_coords_km, aes(x = x_km, y = y_km)) +
  geom_point(color = "red") +
  geom_text(aes(label = rownames(grid_coords_km)), size = 1.5, hjust = -0.2) +
  coord_fixed() +
  theme_minimal() +
  xlab("x in km") +
  ylab("y in km") +
  ggtitle("Transformed coordinates in km")

p_coords <- grid.arrange(p1, p2, ncol = 2)

# save plot
filename <- paste(im_folder, "optim/comephore/coords_transformation.png",
                  sep = "")
ggsave(plot = p_coords, filename = filename, width = 20, height = 15,
       units = "cm")

# get distance matrix
grid_coords_m <- grid_coords_m[, c("x_m", "y_m")]
grid_coords_km <- grid_coords_km[, c("x_km", "y_km")]
colnames(grid_coords_m) <- c("Longitude", "Latitude")
colnames(grid_coords_km) <- c("Longitude", "Latitude")
dist_mat <- get_dist_mat(grid_coords_m, latlon = FALSE)

# Spatial chi
df_dist <- reshape_distances(dist_mat)
# df_dist$value <- round(df_dist$value / 1000, 1) * 1000  # / 1000 in km
df_dist_km <- df_dist
df_dist_km$value <- round(df_dist$value / 1000, 1)


# Spatial chi WLSE #############################################################
h_vect <- sort(unique(df_dist_km$value))
h_vect <- h_vect[h_vect > 0]  # remove 0
hmax <- h_vect[10] # 10th value

hmax <- 7

q_no0_spa <- 0.97
chispa_df <- spatial_chi_alldist(df_dist_km, data_rain = comephore,
                            quantile = q_no0_spa, hmax = hmax, zeros = FALSE)

# chispa_df <- spatial_chi_extremogram(df_dist_km, comephore, q_no0_spa,
#                         hmax = hmax, zeros = FALSE)

etachispa_df <- data.frame(chi = eta(chispa_df$chi),
                           lagspa = log(chispa_df$lagspa))


chispa_plot <- ggplot(chispa_df, aes(lagspa, chi)) +
  btf_theme +
  geom_point(col = btfgreen) +
  xlab(TeX(r"($h$)")) +
  ylab(TeX(r"($\widehat{\chi}(h, 0)$)")) +
  ylim(0, 1)

chispa_plot

# save plot
filename <- paste(im_folder, "WLSE/comephore/full_spatial_chi_", q_no0_spa,
                  ".png", sep = "")
ggsave(filename, width = 20, height = 15, units = "cm")

# WLSE
wlse_spa <- get_estimate_variospa(chispa_df, weights = "exp", summary = TRUE)
print(wlse_spa)

c1 <- wlse_spa[[1]]
beta1 <- wlse_spa[[2]]
alpha1 <- wlse_spa[[3]]

chispa_eta_estim <- ggplot(etachispa_df, aes(lagspa, chi)) +
  btf_theme +
  geom_point(col = btfgreen) +
  xlab(TeX(r"($\log(h)$)")) +
  ylab(TeX(r"($\zeta(\widehat{\chi}(h, 0))$)")) +
  geom_line(aes(x = lagspa, y = alpha1 * lagspa + c1), alpha = 0.6,
            color = "darkred", linewidth = 0.5)

chispa_eta_estim

# save plot
filename <- paste(im_folder, "WLSE/comephore/full_spatial_chi_eta_estim_", q_no0_spa,
                  ".png", sep = "")
ggsave(filename, width = 20, height = 15, units = "cm")


# Temporal chi WLSE ############################################################

tmax <- 10
q_no0_temp <- 0.9 # quantile for temporal chi

# Temporal chi with spatial lag fixed at 0
# compute chiplot of every site with itself but lagged in time
chimat_dtlag <- temporal_chi(comephore, quantile = q_no0_temp, tmax = tmax,
                             mean = FALSE, zeros = FALSE)

chi_df_dt <- data.frame(chimat_dtlag)
colnames(chi_df_dt) <- c(0:tmax) # temporal lags from 1 to tmax
rownames(chi_df_dt) <- c(1:nsites) # stations
# remove lag 0
chi_df_dt <- chi_df_dt[, -1] # remove lag 0

# boxplot all pixels values for chi temp
# Reshape data using gather function
df_gathered <- chi_df_dt %>% gather(key = "variable", value = "value")
df_gathered$group <- factor(as.integer(df_gathered$variable),
                            levels = seq(0, tmax))

# Plot boxplots
chitemp <- ggplot(df_gathered, aes(x = group, y = value)) +
  geom_boxplot(fill = btfgreen, alpha = 0.4, outlier.color = btfgreen) +
  btf_boxplot_theme +
  xlab(TeX(r"($\tau$ (hours))")) +
  ylab(TeX(r"($\widehat{\chi}(0,\tau)$)"))

chitemp

# save plot
filename <- paste(im_folder, "WLSE/comephore/full_temporal_chi_boxplot_",
                q_no0_temp,  ".png", sep = "")
ggsave(filename, width = 20, height = 15, units = "cm")

# Mean of chi
chimat_dt_mean <- temporal_chi(comephore, tmax, quantile = q_no0_temp,
                               mean = TRUE, zeros = FALSE)
df_chi <- data.frame(lag = c(0:tmax), chi = chimat_dt_mean)

df_chi_not0 <- df_chi[df_chi$lag > 0, ]
wlse_temp <- get_estimate_variotemp(df_chi_not0, weights = "exp", summary = TRUE)
print(wlse_temp)
c2 <- as.numeric(wlse_temp[[1]])
beta2 <- as.numeric(wlse_temp[[2]])
alpha2 <- as.numeric(wlse_temp[[3]])
print(beta2)

dftemp <- data.frame(lag = log(df_chi_not0$lag), chi = eta(df_chi_not0$chi))

# remove first row
# dftemp <- dftemp[-1, ]

chitemp_eta_estim <- ggplot(dftemp, aes(x = lag, y = chi)) +
  geom_point(color = btfgreen) +
  btf_theme +
  xlab(TeX(r"($\log(\tau)$)")) +
  ylab(TeX(r"($\zeta(\widehat{\chi}(0,\tau))$)")) +
  geom_line(aes(x = lag, y = alpha2 * lag + c2),
            alpha = 0.5, color = "darkred", linewidth = 0.5)

chitemp_eta_estim

# save plot
filename <- paste(im_folder, "WLSE/comephore/full_temporal_chi_eta_estim_exp_",
                q_no0_temp, ".png", sep = "")
ggsave(filename, width = 20, height = 15, units = "cm")

# Result WLSE
df_result <- data.frame(beta1 =  beta1,
                        beta2 = beta2,
                        alpha1 = alpha1,
                        alpha2 = alpha2)

colnames(df_result) <- c("beta1", "beta2", "alpha1", "alpha2")

kable(df_result, format = "latex") %>%
  kable_styling(bootstrap_options = c("striped", "hover", "condensed",
  "responsive"), latex_options = "H")


# CHOOSE EXTREME EPISODE FOR R-PARETO ##########################################

q <- 0.97 # quantile

# get central site from sites_coords
# comephore_subset <- comephore[rownames(comephore) >= "2020-01-01", ]
comephore_subset <- comephore
set_st_excess <- get_spatiotemp_excess(comephore_subset, quantile = q,
                                      remove_zeros = TRUE)



# # get only year 2024
# comephore_subset <- comephore[rownames(comephore) >= "2008-01-01", ]
plot(comephore_subset[, 50], type = "l", col = "blue",
     main = "Comephore data for pixel 1 in 2024", xlab = "Time", ylab = "Rain (mm)")
list_s <- set_st_excess$list_s
list_t <- set_st_excess$list_t
list_u <- set_st_excess$list_u
# unique(list_u) # check unique excess values
# Spatio-temporal neighborhood parameters
min_spatial_dist <- 5 # in km
delta <- 30 # in hours
episode_size <- delta # size of the episode
s0t0_set <- get_s0t0_pairs(grid_coords_km, comephore_subset,
                            min_spatial_dist = min_spatial_dist,
                            episode_size = episode_size,
                            set_st_excess = set_st_excess,
                            n_max_episodes = 10000,
                            latlon = FALSE)

selected_points <- s0t0_set

# check that for all s0, t0 we have an excess above corresponding threshold
for (i in 1:length(selected_points$s0)) {
  s0 <- selected_points$s0[i]
  t0 <- selected_points$t0[i]
  u_s0 <- selected_points$u_s0[i]
  # check that the excess is above the threshold
  if (comephore_subset[t0, s0] <= u_s0) {
    stop(paste("Excess is not above threshold for s0 =", s0, "and t0 =", t0))
  }
}


library(dplyr)
library(lubridate)
library(readr)
time_lookup <- tibble(
  t0 = unlist(list_t),
  t0_date = parse_date_time(names(list_t), orders = c("ymd HMS", "ymd HM", "ymd"), tz = "UTC"),
  day = day(parse_date_time(names(list_t), orders = c("ymd HMS", "ymd HM", "ymd"), tz = "UTC")),
  month = month(parse_date_time(names(list_t), orders = c("ymd HMS", "ymd HM", "ymd"), tz = "UTC")),
  year = year(parse_date_time(names(list_t), orders = c("ymd HMS", "ymd HM", "ymd"), tz = "UTC")),
  hour = hour(parse_date_time(names(list_t), orders = c("ymd HMS", "ymd HM", "ymd"), tz = "UTC")),
  minute = minute(parse_date_time(names(list_t), orders = c("ymd HMS", "ymd HM", "ymd"), tz = "UTC")),
  second = second(parse_date_time(names(list_t), orders = c("ymd HMS", "ymd HM", "ymd"), tz = "UTC"))
)


describe(unlist(list_t)) # check t0 values
unique(time_lookup$t0_date) # check unique t0 values
head(time_lookup)
selected_points$t0_date <- time_lookup$t0_date[match(selected_points$t0, time_lookup$t0)]

# if there is "Y-M-D" change to "Y-M-D 00:00:00"
selected_points$t0_date <- ifelse(
  nchar(format(selected_points$t0_date, "%Y-%m-%d %H:%M:%S")) == 10,
  paste0(format(selected_points$t0_date, "%Y-%m-%d"), " 00:00:00"),
  format(selected_points$t0_date, "%Y-%m-%d %H:%M:%S")
)

datetimes <- selected_points$t0_date



# plot histogram of t0 months in string character format in english
library(ggplot2)
library(lubridate)
library(dplyr)

datetimes <- as.POSIXct(datetimes, format = "%Y-%m-%d %H:%M:%S", tz = "UTC")
# Extract months from datetimes
df_months <- data.frame(month = month(datetimes, label = TRUE, locale = "en_US.UTF-8"))
# in english
df_months$month <- factor(df_months$month, levels = month.abb)
ggplot(df_months, aes(x = month)) +
  geom_bar(fill = btfgreen, alpha = 0.5, color = "#5f5d5d") +
  btf_theme +
  xlab("Month") +
  ylab("Count") 

# save plot
filename <- paste(im_folder, "optim/comephore/episodes/t0_months_histogram_q",
                  q * 100, "_delta", delta, "_dmin", min_spatial_dist,
                  ".png", sep = "")

ggsave(filename, width = 20, height = 15, units = "cm")

datetimes <- unique(datetimes)
datetimes_df <- data.frame( datetime = datetimes)

filename <- paste0(data_folder, "comephore/t0_episodes",
           "_q", q * 100, "_delta", delta,
           "_dmin", min_spatial_dist, ".csv")

write.csv(datetimes_df, file = filename, row.names = FALSE)

datetimes_df <- read_csv(filename, col_types = cols(
  datetime = col_datetime(format = "%Y-%m-%d %H:%M:%S")
))


# Threshold histogram
df_threshold <- data.frame(u_s0 = selected_points$u_s0)
breaks <- seq(floor(min(df_threshold$u_s0)), ceiling(max(df_threshold$u_s0)), by = 0.1)

ggplot(df_threshold, aes(x = u_s0)) +
  geom_histogram(breaks = breaks, fill = btfgreen, alpha = 0.5, color = "#5f5d5d") +
  btf_theme +
  xlab(TeX(paste0("Threshold for quantile $q = ", q, "$"))) +
  ylab("Count")
filename <- paste(im_folder, "optim/comephore/5km_threshold_histogram_q",
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
u_list <- selected_points$u_s0
list_episodes_points <- get_extreme_episodes(selected_points, comephore_subset,
                                     episode_size = episode_size, unif = FALSE,
                                     beta = 0)

list_episodes <- list_episodes_points$episodes
episode <- list_episodes[[1]]
u <- u_list[[1]]
# check the episode
head(episode)
library(ggplot2)
library(reshape2)  # for melting wide data to long format
 
# Convert matrix to data frame
index <- 1
sort(t0_list)
which(t0_list == t0_list[index])
episode_test <- list_episodes[[index]]
df_episode <- as.data.frame(episode_test)
df_episode$Time <- 0:(nrow(df_episode)-1)  # Add a time column
s0_list[index]
u_episode <- selected_points$u_s0[index]
# t0_episode <- t0_list[index]
# Convert from wide to long format
df_long <- melt(df_episode, id.vars = "Time")
head(df_long)
colnames(df_long) <- c("Time", "Pixel", "Value")
ggplot(df_long, aes(x = Time, y = Value, group = Pixel)) +
  geom_line(color = btfgreen) +
  geom_hline(yintercept = u_episode, color = "red", linetype = "dashed") +
  theme_minimal()

filename <- paste(im_folder, "optim/comephore/extreme_episode", index, "_min", min_spatial_dist,
                  "km_q", q*100, "_delta_", delta, ".png", sep = "")
# filename <- "test.png"
ggsave(filename, width = 20, height = 15, units = "cm")


list_episodes_unif_points <- get_extreme_episodes(selected_points, comephore_subset,
                                      episode_size = episode_size, unif = TRUE)

list_episodes_unif <- list_episodes_unif_points$episodes

s0_list <- selected_points$s0
t0_list <- selected_points$t0
u_list <- selected_points$u_s0

library(parallel)

tau_vect <- 0:10
tmax <- max(tau_vect)
df_coords <- as.data.frame(grid_coords_km)
# Compute the lags and excesses for each conditional point
list_results <- mclapply(1:length(s0_list), function(i) {
  s0 <- s0_list[i]
  s0_coords <- df_coords[s0, ]
  # t0 <- t0_list[i]
  u <- u_list[i]
  episode <- list_episodes[[i]] / u
  ind_t0_ep <- 0 # index of t0 in the episode
  lags <- get_conditional_lag_vectors(df_coords, s0_coords, ind_t0_ep,
                                tau_vect, latlon = FALSE)
  # lags$hnorm <- lags$hnorm / 1000 # convert to km
  excesses <- empirical_excesses_rpar(episode, 1, lags, t0 = ind_t0_ep)
  lags$tau <- lags$tau # convert tau from hours to seconds
  # lags$hnorm <- lags$hnorm / 1000 # convert to km
  # lags$hx <- lags$hx / 1000 # convert hx from m to km
  # lags$hy <- lags$hy / 1000 # convert hy from m
  list(lags = lags, excesses = excesses)
}, mc.cores = detectCores() - 1)

list_lags <- lapply(list_results, `[[`, "lags")
list_excesses <- lapply(list_results, `[[`, "excesses")

# save workspace
# save.image("workspace.RData")
load("workspace.RData")

# load("workspace.RData")
s0 <- s0_list[index]
s0_coords <- df_coords[s0, ]
excesses <- list_excesses[[10]]
excesses$kij
excesses[1:20, ]
df_lags <- list_lags[[1]]
tail(df_lags)

sum_excess_vect <- c()
for (i in 1:length(list_excesses)) {
  excesses <- list_excesses[[i]]
  sum_excess_vect <- c(sum_excess_vect, sum(excesses$kij))
}


# show distribution of sum of excesses
library(ggplot2)
df_excesses <- data.frame(sum_excess = sum_excess_vect)
ggplot(df_excesses, aes(x = sum_excess)) +
  geom_histogram(bins = 30, fill = btfgreen, alpha = 0.5, color = "#5f5d5d") +
  btf_theme +
  xlab("Sum of excesses") +
  ylab("Count") +
  ggtitle("Distribution of sum of excesses")

ggsave(filename = paste0(im_folder, "optim/comephore/sum_excesses_histogram_q_",
                         q * 100, "_min_spatial_dist_", min_spatial_dist,
                         "km_tmax_", tmax, "h_delta_", delta, ".png"),
       width = 20, height = 15, units = "cm")


# ADD WIND DATA ################################################################

# wind_per_episode <- Map(compute_wind_episode_vector_mean, list_episodes,
#                MoreArgs = list(wind_df = wind_mtp))

# wind_ep_df <- do.call(rbind, wind_per_episode)
# head(wind_ep_df)
# tail(wind_ep_df)

# # Plot mean vector wind rose
# wind_ep_df$DD_vector_mean <- as.numeric(wind_ep_df$DD_vector_mean)
# ggplot(wind_ep_df, aes(x = DD_vector_mean, y = FF_vector_mean)) +
#   geom_point(color = btfgreen, alpha = 0.5) +
#   coord_polar() +
#   xlab("Wind direction (degrees)") +
#   ylab("Wind speed (m/s)") +
#   ggtitle("Mean wind vector per episode") +
#   theme_minimal()


# # Bin wind direction into 30° sectors
# wind_ep_df <- wind_ep_df %>%
#   mutate(dir_bin = cut(DD_vector_mean,
#                        breaks = seq(0, 360, by = 45),  # 8 bins for 360°
#                        labels = c("N", "NE", "E", "SE", "S", "SW", "W", "NW"),
#                        include.lowest = TRUE))

# # Aggregate mean speed by direction bin
# agg_df <- wind_ep_df %>%
#   group_by(dir_bin) %>%
#   summarise(mean_speed = mean(FF_vector_mean, na.rm = TRUE))

# # Plot wind rose style bar chart
# ggplot(agg_df, aes(x = dir_bin, y = mean_speed, fill = mean_speed)) +
#   geom_bar(stat = "identity", width = 1, color = "black") +
#   coord_polar(start = -pi/8) +  # Align North to top
#   scale_fill_gradient(
#     low = "#cbe6d7",
#     high = btfgreen,
#     name = "Force (m/s)"  # <-- Legend title here
#   ) +
#   labs(
#     title = "",
#     x = "Direction",
#     y = NULL
#   ) +
#   theme_minimal() +
#   theme(axis.text.x = element_text(size = 10))

# # save plot
# filename <- paste(im_folder, "optim/comephore/wind_rose_mean_vector_per_episode.png", sep = "")
# ggsave(filename, width = 20, height = 15, units = "cm")

# ADD ADVECTION ESTIMATES ######################################################

adv_filename <- paste0(data_folder, "comephore/adv_estim/advection_results.csv")
adv_df <- read.csv(adv_filename, sep = ",")
head(adv_df)

# get only matching episodes from selected_points
matching_indices <- match(selected_points$t0_date, adv_df$t0)
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
  t0_date <- selected_episodes$t0_date[i]
  adv_row <- adv_df[adv_df$t0 == t0_date, ]
  if (nrow(adv_row) == 0) {
    print(i)
    print(paste("No advection data found for t0_date =", t0_date))
  } else {
    # if there are multiple rows, take the first one
    adv_row <- adv_row[1, ]
    adv_x <- adv_row$mean_dx_mps
    adv_y <- adv_row$mean_dy_mps
    selected_episodes$adv_x[i] <- adv_x
    selected_episodes$adv_y[i] <- adv_y
  }
}

head(selected_episodes)
# Remove episodes with NA advection values
ind_NA_adv <- which(is.na(selected_episodes$adv_x) | is.na(selected_episodes$adv_y))
selected_episodes_nona <- selected_episodes[-ind_NA_adv, ]

wind_df <- selected_episodes[, c("adv_x", "adv_y")]
colnames(wind_df) <- c("vx", "vy")
length(wind_df$vx) # should be the same as number of episodes
# wind_df$vx <- - wind_df$vx  # invert x component
# wind_df$vy <- - wind_df$vy  # invert y component
# convert m/s to km/h
wind_df$vx <- wind_df$vx * 3.6  # convert to km/h
wind_df$vy <- wind_df$vy * 3.6  # convert to km/h
# OPTIMIZATION #################################################################

# Compute the wind vector for each episode (-FF because it's the wind direction)
# sin and cos are shifted because 0 degree means North
# wind_ep_df$vx <- -wind_ep_df$FF_vector_mean * sin(wind_ep_df$DD_vector_mean * pi / 180)
# wind_ep_df$vy <- -wind_ep_df$FF_vector_mean * cos(wind_ep_df$DD_vector_mean * pi / 180)

# if values of vx or vy are really close to 0, set them to 0
# wind_ep_df$vx[abs(wind_ep_df$vx) < 1e-08] <- 0
# wind_ep_df$vy[abs(wind_ep_df$vy) < 1e-08] <- 0

# wind_df <- wind_ep_df[, c("vx", "vy")]
# which episode has NA wind values
ind_NA <- which(is.na(wind_df$vx))
ind_NA <- ind_NA_adv
wind_opt <- wind_df
if (any(ind_NA > 0)) {
  # remove these episodes
  # wind_opt <- wind_df[-ind_NA, ]
  episodes_opt <- list_episodes[-ind_NA]
  lags_opt <- list_lags[-ind_NA]
  excesses_opt <- list_excesses[-ind_NA]
} else {
  wind_opt <- wind_df
  episodes_opt <- list_episodes
  lags_opt <- list_lags
  excesses_opt <- list_excesses
}


eta1 <- 1.5
eta2 <- 1
i <- 1
excesses <- list_excesses[[i]]
lags <- list_lags[[i]]
head(lags)
# Choisir advection (varie selon wind_df)
if (!all(is.na(wind_df))) {
  adv <- as.numeric(wind_df[i, ])
  # if adv_x or adv_y are really close
} else {
  adv <- c(0, 0)
}

adv_x <- adv[1]  # vx
adv_y <- adv[2]  # vy
params_adv <- c(beta1, beta2, alpha1, alpha2, adv_x, adv_y)  # add adv to params
adv <- params_adv[5:6]  # adv_x, adv_y
# Obtenir la table chi complète
chi_df <- theoretical_chi(params = params_adv, df_lags = lags,
                          latlon = FALSE, directional = TRUE,
                          convert_in_hours = F, convert_in_km = F)
hnorm <- chi_df$hnorm
hnormV <- chi_df$hnormV
summary(hnorm - hnormV)
library(ggplot2)

# 1. Vecteurs de déplacement (hx, hy)
ggplot(chi_df, aes(x = hx, y = hy)) +
  geom_point(alpha = 0.6) +
  labs(title = "Vecteurs de déplacement hx/hy",
       x = "hx (km)", y = "hy (km)") +
  theme_minimal()

# 2. Variogramme (vario) vs distance and for different time lags as factor
ggplot(chi_df, aes(x = hnormV, y = vario, color = factor(tau))) +
  geom_point(alpha = 0.6) +
  labs(title = "",
       x = "Distance (km)", y = "Variogramme",
       color = "Time Lag (tau)") +
  theme_minimal()

# 3. Chi vs distance
ggplot(chi_df, aes(x = hnormV, y = chi, color = factor(tau))) +
  geom_point(alpha = 0.6) +
  labs(title = "Chi vs Distance", x = "Distance (km)", y = "Chi") +
  theme_minimal()

# 4. Chi vs vario (fonction de transformation)
ggplot(chi_df, aes(x = vario, y = chi)) +
  geom_point(color = "purple") +
  labs(title = "Chi vs Variogramme", x = "Variogramme", y = "Chi") +
  theme_minimal()

print(adv) # Advection utilisée pour cet épisode
summary(chi_df$hx)
summary(chi_df$hy)
summary(chi_df$chi)
summary(chi_df$vario)

ggplot(chi_df, aes(x = s1x, y = s1y)) +
  geom_segment(aes(xend = s2x, yend = s2y), arrow = arrow(length = unit(0.1, "cm")), alpha = 0.3) +
  labs(title = "Displacement between points",
       x = "s1x (km)", y = "s1y (km)") +
  coord_fixed()

for (i in c(1, 210)) {
  adv <- as.numeric(wind_df[i, ])
  adv_x <- adv[1]  # vx
  adv_y <- adv[2]  # vy
  params_adv <- params_adv <- c(beta1, beta2, alpha1, alpha2, adv_x, adv_y)
  chi_df <- theoretical_chi(params_adv, list_lags[[i]],
                            latlon = FALSE, directional = TRUE,
                            convert_in_hours = FALSE, convert_in_km = FALSE)
  
  chi_df$episode <- i
  if (i == 1) {
    chi_all <- chi_df
  } else {
    chi_all <- rbind(chi_all, chi_df)
  }
}

ggplot(chi_all, aes(x = hnormV, y = chi, color = factor(episode))) +
  geom_point(alpha=0.5) +
  labs(title = "Comparaison de Chi entre épisodes", x = "Distance (km)", y = "Chi", color = "Épisode") +
  theme_minimal()


init_param <- c(beta1, beta2, alpha1, alpha2, 1, 1)
nrow(wind_opt)
adv_df <- data.frame(vx = 0, vy = 0)
length(episodes_opt) # should be the same
result <- optim(par = init_param, fn = neg_ll_composite,
        list_lags = lags_opt, list_episodes = episodes_opt,
        list_excesses = excesses_opt, hmax = 5,
        wind_df = wind_opt,
        latlon = FALSE,
        directional = TRUE,
        fixed_eta1 = FALSE,
        fixed_eta2 = FALSE,
        convert_in_hours = FALSE,
        convert_in_km = FALSE,
        method = "L-BFGS-B",
        lower = c(1e-08, 1e-08, 1e-08, 1e-08, 1e-08, 1e-08),
        upper = c(10, 10, 1.999, 1.999, 10, 10),
        control = list(maxit = 10000,
                      trace = 1),
        hessian = T)

result
test_init <-  result$par

# From results get CI 
# Compute the standard errors from the Hessian
hessian <- result$hessian
if (is.null(hessian)) {
  stop("Hessian is NULL, cannot compute standard errors.")
} 
# Supposons que tu connais les indices des paramètres fixés
# fixed_indices <- c(6)

# # Tu gardes les lignes/colonnes des paramètres variables
# var_indices <- setdiff(seq_len(nrow(hessian)), fixed_indices)

# # Sous-matrice hessienne restreinte aux paramètres libres
# hessian <- hessian[var_indices, var_indices]

# Calcul des erreurs standard uniquement pour les paramètres libres

inv_hessian <- solve(hessian)

diag_var <- diag(inv_hessian)
diag_var[diag_var <= 0] <- 0

# Calcul des erreurs standard pour les paramètres libres
se <- sqrt(diag_var)

# Create a data frame with the results
df_results <- data.frame(
  Parameter = c("beta1", "beta2", "alpha1", "alpha2", "vx", "vy"),
  Estimate = result$par,
  StdError = se,
  LowerCI = result$par - 1.96 * se,
  UpperCI = result$par + 1.96 * se
)

# do a latex table with kable and 


# > result
# $par
# [1] 0.1129827 1.0965668 1.0217060 0.6181883 1.0000000 1.0000000

# $value
# [1] 106061.9


# EMPIRICAL VARIOGRAM ##########################################################
# library(sp)
# library(spacetime)
# library(gstat)
# library(reshape2)

# coords_pixels$pixel <- rownames(coords_pixels)

# df_comephore$date <- as.POSIXct(df_comephore$date)

# df_long <- melt(df_comephore, id.vars = "date", variable.name = "pixel", value.name = "rain")
# head(df_long)
# colnames(df_long) <- c("date", "pixel", "rain")
# # df_long <- df_long[df_long$rain > 0, ]

# # if(nrow(df_long) == 0) stop("Aucune valeur de pluie non nulle détectée.")

# df_long <- merge(df_long, coords_pixels, by = "pixel")

# coordinates(df_long) <- ~ Longitude + Latitude

# sp_pts <- SpatialPoints(coords_pixels[, c("Longitude", "Latitude")])

# time_pts <- sort(unique(df_long$date))


# library(tidyr)
# df_full <- expand.grid(date = time_pts, pixel = coords_pixels$pixel)
# df_full <- merge(df_full, df_long@data, by = c("date", "pixel"), all.x = TRUE)

# df_full$rain[is.na(df_full$rain)] <- 0

# df_full <- merge(df_full, coords_pixels, by = "pixel")

# sp_pts_full <- SpatialPoints(df_full[, c("Longitude", "Latitude")])

# stfdf <- STFDF(sp = sp_pts, time = time_pts,
#                data = data.frame(rain = df_full$rain))
# unique(diff(time(stfdf)))

# vg_st <- variogramST(rain ~ 1, data = stfdf, assumeRegular = TRUE,
#                      width = 1000, tlags = 0)

# vg_st
# plot(vg_st)


# VARIOGRAM PLOTS ##############################################################

# compute variogram with parameters
df_result <- data.frame(beta1 = result$par[1],
                        beta2 = result$par[2],
                        alpha1 = result$par[3],
                        alpha2 = result$par[4],
                        eta1 = result$par[5],
                        eta2 = result$par[6])

beta1 <- df_result$beta1
beta2 <- df_result$beta2
alpha1 <- df_result$alpha1
alpha2 <- df_result$alpha2
eta1 <- df_result$eta1
eta2 <- df_result$eta2

library(dplyr)
# Compute advection for each row
adv_df_ep <- wind_df %>%
  mutate(
    adv_x = (abs(vx)^eta1) * sign(vx) * eta2,
    adv_y = (abs(vy)^eta1) * sign(vy) * eta2
  )

# remove NA values
adv_df_ep <- adv_df_ep[!is.na(adv_df_ep$adv_x), ]
head(adv_df_ep)

episode <- list_episodes[[1]]

lags <- list_lags[[1]]

wind_ep_1 <- wind_df[1, ]

adv_ep_1 <- adv_df_ep[1, ]

head(df_comephore)

# Vecteur d'advection pour l'épisode 1
vx <- wind_ep_1$vx
vy <- wind_ep_1$vy

adv_x <- (abs(vx)^eta1) * sign(vx) * eta2
adv_y <- (abs(vy)^eta1) * sign(vy) * eta2

lags <- lags %>%
  mutate(
    hx = s2x - s1x,
    hy = s2y - s1y,
    hnorm = sqrt(hx^2 + hy^2)
  )
# Recalage des lags
lags_recal <- lags %>%
  mutate(
    hx_recal = hx - adv_x * tau,
    hy_recal = hy - adv_y * tau,
    hnorm_recal = sqrt(hx_recal^2 + hy_recal^2),
    gamma_model = beta1 * (abs(hx_recal)^alpha1) + beta1 * (abs(hy_recal)^alpha1) + beta2 * (abs(tau)^alpha2)
  )


ggplot(lags_recal %>% filter(tau == 10), aes(x = hx_recal, y = hy_recal, fill = gamma_model)) +
  geom_tile() +
  coord_equal() +
  labs(title = "Directional variogram for tau = 0", x = "hx", y = "hy")



set.seed(42)  # Pour reproductibilité

# Tirer 10 épisodes aléatoires
episode_ids <- c(1, sample(seq_along(list_lags), 10))
tau_value <- 0
lags_recal_all <- purrr::map_dfr(episode_ids, function(i) {
  # Données de l'épisode i
  lags <- list_lags[[i]]
  wind <- wind_df[i, ]  # même index

  # Calcul advection (si absents)
  adv_x <- (abs(wind$vx)^eta1) * sign(wind$vx) * eta2
  adv_y <- (abs(wind$vy)^eta1) * sign(wind$vy) * eta2
  
  # Recalcul des hx/hy recalés
  lags_recal <- lags %>%
    mutate(
      hx_recal = s1x - s2x + adv_x * tau,
      hy_recal = s1y - s2y + adv_y * tau,
      gamma_model = beta1 * (abs(hx_recal)^alpha1) +
                    beta1 * (abs(hy_recal)^alpha1) +
                    beta2 * (abs(tau)^alpha2),
      episode_id = paste0("E", i)
    )
  p <- ggplot(lags_recal %>% filter(tau == tau_value),
             aes(x = hx_recal, y = hy_recal, fill = gamma_model)) +
    geom_tile() +
    coord_equal() +
      labs(
    title = paste("Episode", i, "at", tau_value, "hour(s)"),
    x = expression(h[x] - tau * V),
    y = expression(h[y] - tau * V),
    fill = expression(gamma(h, tau))
  ) +
  scale_fill_viridis_c()

  #save plot
  foldername <- paste0(im_folder, "optim/comephore/vario_estim/q_", q * 100, 
    "_mindist_", min_spatial_dist, "km_delta_", delta, "h/")
  filename <- paste(foldername, "variogram_ep_", i, "_tau", tau_value, ".png", sep = "")
  ggsave(filename, plot = p, width = 20, height = 15, units = "cm")
  return(lags_recal)
})

################################################################################

library(tidyverse)
library(viridis)

episode_id <- 348
tau_values <- 0:5
pixel_size_km <- 1  # Taille du pixel = 1 km

lags <- list_lags[[episode_id]]
wind <- wind_df[episode_id, ]

adv_x <- (abs(wind$vx)^eta1) * sign(wind$vx) * eta2
adv_y <- (abs(wind$vy)^eta1) * sign(wind$vy) * eta2

# Calcul hx/hy recalés pour tous les tau
lags_all <- map_dfr(tau_values, function(tau) {
  lags %>%
    mutate(
      tau = tau,
      hx_recal = s1x - s2x + adv_x * tau,
      hy_recal = s1y - s2y + adv_y * tau,
      gamma_model = beta1 * (abs(hx_recal)^alpha1) +
                    beta1 * (abs(hy_recal)^alpha1) +
                    beta2 * (abs(tau)^alpha2)
    )
})

# Bornes globales
hx_min <- floor(min(lags_all$hx_recal))
hx_max <- ceiling(max(lags_all$hx_recal))
hy_min <- floor(min(lags_all$hy_recal))
hy_max <- ceiling(max(lags_all$hy_recal))

# Nombre pixels (arrondi au supérieur)
epsilon <- 10  # petit dépassement pour inclure la borne max

n_pixels_x <- ceiling((hx_max - hx_min) / pixel_size_km)
n_pixels_y <- ceiling((hy_max - hy_min) / pixel_size_km)

hx_bins <- seq(hx_min, hx_min + n_pixels_x * pixel_size_km + epsilon, by = pixel_size_km)
hy_bins <- seq(hy_min, hy_min + n_pixels_y * pixel_size_km + epsilon, by = pixel_size_km)


lags_all_binned <- lags_all %>%
  mutate(
    hx_bin = cut(hx_recal, breaks = hx_bins, include.lowest = TRUE, labels = FALSE),
    hy_bin = cut(hy_recal, breaks = hy_bins, include.lowest = TRUE, labels = FALSE)
  ) %>%
  filter(!is.na(hx_bin) & !is.na(hy_bin)) %>%
  mutate(
    hx_center = hx_bins[hx_bin] + pixel_size_km / 2,
    hy_center = hy_bins[hy_bin] + pixel_size_km / 2
  )

# Grille complète pour chaque tau
full_grid <- expand.grid(
  hx_center = hx_bins[-length(hx_bins)] + pixel_size_km / 2,
  hy_center = hy_bins[-length(hy_bins)] + pixel_size_km / 2,
  tau = tau_values
)

# Moyenne gamma par pixel et tau
lags_mean <- lags_all_binned %>%
  group_by(tau, hx_center, hy_center) %>%
  summarise(gamma_mean = mean(gamma_model, na.rm = TRUE), .groups = "drop")

lags_mean_complete <- full_grid %>%
  left_join(lags_mean, by = c("tau", "hx_center", "hy_center"))

# Palette commune
gamma_min <- quantile(lags_mean_complete$gamma_mean, 0.01, na.rm = TRUE)
gamma_max <- quantile(lags_mean_complete$gamma_mean, 0.99, na.rm = TRUE)

plots <- map(tau_values, function(tau_val) {
  ggplot(lags_mean_complete %>% filter(tau == tau_val),
         aes(x = hx_center, y = hy_center, fill = gamma_mean)) +
    geom_raster() +
    coord_equal() +
    scale_fill_viridis_c(limits = c(gamma_min, gamma_max), oob = scales::squish, na.value = "transparent") +
    labs(
      title = paste("Épisode", episode_id, "- τ =", tau_val),
      x = expression(h[x] - tau * V),
      y = expression(h[y] - tau * V),
      fill = expression(gamma(h, tau))
    ) +
    theme_minimal()
})

# Sauvegarde des plots
foldername <- paste0(im_folder, "optim/comephore/vario_estim/q_", q * 100, 
                     "_mindist_", min_spatial_dist, "km_delta_", delta, "h/")
dir.create(foldername, recursive = TRUE, showWarnings = FALSE)

for (i in seq_along(plots)) {
  filename <- paste0(foldername, "variogram_ep_", episode_id, "_tau_", tau_values[i], ".png")
  ggsave(filename, plot = plots[[i]], width = 20, height = 15, units = "cm")
}



# Plot estimation of advection

adv_df_ep <- wind_df %>%
  mutate(
    adv_x = (abs(vx)^eta1) * sign(vx) * eta2,
    adv_y = (abs(vy)^eta1) * sign(vy) * eta2
  )

adv_df_ep <- adv_df_ep %>%
  mutate(
    norm = sqrt(adv_x^2 + adv_y^2),
    unit_x = adv_x / norm,
    unit_y = adv_y / norm
  )

ggplot(adv_df_ep, aes(x = adv_x, y = adv_y)) +
  geom_segment(aes(xend = adv_x + 0.2 * unit_x, yend = adv_y + 0.2 * unit_y),
               arrow = arrow(length = unit(0.2, "cm")), color = btfgreen) +
  coord_fixed() +
  labs(title = "Normalized advection directions",
       x = "Adv_x (m/s)", y = "Adv_y (m/s)") +
  theme_minimal()

ggplot(adv_df_ep, aes(x = 0, y = 0)) +
  geom_segment(aes(xend = adv_x, yend = adv_y),
               arrow = arrow(length = unit(0.2, "cm")), color = btfgreen) +
  coord_fixed() +
  labs(title = "",
     x = expression(v[x]~"(m/s)"),
    y = expression(v[y]~"(m/s)")
  ) +
  ylim(c(-500, 500)) +
  xlim(c(-500, 1500)) +
  theme_minimal()

# save plot
filename <- paste(im_folder, "optim/comephore/advection_vector_field.pdf", sep = "")
ggsave(filename, width = 20, height = 15, units = "cm")



library(dplyr)

adv_df_ep <- wind_df %>%
  mutate(
    adv_x = (abs(vx)^eta1) * sign(vx) * eta2,
    adv_y = (abs(vy)^eta1) * sign(vy) * eta2,
    speed = sqrt(adv_x^2 + adv_y^2) / 3.6,  # convert km/h to m/s
    direction = (atan2(adv_y, adv_x) * 180 / pi) %% 360  # direction en degrés
  )


custom_breaks <- c(0, 0.5, 1, 2, 4, 6, 8, 10, Inf)  # vitesse in m/s
custom_labels <- c(
  ">0-0.5 m/s",
  "0.5-1 m/s",
  "1-2 m/s",
  "2-4 m/s",
  "4-6 m/s",
  "6-8 m/s",
  "8-10 m/s",
  ">10 m/s"
)

adv_df_ep <- adv_df_ep %>%
  mutate(
    speed_bin = cut(speed,
                    breaks = custom_breaks,
                    labels = custom_labels,
                    include.lowest = TRUE,
                    right = FALSE)
  )


# get cardinals
adv_df_ep <- adv_df_ep %>%
  mutate(
    dir_bin = cut(direction,
                  breaks = seq(0, 360, by = 45),  # 8 bins for 360°
                  labels = c("N", "NE", "E", "SE", "S", "SW", "W", "NW"),
                  include.lowest = TRUE)
  )


# where are na
na_indices <- which(is.na(adv_df_ep$dir_bin) | is.na(adv_df_ep$speed_bin))

# remove na
adv_df_ep <- adv_df_ep %>%
  filter(!is.na(dir_bin) & !is.na(speed_bin))

# remove 0 advection values
adv_df_ep <- adv_df_ep %>%
  filter(!(adv_x == 0 & adv_y == 0))
# number of zero advection values
n_zero_adv <- nrow(wind_df) - nrow(adv_df_ep)

custom_palette <- c(
  "#dceef8",  # très clair bleu
  "#b3cde3",  # bleu pastel
  "#89abc7",  # bleu moyen
  "#5e8db1",  # bleu soutenu
  "#4f74b1",  # bleu foncé
  "#8856a7",  # violet clair
  "#6a3d9a",  # violet soutenu
  "#4d004b"   # violet très foncé
)


# Plot wind rose with speed bins
ggplot(adv_df_ep, aes(x = dir_bin, fill = speed_bin)) +
  geom_bar(width = 0.7) +
  scale_fill_manual(
    values = custom_palette,
    name = "Speed (m/s)"
  ) +
  coord_polar(start = -pi/8) +
  labs(title = "",
       x = "Direction",
       y = "Count") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(size = 10),
    legend.key.height = unit(0.6, "cm"),
    legend.key.width = unit(0.4, "cm")
  )

# save plot
filename <- paste(im_folder, "optim/comephore/adv_estimates/wind_rose_speed_per_direction.pdf", sep = "")
ggsave(filename, width = 20, height = 15, units = "cm")


colnames(selected_episodes) <- c("s0", "t0", "u_s0", "t0_date", "adv_x_ms", "adv_y_ms")
# add adv estimates in selected_episodes
selected_episodes$adv_hat_x_ms <- adv_df_ep$adv_x
selected_episodes$adv_hat_y_ms <- adv_df_ep$adv_y

foldername <- paste0(data_folder, "comephore/optim_adv/")
dir.create(foldername, recursive = TRUE, showWarnings = FALSE)
wind_df_filename <- paste(data_folder, "/comephore/optim_adv/adv_df.csv", sep = "")
write.csv(selected_episodes, wind_df_filename, row.names = FALSE)
