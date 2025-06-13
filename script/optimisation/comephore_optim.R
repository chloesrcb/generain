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

load("workspace.RData")
# library(generain)

# LOAD DATA ####################################################################
filename_com <- paste0(data_folder, "comephore/zoom_5km.csv")
comephore_raw <- read.csv(filename_com, sep = ",")
filename_loc <- paste0(data_folder, "comephore/loc_px_zoom_5km.csv")
loc_px <- read.csv(filename_loc, sep = ",")

df_comephore <- comephore_raw
head(df_comephore)

# Take only data after 2007
colnames(df_comephore)[1] <- "date"
df_comephore <- df_comephore[df_comephore$date >= "2008-01-01", ]

# put date in index
rownames(df_comephore) <- df_comephore$date
comephore <- df_comephore[-1] # remove dates column
# Get distances matrix
dist_mat <- get_dist_mat(loc_px)
nsites <- nrow(loc_px)


# # get wind data
filename_wind <- paste0(data_folder, "wind/data_gouv/wind_mtp.csv")
wind_mtp <- read.csv(filename_wind)

# Convert datetime to POSIXct
wind_mtp$datetime <- as.POSIXct(wind_mtp$datetime,
                                format = "%Y-%m-%d %H:%M:%S", tz = "UTC")

# Apply function to the DD column
wind_mtp$cardDir <- sapply(wind_mtp$DD, convert_to_cardinal)
wind_mtp$cardDir <- as.character(wind_mtp$cardDir)  # Ensure it's character
wind_mtp$cardDir[is.na(wind_mtp$DD)] <- NA
# Check if NA values are properly handled
summary(wind_mtp)

head(wind_mtp$cardDir)

# DATA WITHOUT 0 ###############################################################

# Remove rows with all columns equal to 0
# comephore_no0 <- comephore[rowSums(comephore) != 0, ]

# nrow(comephore)
# nrow(comephore_no0)

# CHOOSE QUANTILE ##############################################################

# get a matrix of high quantiles for all pair
# q <- 0.99 # quantile
# list_count_quant <- quantile_matrix(q, comephore, qlim = TRUE, zeros = TRUE,
#                                     count_min = 150) # with removing zeros
# quant_mat <- list_count_quant[1][[1]]
# count_mat <- list_count_quant[2]
colnames(comephore)
head(comephore)
library(ggplot2)
par(mfrow = c(1, 2))
comephore_pair <- comephore[, c(12, 2)]
comephore_pair_no0 <- comephore_pair[rowSums(comephore_pair) > 0, ]
chiplot(comephore_pair_no0, xlim = c(0.85, 1), ylim1 = c(0.5, 1), which = 1,
        qlim = c(0.85, 0.999))
abline(v = 0.96, col = "red", lty = 2)

# Same with all zeros
comephore_pair <- comephore[, c(12, 2)]
chiplot(comephore_pair, xlim = c(0.99, 1), ylim1 = c(0.5, 1), which = 1,
        qlim = c(0.99, 0.9999))
abline(v = 0.998, col = "red", lty = 2)

# count conjoint excesses
q <- 0.97
# uniformize the data
n <- nrow(comephore_pair_no0)
data_unif <- cbind(rank(comephore_pair_no0[, 1]) / (n + 1),
                   rank(comephore_pair_no0[, 2]) / (n + 1))

count_excesses <- sum(data_unif[, 1] > q & data_unif[, 2] > q)
print(count_excesses)

# count conjoint excesses
q <- 0.9975
# uniformize the data
n <- nrow(comephore_pair)
data_unif <- cbind(rank(comephore_pair[, 1]) / (n + 1),
                          rank(comephore_pair[, 2]) / (n + 1))

count_excesses <- sum(data_unif[, 1] > q & data_unif[, 2] > q)
print(count_excesses)

q_no0_spa <- 0.97
# get the corresponding quantile in the data without 0 and with 0
q_value_no0 <- quantile(c(comephore_pair_no0[,1], comephore_pair_no0[,2]), 
                        probs = q_no0_spa, na.rm = TRUE)

q_equiv_prob <- ecdf(c(comephore_pair[,1], comephore_pair[,2]))(q_value_no0)
print(q_equiv_prob)

q_equiv <- quantile(c(comephore[,1], comephore[,15]), probs = q_equiv_prob, 
              na.rm = TRUE)
print(q_equiv)  # quantile and threshold

q_spa <- round(q_equiv_prob, 4)

# Temporal chi
par(mfrow = c(1,2))
rain_nolag <- comephore$p102[1:(length(comephore$p102) - 5)]
rain_lag <- comephore$p102[6:length(comephore$p102)]
comephore_pair <- cbind(rain_nolag, rain_lag)
comephore_pair_no0 <- comephore_pair[rowSums(comephore_pair) > 0, ]
chiplot(comephore_pair_no0, xlim = c(0.8, 0.99), ylim1 = c(0, 1), which = 1,
        qlim = c(0.8, 0.99))
abline(v = 0.94, col = "red", lty = 2)

comephore_pair <- cbind(rain_nolag, rain_lag)
# comephore_pair <- comephore_pair[rowSums(comephore_pair) > 0, ]
chiplot(comephore_pair, xlim = c(0.98, 1), ylim1 = c(0, 1), which = 1,
        qlim = c(0.98, 0.999))
abline(v = 0.995, col = "red", lty = 2)

# count conjoint excesses
q_no0_spa <- 0.97

q_no0_temp <- 0.94
# get the corresponding quantile in the data without 0 and with 0
q_value_no0 <- quantile(c(comephore_pair_no0[,1], comephore_pair_no0[,2]),
                        probs = q_no0_temp, na.rm = TRUE)

q_equiv_prob <- ecdf(c(comephore_pair[,1], comephore_pair[,2]))(q_value_no0)
print(q_equiv_prob)

q_temp <- round(q_equiv_prob, 4)
q_equiv <- quantile(c(comephore[,1], comephore[,15]), probs = q_equiv_prob, 
              na.rm = TRUE)
print(q_equiv)  # quantile and threshold

# BUHL WLSE ####################################################################

# Temporal chi
tmax <- 10
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
filename <- paste(im_folder, "WLSE/comephore/full_temporal_chi_boxplot_", q_no0_temp,
                 ".png", sep = "")
ggsave(filename, width = 20, height = 15, units = "cm")


# Mean of chi
tmax = 10
q_no0_temp = 0.9
chimat_dt_mean <- temporal_chi(comephore, tmax, quantile = q_no0_temp,
                               mean = TRUE, zeros = F)
df_chi <- data.frame(lag = c(0:tmax), chi = chimat_dt_mean)
# chitemp_plot <- ggplot(df_chi, aes(x = lag, y = chi)) +
#   geom_point(color = btfgreen) +
#   btf_theme +
#   xlab(TeX(r"($\tau$ (hours))")) +
#   ylab(TeX(r"($\widehat{\chi}(0,\tau)$)"))

# chitemp_plot

df_chi_not0 <- df_chi[df_chi$lag > 0, ]
wlse_temp <- get_estimate_variotemp(df_chi_not0, weights = "exp", summary = TRUE)
print(wlse_temp)
c2 <- wlse_temp[[1]]
beta2 <- wlse_temp[[2]]
alpha2 <- wlse_temp[[3]]
print(beta2)

dftemp <- data.frame(lag = log(df_chi_not0$lag), chi = eta(df_chi_not0$chi))

# remove first row
# dftemp <- dftemp[-1, ]

chitemp_eta_estim <- ggplot(dftemp, aes(x = lag, y = chi)) +
  geom_point(color = btfgreen) +
  btf_theme +
  xlab(TeX(r"($\log(\tau)$)")) +
  ylab(TeX(r"($\eta(\widehat{\chi}(0,\tau))$)")) +
  geom_line(aes(x = lag, y = alpha2 * lag + c2),
            alpha = 0.5, color = "darkred", linewidth = 0.5)

chitemp_eta_estim

# save plot
filename <- paste(im_folder, "WLSE/comephore/full_temporal_chi_eta_estim_",
                q_no0_temp, ".png", sep = "")
ggsave(filename, width = 20, height = 15, units = "cm")


# Get coords
sites_coords <- loc_px[, c("Longitude", "Latitude")]
rownames(sites_coords) <- loc_px$pixel_name



library(sf)
library(dplyr)
library(ggplot2)
library(gridExtra)

sites_coords_sf <- st_as_sf(sites_coords, coords = c("Longitude", "Latitude"), crs = 4326)

sites_coords_sf <- st_transform(sites_coords_sf, crs = 2154)

coords_m <- st_coordinates(sites_coords_sf)

grid_coords <- sites_coords
grid_coords$x_km <- (coords_m[, "X"] - min(coords_m[, "X"])) / 1000
grid_coords$y_km <- (coords_m[, "Y"] - min(coords_m[, "Y"])) / 1000


p1 <- ggplot(grid_coords, aes(x = Longitude, y = Latitude)) +
  geom_point(color = "blue", size = 2) +
  geom_text(aes(label = rownames(grid_coords)), hjust = -0.2, size = 1.5) +
  coord_fixed() +
  ggtitle("GPS coordinates WGS84") +
  theme_minimal()

p2 <- ggplot(grid_coords, aes(x = x_km, y = y_km)) +
  geom_point(color = "red") +
  geom_text(aes(label = rownames(sites_coords)), size = 1.5, hjust = -0.2) +
  coord_fixed() +
  theme_minimal() +
  ggtitle("Transformed coordinates in km")

a <- grid.arrange(p1, p2, ncol = 2)

# save plot
filename <- paste(im_folder, "optim/comephore/coords_transformation.png", sep = "")
ggsave(plot= a, filename = filename, width = 20, height = 15, units = "cm")

# get distance matrix
grid_coords <- grid_coords[, c("x_km", "y_km")]
colnames(grid_coords) <- c("Longitude", "Latitude")
dist_mat <- get_dist_mat(grid_coords, latlon = FALSE)

# Spatial chi
df_dist <- reshape_distances(dist_mat)
# df_dist$value <- ceiling(df_dist$value / 100) * 100  # / 1000 in km
df_dist_km <- df_dist
df_dist_km$value <- ceiling(round(df_dist$value, 1) * 1000 ) /1000
# df_dist_km$value <- df_dist$value / 1000

h_vect <- sort(unique(df_dist_km$value))
h_vect <- h_vect[h_vect > 0]  # remove 0
hmax <- h_vect[10] # 10th value

hmax <- 7

q_no0_spa <- 0.98
chispa_df <- spatial_chi_alldist(df_dist_km, data_rain = comephore,
                              quantile = q_no0_spa, hmax = hmax, zeros = FALSE)

spatial_chi_plot(chispa_df, hmax = hmax, q = q_no0_spa)
# chispa_df$lagspa <- chispa_df$lagspa 
etachispa_df <- data.frame(chi = eta(chispa_df$chi),
                           lagspa = log(chispa_df$lagspa))

chispa_df <- spatial_chi_extremogram(df_dist_km, comephore, q_no0_spa,
                        hmax = hmax, zeros = FALSE)


chispa_plot <- ggplot(chispa_df, aes(lagspa, chi)) +
  btf_theme +
  geom_point(col = btfgreen) +
  xlab(TeX(r"($h$)")) +
  ylab(TeX(r"($\widehat{\chi}(h, 0)$)")) +
  ylim(0, 1)

chispa_plot

# save plot
filename <- paste(im_folder, "WLSE/comephore/full_spatial_chi_", q,
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
  ylab(TeX(r"($\eta(\widehat{\chi}(h, 0))$)")) +
  geom_line(aes(x = lagspa, y = alpha1 * lagspa + c1), alpha = 0.6,
            color = "darkred", linewidth = 0.5)

chispa_eta_estim

# save plot
filename <- paste(im_folder, "WLSE/comephore/full_spatial_chi_eta_estim_", q,
                  ".png", sep = "")
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
save.image("workspace.RData")

# CHOOSE EXTREME EPISODE FOR R-PARETO ##########################################

q <- 0.998 # quantile
# get central site from sites_coords

set_st_excess <- get_spatiotemp_excess(comephore, q)

list_s <- set_st_excess$list_s
unique(unlist(list_s))
list_t <- set_st_excess$list_t

comephore[list_t[[1]], list_s[[1]]]

list_u <- set_st_excess$list_u
list_u[[1]]

# Spatio-temporal neighborhood parameters
min_spatial_dist <- 3 # in km
delta <- 30 # in hours
episode_size <- delta # size of the episode
s0t0_set <- get_s0t0_pairs(grid_coords, comephore,
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
  if (comephore[t0, s0] <= u_s0) {
    stop(paste("Excess is not above threshold for s0 =", s0, "and t0 =", t0))
  }
}


# sites_coords <- grid_coords
# # Step 1: Calculate the centroid
# centroid_lon <- mean(sites_coords$Longitude)
# centroid_lat <- mean(sites_coords$Latitude)

# # Step 2: Calculate Euclidean distance to the centroid
# sites_coords$dist_to_centroid <- sqrt((sites_coords$Longitude - centroid_lon)^2 +
#                                       (sites_coords$Latitude - centroid_lat)^2)

# # Step 3: Find the most central site
# most_central_site <- sites_coords[which.min(sites_coords$dist_to_centroid), ]

# # Output the most central site
# print(most_central_site)

# # remove dist_to_centroid column
# sites_coords <- sites_coords[, -ncol(sites_coords)]

# central_site <- rownames(most_central_site)
# central_site_coords <- most_central_site[1:2]

# # Get the selected episodes
# selected_points <- select_extreme_episodes(sites_coords,
#                                         data = comephore,
#                                         min_spatial_dist = min_spatial_dist,
#                                         episode_size = episode_size,
#                                         n_max_episodes = 10000,
#                                         quantile = q,
#                                         time_ext = 0, latlon = FALSE)

n_episodes <- length(selected_points$s0)
print(n_episodes)
length(unique(selected_points$s0)) # can be same s0
length(unique(selected_points$t0)) # never same t0?
print(min(selected_points$u_s0)) # min threshold
t0_list <- selected_points$t0
s0_list <- selected_points$s0
list_episodes_points <- get_extreme_episodes(selected_points, comephore,
                                     episode_size = episode_size, unif = FALSE,
                                     beta = 0)

list_episodes <- list_episodes_points$episodes
episode <- list_episodes[[1]]
# check the episode
head(episode)
library(ggplot2)
library(reshape2)  # for melting wide data to long format

# Convert matrix to data frame
index <- 10
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
  # labs(title = "Extreme Episode", x = "Time", y = "Value") +
  # annotate("text", x = t0_episode, y = 0, label = expression(t[0]),
  #          color = "red", vjust = 0, size = 7) +
  theme_minimal()

filename <- paste(im_folder, "optim/comephore/extreme_episode", index, "_min", min_spatial_dist,
                  "km_max", tmax, "h_delta_", delta, ".png", sep = "")
# filename <- "test.png"
ggsave(filename, width = 20, height = 15, units = "cm")


library(ggplot2)
library(reshape2)

unique_s0_list <- unique(s0_list)

for (target_s0 in unique_s0_list) {
  matching_indices <- which(s0_list == target_s0)
  
  df_list <- list()
  
  for (i in matching_indices) {
    episode_test <- list_episodes[[i]]
    df_episode <- as.data.frame(episode_test)
    df_episode$Time <- 0:(nrow(df_episode) - 1)

    df_long <- melt(df_episode, id.vars = "Time", variable.name = "Pixel", value.name = "Value")
    colnames(df_long) <- c("Time", "Pixel", "Value")
    df_filtered <- subset(df_long, Pixel == target_s0) 
    df_filtered$Episode <- as.factor(i)
    df_list[[i]] <- df_filtered
  }
  
  df_all <- do.call(rbind, df_list)
  u_episode <- selected_points$u_s0[matching_indices[1]]
  p <- ggplot(df_all, aes(x = Time, y = Value, group = Episode)) +
    geom_line(color = btfgreen, alpha = 0.7) +
    labs(x = "Time", y = "Rainfall (mm)") +
    theme_minimal() +
    geom_hline(yintercept = u_episode, color = "red", linetype = "dashed")

  # Find the min and max values of the 'Value' for setting the position of the annotation
  y_min <- min(df_all$Value, na.rm = TRUE)
  y_max <- max(df_all$Value, na.rm = TRUE)
  
  # Set the annotation's y position dynamically based on the min value
  annotation_y <- y_min - 0.1 * (y_max - y_min)  # Adjust the multiplier for optimal positioning

  # Add vertical line at t0 and annotation for t0
  # p <- p + geom_vline(xintercept = delta, color = "#ff0000a6", linetype = "dashed") +
  #   annotate("text", x = delta + 0.5, y = annotation_y, label = expression(t[0]),
  #            color = "red", vjust = 0, size = 5)
  
  
  filename <- paste0(im_folder, "optim/comephore/excess_episodes/extreme_episodes_site_", target_s0, 
                     "_min", min_spatial_dist, "km", "_delta_", delta, ".png")
  ggsave(filename, plot = p, width = 20, height = 15, units = "cm")

}


length(list_episodes)
# Verif first episode
episode <- list_episodes[[1]]
# find col and row of the max value
max_value <- max(episode)
max_value_coords <- which(episode == max_value, arr.ind = TRUE)
s0_col <- colnames(episode)[max_value_coords[2]]
t0_row <- rownames(episode)[max_value_coords[1]]

# selected_points$s0[1]
# selected_points$t0[1]
# rownames(comephore)[selected_points$t0[1]]

# histogram of the number of episodes per s0
df_hist <- data.frame(table(selected_points$s0))
colnames(df_hist) <- c("s0", "count")
ggplot(df_hist, aes(x = count)) +
  geom_histogram(binwidth = 1, fill = btfgreen, color = "black") +
  btf_theme +
  xlab("Site") +
  ylab("Count")

s0_max_ep <- df_hist$s0[which.max(df_hist$count)]
beta <- 0
list_overlaps <- c()
for (s in unique(selected_points$s0)) {
  overlaps <- check_intervals_overlap(s, selected_points, delta, beta)
  list_overlaps <- c(list_overlaps, overlaps)
}
any(overlaps)


# get month and year for p236
# list_date <- lapply(1:length(list_p236), function(ep) {
#   list_date_t0 <- rownames(list_p236[[ep]])[1]
#   date_t0 <- as.Date(list_date_t0)
#   day <- format(date_t0, "%d")
#   month <- format(date_t0, "%m")
#   year <- format(date_t0, "%Y")
#   list(day = day, month = month, year = year)
# })

# list_date <- do.call(rbind, list_date)
# rownames(list_date) <- c(1:47)
# colnames(list_date) <- c("day", "month", "year")
# year = "2010"
# list_date[list_date[,3] == year, ]
# # 34, 47
list_episodes_unif_points <- get_extreme_episodes(selected_points, comephore,
                                      episode_size = episode_size, unif = TRUE)

list_episodes_unif <- list_episodes_unif_points$episodes
# list_episodes[[100]] # check the episode
# nrow(list_episodes[[1]]) # size of episode

# boxplot(list_episodes[[1]])
s0_list <- selected_points$s0
t0_list <- selected_points$t0
u_list <- selected_points$u_s0


# filename <- paste0(im_folder, "optim/comephore/timeseries_p99_centre_zoom_t0.png")
# png(filename)
# plot(df_comephore$p99, type = "l", col = btfgreen, xlim = c(59270, 59340),
#      xlab = "Time", ylab = "Rainfall (mm)", ylim = c(0, 120))
# abline(v = 59303, col = "#cc3126", lty = 1, lwd = 1)
# abline(v = 59303 + 12, col = "#be8ba5", lty = 2)
# abline(v= 59303 - 12, col = "#be8ba5", lty = 2)
# dev.off()
library(parallel)

tau_vect <- 0:10
tmax <- max(tau_vect)
sites_coords <- as.data.frame(sites_coords)
# Compute the lags and excesses for each conditional point
list_results <- mclapply(1:length(s0_list), function(i) {
  s0 <- s0_list[i]
  s0_coords <- sites_coords[s0, ]
  # t0 <- t0_list[i]
  episode <- list_episodes_unif[[i]]
  u <- u_list[i]
  ind_t0_ep <- 0 # index of t0 in the episode
  lags <- get_conditional_lag_vectors(sites_coords, s0_coords, ind_t0_ep,
                                tau_vect, latlon = FALSE)
  # lags$hnorm <- lags$hnorm / 1000 # convert to km
  excesses <- empirical_excesses(episode, q, lags, type = "rpareto",
                                t0 = ind_t0_ep)
  list(lags = lags, excesses = excesses)
}, mc.cores = detectCores() - 1)

list_lags <- lapply(list_results, `[[`, "lags")
list_excesses <- lapply(list_results, `[[`, "excesses")


# params <- c(beta1, beta2, alpha1, alpha2, 0.1, 0.1)
# chi <- theoretical_chi(params, excesses)
# list_lags[[1]]
# list_excesses[[1]]$kij

s0 <- s0_list[1]
s0_coords <- sites_coords[s0, ]
excesses <- list_excesses[[1]]
# sort excesses by hnorm from smallest to largest
excesses$kij
excesses[1:20, ]
df_lags <- list_lags[[1]]
head(df_lags)

# ADD WIND DATA ################################################################

list_episodes_points <- get_extreme_episodes(selected_points, comephore,
                                    episode_size = episode_size, unif = FALSE)

list_episodes <- list_episodes_points$episodes

wind_per_episode <- Map(compute_wind_episode, list_episodes, s0_list, u_list,
               MoreArgs = list(wind_df = wind_mtp, speed_time = 0))

wind_ep_df <- do.call(rbind, wind_per_episode)
head(wind_ep_df)

# show rows where DD_t0 is different than DD_excess
wind_ep_df[wind_ep_df$DD_t0 != wind_ep_df$DD_excess, ]

sum(is.na(wind_ep_df$DD_deg))
sum(is.na(wind_ep_df$DD_t0))
# sum(is.na(wind_ep_df$DD_excess))
sum(is.na(wind_ep_df$FF))

# wind_ep_df$cardDirt0 <- sapply(wind_ep_df$DD_t0, convert_to_cardinal)

wind_ep_df$cardDir <- factor(wind_ep_df$cardDir,
                       levels = c("N", "NE", "E", "SE", "S", "SW", "W", "NW"))


# Define wind speed categories and reorder levels to change legend order
wind_ep_df <- wind_ep_df %>%
  mutate(FF_interval = factor(cut(FF,
                            breaks = c(0, 2, 5, 10, Inf),
                            labels = c("<2", "2-5", "5-10", ">10"),
                            right = FALSE),
                    levels = c(">10", "5-10", "2-5", "<2")))  # Reverse order

# Define wind direction order
wind_ep_df$cardDir <- factor(wind_ep_df$cardDir,
                    levels = c("N", "NE", "E", "SE", "S", "SW", "W", "NW"))

# Custom color palette
custom_colors <- c("#0335258e", "#1a755a8e", "#5b99868e", "#98d6c48e")

# Plot wind rose
ggplot(wind_ep_df, aes(x = cardDir, fill = FF_interval)) +
  geom_bar(position = "stack", width = 1, color = "black") +
  coord_polar(start = 15 * pi / 8) +
  scale_fill_manual(values = custom_colors, name = "Force (m/s)") +
  labs(x = "Direction", y = "Count", title = "") +
  theme_minimal() +
  theme(axis.text.x = element_text(size = 12, face = "bold"),
        legend.title = element_text(size = 12, face = "bold"),
        legend.text = element_text(size = 10))

# Save plot
end_filename <- paste(q*1000, "q_", n_episodes, "ep_", min_spatial_dist, "km_",
                      tmax, "h_delta_", delta, ".png", sep = "")
filename <- paste(im_folder, "wind/datagouv/excess_episodes/wind_card_dir_", 
                  end_filename, sep = "")
ggsave(filename, width = 20, height = 15, units = "cm")


# Ensure cardDir is a factor with correct order
wind_ep_df$cardDir_t0 <- factor(wind_ep_df$cardDir_t0,
                        levels = c("N", "NE", "E", "SE", "S", "SW", "W", "NW"))

# Plot wind rose
ggplot(wind_ep_df, aes(x = cardDir_t0, fill = FF_interval)) +
  geom_bar(position = "stack", width = 1, color = "black") +
  coord_polar(start = 15 * pi / 8) +
  scale_fill_manual(values = custom_colors, name = "Force (m/s)") +
  labs(x = "Direction in t0", y = "Count", title = "") +
  theme_minimal() +
  theme(axis.text.x = element_text(size = 12, face = "bold"),
        legend.title = element_text(size = 12, face = "bold"),
        legend.text = element_text(size = 10))

# save the wind data
filename <- paste(im_folder, "wind/datagouv/excess_episodes/wind_card_dir_t0_", end_filename,
                  sep = "")
ggsave(filename, width = 20, height = 15, units = "cm")

# Plot wind rose for DD_t0
ggplot(wind_ep_df, aes(x = DD_t0)) +
  geom_bar(position = "stack", width = 3, color = "#5048489d", fill = btfgreen,
        alpha = 0.8) +
  coord_polar(start = 15.85 * pi / 8) +
  labs(x = "Direction in t0 (degree)", y = "Count", title = "") +
  theme_minimal() +
  theme(axis.text.x = element_text(size = 12, face = "bold"),
        legend.title = element_text(size = 12, face = "bold"),
        legend.text = element_text(size = 10))

# save the wind data
filename <- paste(im_folder, "wind/datagouv/excess_episodes/wind_dir_t0_", end_filename,
                    sep = "")
ggsave(filename, width = 20, height = 15, units = "cm")


# Ensure cardDir is a factor with correct order
wind_ep_df$cardDir_excess <- factor(wind_ep_df$cardDir_excess,
                        levels = c("N", "NE", "E", "SE", "S", "SW", "W", "NW"))

# Plot wind rose
ggplot(wind_ep_df, aes(x = cardDir_excess, fill = FF_interval)) +
  geom_bar(position = "stack", width = 1, color = "black") +
  coord_polar(start = 15 * pi / 8) +
  scale_fill_manual(values = custom_colors, name = "Force (m/s)") +
  labs(x = "Direction", y = "Count", title = "") +
  theme_minimal() +
  theme(axis.text.x = element_text(size = 12, face = "bold"),
        legend.title = element_text(size = 12, face = "bold"),
        legend.text = element_text(size = 10))
# save the wind data
filename <- paste(im_folder, "wind/datagouv/excess_episodes/wind_card_dir_excess_", end_filename,
                  sep = "")
ggsave(filename, width = 20, height = 15, units = "cm")


# Plot wind rose for DD_excess
ggplot(wind_ep_df, aes(x = DD_excess)) +
  geom_bar(position = "stack", width = 3, color = "#5048489d", fill = btfgreen,
        alpha = 0.8) +
  coord_polar(start = 15.85 * pi / 8) +
  labs(x = "Direction (degree)", y = "Count", title = "") +
  theme_minimal() +
  theme(axis.text.x = element_text(size = 12, face = "bold"),
        legend.title = element_text(size = 12, face = "bold"),
        legend.text = element_text(size = 10))
# save the wind data
# save the wind data
filename <- paste(im_folder, "wind/datagouv/excess_episodes/wind_dir_excess_", end_filename,
                    sep = "")
ggsave(filename, width = 20, height = 15, units = "cm")



# OPTIMIZATION #################################################################

list_episodes_points <- get_extreme_episodes(selected_points, comephore,
                                    episode_size = episode_size, unif = FALSE)

list_episodes <- list_episodes_points$episodes

# Compute the wind vector for each episode (-FF because it's the wind direction)
# sin and cos are shifted because 0 degree means North
wind_ep_df$vx <- -wind_ep_df$FF * sin(wind_ep_df$DD_t0 * pi / 180)
wind_ep_df$vy <- -wind_ep_df$FF * cos(wind_ep_df$DD_t0 * pi / 180)
# if values of vx or vy are really close to 0, set them to 0
wind_ep_df$vx[abs(wind_ep_df$vx) < 1e-08] <- 0
wind_ep_df$vy[abs(wind_ep_df$vy) < 1e-08] <- 0

mean_vx <- mean(wind_ep_df$vx, na.rm = TRUE)
mean_vy <- mean(wind_ep_df$vy, na.rm = TRUE)

head(wind_ep_df)

wind_df <- wind_ep_df[, c("DD_t0", "vx", "vy")]
colnames(wind_df) <- c("dir", "vx", "vy")
# which episode has NA wind values
ind_NA <- which(is.na(wind_df$vx))

if (any(ind_NA > 0)) {
  # remove these episodes
  wind_opt <- wind_df[-ind_NA, ]
  episodes_opt <- list_episodes[-ind_NA]
  lags_opt <- list_lags[-ind_NA]
  excesses_opt <- list_excesses[-ind_NA]
} else {
  wind_opt <- wind_df
  episodes_opt <- list_episodes
  lags_opt <- list_lags
  excesses_opt <- list_excesses
}
save.image("workspace.RData")

# wind_df_test <- wind_df
# wind_df_test$vx <- 0.1
# wind_df_test$vy <- 0.1
# load("workspace.RData")
init_param <- c(beta1, beta2, alpha1, alpha2, 1, 1)
# init_param <- c(beta1, beta2, alpha1, alpha2, mean_vx, mean_vy)
# init_param <- c(0.01, 0.2, 1.5, 1, 0.2, 0.1)
# q <- 1
result <- optim(par = init_param, fn = neg_ll_composite,
        list_lags = lags_opt, list_episodes = episodes_opt,
        list_excesses = excesses_opt, hmax = NA,
        wind_df = wind_df,
        latlon = FALSE,
        directional = TRUE,
        method = "L-BFGS-B",
        fixed_eta1 = TRUE,
        fixed_eta2 = TRUE,
        lower = c(1e-08, 1e-08, 1e-08, 1e-08, 1e-08, 1e-08),
        upper = c(10, 10, 1.999, 1.999, 10, 10),
        control = list(maxit = 10000,
                      trace = 1,
                      parscale = c(1, 1, 1, 1, 1, 1)),
        hessian = F)
result

nll_list <- c()
beta1_list <- seq(0.0005, 0.1, length.out = 50)
for (beta in beta1_list) {
  params <- c(beta, beta2, alpha1, alpha2, 1, 1)
  nll <- neg_ll_composite(params, list_lags = lags_opt,
                          list_episodes = episodes_opt,
                          list_excesses = excesses_opt,
                          hmax = NA,
                          wind_df = wind_df,
                          latlon = FALSE,
                          directional = TRUE)
  nll_list <- c(nll_list, nll)
}


min_beta1 <- beta1_list[which.min(nll_list)]
min_nll <- min(nll_list)

df_plot <- data.frame(beta1 = beta1_list, nll = nll_list)
# Plot
p <- ggplot(df_plot, aes(x = beta1, y = nll)) +
  geom_line(color = "#4b4b8a", linewidth = 1.2) +
  geom_point(aes(x = min_beta1, y = min_nll), color = "red", size = 3) +
  geom_vline(xintercept = min_beta1, linetype = "dashed", color = "red") +
  geom_label(aes(x = min_beta1, y = min_nll, 
                 label = sprintf("minimum in %.4f", min_beta1)),
             vjust = -1.2, hjust = 0, size = 4, fill = "white", color = "red") +
  labs(x = expression(beta[1]), y = "- Log-likelihood") +
  theme_minimal()

p

# Save the plot
end_filename <- paste(q * 1000, "q_", min_spatial_dist, "km_", tmax,
                      "h_delta_", delta, ".pdf", sep = "")
filename <- paste0(im_folder, "optim/comephore/neg_ll_beta1_", end_filename, sep = "")
ggsave(plot= p, filename = filename, width = 20, height = 15, units = "cm")


nll_list <- c()
beta2_list <- seq(0.1, 0.8, length.out = 50)
for (beta in beta2_list) {
  params <- c(beta1, beta, alpha1, alpha2, 1, 1)
  nll <- neg_ll_composite(params, list_lags = lags_opt,
                          list_episodes = episodes_opt,
                          list_excesses = excesses_opt,
                          hmax = NA,
                          wind_df = wind_df,
                          latlon = FALSE,
                          directional = TRUE)
  nll_list <- c(nll_list, nll)
}


min_beta2 <- beta2_list[which.min(nll_list)]
min_nll <- min(nll_list)

df_plot <- data.frame(beta2 = beta2_list, nll = nll_list)
# Plot
p <- ggplot(df_plot, aes(x = beta2, y = nll)) +
  geom_line(color = "#4b4b8a", linewidth = 1.2) +
  geom_point(aes(x = min_beta2, y = min_nll), color = "red", size = 3) +
  geom_vline(xintercept = min_beta2, linetype = "dashed", color = "red") +
  geom_label(aes(x = min_beta2, y = min_nll, 
                 label = sprintf("minimum in %.4f", min_beta2)),
             vjust = -1.2, hjust = 0, size = 4, fill = "white", color = "red") +
  labs(x = expression(beta[2]), y = "- Log-likelihood") +
  theme_minimal()

p

# Save the plot
end_filename <- paste(q * 1000, "q_", min_spatial_dist, "km_", tmax,
                      "h_delta_", delta, ".pdf", sep = "")
filename <- paste0(im_folder, "optim/comephore/neg_ll_beta2_", end_filename, sep = "")
ggsave(plot= p, filename = filename, width = 20, height = 15, units = "cm")




# Séquence des beta > 0 (échelle naturelle)
beta1_seq <- seq(1e-4, 0.2, length.out = 50)
beta2_seq <- seq(0.1, 1, length.out = 50)

# Matrice pour stocker les valeurs NLL
nll_grid <- matrix(NA, nrow = length(beta1_seq), ncol = length(beta2_seq))

# Paramètres fixes
alpha1 <- 1.5
alpha2 <- 0.8
eta1 <- 1
eta2 <- 1

for (i in seq_along(beta1_seq)) {
  for (j in seq_along(beta2_seq)) {
    params_test <- c(beta1_seq[i], alpha1, beta2_seq[j], alpha2, eta1, eta2)
    
    nll_val <- tryCatch({
      neg_ll_composite(
        params = params_test,
        list_episodes = list_episodes,
        list_excesses = list_excesses,
        list_lags = list_lags,
        wind_df = NA,
        hmax = hmax,
        latlon = TRUE,
        directional = TRUE,
        rpar = TRUE
      )
    }, error = function(e) NA)
    
    nll_grid[i, j] <- nll_val
  }
}



library(reshape2)
library(ggplot2)

df_plot <- melt(nll_grid)
df_plot$beta1 <- rep(beta1_seq, each = length(beta2_seq))
df_plot$beta2 <- rep(beta2_seq, times = length(beta1_seq))

ggplot(df_plot, aes(x = beta1, y = beta2, fill = value)) +
  geom_tile() +
  scale_fill_viridis_c(name = "NLL", na.value = "white") +
  labs(
    x = expression(beta[1]),
    y = expression(beta[2]),
    title = "Profil de -log vraisemblance"
  ) +
  theme_minimal()

















# Check the convergence
if (result$convergence != 0) {
  warning("The optimization did not converge.")
}

# Extract the results
df_result <- data.frame(quantile = q,
                        delta = delta,
                        minDist = min_spatial_dist,
                        beta1 =  round(result$par[1], 3),
                        beta2 = round(result$par[2], 3),
                        alpha1 = round(result$par[3], 3),
                        alpha2 = round(result$par[4], 3),
                        eta1 = round(result$par[5], 3),
                        eta2 = round(result$par[6], 3))

kable(df_result, format = "latex") %>%
  kable_styling(bootstrap_options = c("striped", "hover", "condensed",
  "responsive"), latex_options = "H")


# # fix parameters
# neg_ll_composite_fixpar <- function(beta1, beta2, alpha1, alpha2, eta1, eta2, list_lags,
#                     list_excesses, wind_df, hmax = NA, latlon = TRUE) {
#   params <- c(beta1, beta2, alpha1, alpha2, eta1, eta2)
#   nll_composite <- neg_ll_composite(params, list_lags,
#                                       list_excesses, wind_df,
#                                       hmax = hmax, latlon = latlon,
#                                       directional = TRUE)
#   return(nll_composite)
# }

# library(bbmle)
# result <- mle2(neg_ll_composite_fixpar,
#               start = list(beta1 = init_param[1],
#                            beta2 = init_param[2],
#                            alpha1 = init_param[3],
#                            alpha2 = init_param[4],
#                            eta1 = 1,
#                            eta2 = 1),
#               data = list(list_lags = lags_opt,
#                   list_excesses = excesses_opt, hmax = 7,
#                   wind_df = wind_opt),
#               method = "L-BFGS-B",
#               lower = c(1e-08, 1e-08, 1e-08, 1e-08, 1e-08, 1e-08),
#               upper = c(Inf, Inf, 1.999, 1.999, Inf, Inf),
#               control = list(maxit = 10000),
#               fixed = list(beta2 = init_param[2], alpha2 = init_param[4]))
# result


# Vector of test values for eta1
eta1_values <- seq(0.1, 2, length.out = 50)

# Store for negative log-likelihood values
neg_log_likelihoods <- numeric(length(eta1_values))

params <- init_param
# Loop over different values of eta1
for (i in seq_along(eta1_values)) {
  params[5] <- eta1_values[i]
  neg_ll <- neg_ll_composite(params, list_lags = lags_opt,
                              list_episodes = episodes_opt, latlon = TRUE,
                              directional = TRUE,
                              list_excesses = excesses_opt, hmax = 7,
                              wind_df = wind_opt)
  neg_log_likelihoods[i] <- neg_ll
}

# Create data frame for ggplot
df_plot <- data.frame(eta1 =eta1_values, neg_log_likelihood = neg_log_likelihoods)

# Generate the plot using ggplot2
p <- ggplot(df_plot, aes(x = eta1, y = neg_log_likelihood)) +
  geom_line(color = "#4b4b8a", linewidth = 1.2) +
  labs(x = parse(text = "eta[1]"), y = "- Log-likelihood") +
  theme_minimal()
p
# Save the plot
end_filename <- paste(q * 1000, "q_", min_spatial_dist, "km_", tmax,
                      "h_delta_", delta, ".pdf", sep = "")
filename <- paste0(im_folder, "optim/comephore/likelihood_varying_eta1_", end_filename)

ggsave(filename, plot = p, width = 7, height = 5, dpi = 300)


# 3D plot of log-likelihood surface
eta1_values <- seq(0.8, 1.8, length.out = 20)
eta2_values <- seq(0.8, 1.2, length.out = 30)

ll_values <- matrix(nrow = length(eta1_values), ncol = length(eta2_values))

for (i in seq_along(eta1_values)) {
  for (j in seq_along(eta2_values)) {
    params_test <- init_param
    params_test[5] <- eta1_values[i]
    params_test[6] <- eta2_values[j]
    ll_values[i, j] <- generain::neg_ll_composite(params_test,
                                  list_lags = lags_opt,
                                  list_episodes = episodes_opt,
                                  list_excesses = excesses_opt, hmax = 7,
                                  wind_df = wind_opt)
  }
}

library(plotly)

fig <- plot_ly(
  x = eta1_values, y = eta2_values, z = ll_values,
  type = "surface", colorscale = "Viridis"
)

fig <- fig %>%
  layout(
    title = "3D Log-Likelihood Surface",
    scene = list(
      xaxis = list(title = "eta1"),
      yaxis = list(title = "eta2"),
      zaxis = list(title = "- Log-likelihood")
    )
  )

fig



contour(eta1_values, eta2_values, ll_values, nlevels = 40, 
  xlab = expression(eta[1]), ylab = expression(eta[2]),
  main = "Log-likelihood contours for eta1 and eta2")

# Save the plot
end_filename <- paste(q * 1000, "q_", min_spatial_dist, "km_", tmax,
                      "h_delta_", delta, ".png", sep = "")
filename <- paste0(im_folder, 
            "optim/comephore/log_likelihood_contours_eta1_eta2_", 
            end_filename)
png(filename, width = 800, height = 600)
contour(eta1_values, eta2_values, ll_values, nlevels = 40, 
  xlab = expression(eta[1]), ylab = expression(eta[2]),
  main = "Log-likelihood contours for eta1 and eta2")
dev.off()


# 2D plot of log-likelihood surface
eta1_values <- seq(0.1, 1.5, length.out = 30)
beta2_values <- seq(0.1, 3, length.out = 30)

ll_values <- matrix(nrow = length(eta1_values), ncol = length(beta2_values))
for (i in seq_along(eta1_values)) {
  for (j in seq_along(beta2_values)) {
    params_test <- init_param
    params_test[5] <- eta1_values[i]
    params_test[6] <- beta2_values[j]
    ll_values[i, j] <- generain::neg_ll_composite(params_test,
                                  list_lags = lags_opt,
                                  list_episodes = episodes_opt,
                                  list_excesses = excesses_opt, hmax = 7,
                                  wind_df = wind_opt)
  }
}

# Save the plot
end_filename <- paste(q * 1000, "q_", min_spatial_dist, "km_", tmax,
                      "h_delta_", delta, ".png", sep = "")
filename <- paste0(im_folder,
            "optim/comephore/log_likelihood_contours_eta1_beta2_", 
            end_filename)
png(filename, width = 800, height = 600)
contour(eta1_values, beta2_values, ll_values, nlevels = 40, 
  xlab = expression(eta[1]), ylab = expression(beta[2]),
  main = "Log-likelihood contours for eta1 and beta2")
dev.off()

# 3D plot of log-likelihood surface
fig <- plot_ly(
  x = eta1_values, y = beta2_values, z = ll_values,
  type = "surface", colorscale = "Viridis"
)

fig <- fig %>%
  layout(
    title = "3D Log-Likelihood Surface",
    scene = list(
      xaxis = list(title = "eta1"),
      yaxis = list(title = "beta2"),
      zaxis = list(title = "- Log-likelihood")
    )
  )
fig


htmlwidgets::saveWidget(fig, "interactive_plot.html")
# Save the figure as a PNG
orca(fig, "plot.png")

# GENERATE TABLE FROM TEST #####################################################

library(knitr)
library(kableExtra)
library(parallel)

# Configurations
q_values <- c(0.998)
delta_values <- c(12)
min_spatial_dist_values <- c(5)
tau_max <- c(10)
wind_direction_type <- c("t_0")

# Create a parameter grid
param_grid <- expand.grid(
  q = q_values,
  delta = delta_values,
  min_spatial_dist = min_spatial_dist_values,
  tau_max = tau_max,
  wind_dir = wind_direction_type,
  stringsAsFactors = FALSE
)

# Function to process each parameter set
process_params <- function(sites_coords, params) {
  q <- params$q
  delta <- params$delta
  episode_size  <- 2 * delta
  min_spatial_dist <- params$min_spatial_dist
  tau_max <- params$tau_max
  wind_dir <- params$wind_dir

  print(paste("q:", q, "delta:", delta, "minDist:", min_spatial_dist,
              "taumax:", tau_max, "direction:", wind_dir))

  # Select (s0, t0) pairs for extreme episodes
  selected_points <- select_extreme_episodes(sites_coords, comephore, q,
                                             min_spatial_dist,
                                             episode_size = episode_size,
                                             n_max_episodes = 10000,
                                             time_ext = 0)
  n_episodes <- length(selected_points$s0)

  # Get the extreme episodes
  list_episodes_unif_points <- get_extreme_episodes(selected_points,
                                      comephore,
                                      episode_size = episode_size, unif = TRUE)
  list_episodes_unif <- list_episodes_unif_points$episodes

  list_episodes_points <- get_extreme_episodes(selected_points, comephore,
                                    episode_size = episode_size, unif = FALSE)

  list_episodes <- list_episodes_points$episodes

  s0_list <- selected_points$s0
  u_list <- selected_points$u_s0
  tau_vect <- 0:tau_max
  ind_t0_ep <- 0  # index of t0 in the episode

  print("Compute the lags and excesses for each conditional point")

  # Compute lags and excesses in parallel
  list_results <- mclapply(1:length(s0_list), function(i) {
    s0 <- s0_list[i]
    s0_coords <- sites_coords[s0, ]
    episode <- list_episodes_unif[[i]]

    lags <- get_conditional_lag_vectors(sites_coords, s0_coords,
                                ind_t0_ep, tau_vect, latlon = TRUE)
    excesses <- empirical_excesses(episode, q, lags, type = "rpareto",
                                    t0 = ind_t0_ep)
    list(lags = lags, excesses = excesses)
  }, mc.cores = detectCores() - 1)

  list_lags <- lapply(list_results, `[[`, "lags")
  list_excesses <- lapply(list_results, `[[`, "excesses")
  
  # Get wind components for each episode
  wind_per_episode <- Map(compute_wind_episode, list_episodes, s0_list, u_list,
                          MoreArgs = list(wind_df = wind_mtp, speed_time = 0))

  wind_ep_df <- do.call(rbind, wind_per_episode)

  if (wind_dir == "t_0") { # Use wind direction at t0
    wind_ep_df$vx <- -wind_ep_df$FF * sin(wind_ep_df$DD_t0 * pi / 180)
    wind_ep_df$vy <- -wind_ep_df$FF * cos(wind_ep_df$DD_t0 * pi / 180)
  } else { # Use excess wind direction
    wind_ep_df$vx <- -wind_ep_df$FF * sin(wind_ep_df$DD_excess * pi / 180)
    wind_ep_df$vy <- -wind_ep_df$FF * cos(wind_ep_df$DD_excess * pi / 180)
  }

  wind_df <- wind_ep_df[, c("vx", "vy")]

  # Remove NA values
  ind_NA <- which(is.na(wind_df$vx))
  if (length(ind_NA) > 0) {
    wind_df <- wind_df[-ind_NA, ]
    list_lags <- list_lags[-ind_NA]
    list_excesses <- list_excesses[-ind_NA]
  }

  # Initial parameters
  init_param <- c(beta1, beta2, alpha1, alpha2, 1, 1)

  print("Optimization")

  # Optimization
  result <- optim(par = init_param, fn = neg_ll_composite,
                  list_episodes = list_episodes,
                  list_lags = list_lags,
                  list_excesses = list_excesses, hmax = 7,
                  wind_df = wind_df, latlon = TRUE, directional = TRUE,
                  method = "L-BFGS-B",
                  lower = c(1e-08, 1e-08, 1e-08, 1e-08, 1e-08, 1e-08),
                  upper = c(Inf, Inf, 1.999, 1.999, Inf, Inf),
                  control = list(maxit = 10000, trace = 1),
                  hessian = FALSE)

  # Store results
  if (result$convergence != 0) {
    df_result <- data.frame(q, delta, min_spatial_dist, tau = tau_max,
                            wind_dir, beta1 = "NC", beta2 = "NC", alpha1 = "NC",
                            alpha2 = "NC", eta1 = "NC", eta2 = "NC")
  } else {
    df_result <- data.frame(q, delta, min_spatial_dist, tau = tau_max,
                            wind_dir,
                            beta1 = round(result$par[1], 3),
                            beta2 = round(result$par[2], 3),
                            alpha1 = round(result$par[3], 3),
                            alpha2 = round(result$par[4], 3),
                            eta1 = round(result$par[5], 3),
                            eta2 = round(result$par[6], 3))
  }

  df_result$nepisodes <- n_episodes
  return(df_result)
}

# Run in parallel
results_list <- mclapply(seq_len(nrow(param_grid)), function(i) {
  process_params(sites_coords, param_grid[i, ])
}, mc.cores = detectCores() - 1)

# Combine results
final_results <- do.call(rbind, results_list)
# remove tau column
final_results <- final_results[, -which(names(final_results) == "tau")]
# put nepisodes column before beta1 column
final_results <- final_results[, c(1, 2, 3, 4, 11, 5, 6, 7, 8, 9, 10)]

# Create LaTeX table
kable(final_results, format = "latex",
      caption = "Optimization results with wind components as covariates") %>%
  kable_styling(bootstrap_options = c("striped", "hover", "condensed", "responsive"), latex_options = "H")


# final_result_tau10 <- final_results

# ADVECTION ESTIMATION #########################################################

# get first result
df_result <- final_results[7, ]
df_result

# Estimated advection vector with wind components
adv_hat <- abs(wind_df) ^ as.numeric(df_result$eta1) * sign(wind_df) *
                            as.numeric(df_result$eta2)
print(adv_hat)

# remove NA values
adv_hat <- adv_hat[!is.na(adv_hat$vx), ]

# get wind FF and DD from adv_hat
adv_hat$FF_hat_kmh <- sqrt(adv_hat$vx^2 + adv_hat$vy^2) # Wind magnitude in km/h
adv_hat$FF_hat_ms <- adv_hat$FF_hat_kmh / 3.6  # Convert to m/s
adv_hat$FF_hat_mms <- adv_hat$FF_hat_ms * 1000  # Convert to mm/s
adv_hat$DD_hat <- atan2(adv_hat$vy, adv_hat$vx) * 180 / pi
head(adv_hat)

adv_hat$DD_hat <- ifelse(adv_hat$DD_hat < 0, adv_hat$DD_hat + 360,
                        adv_hat$DD_hat)
head(adv_hat)

adv_hat$cardDir <- sapply(adv_hat$DD_hat, convert_to_cardinal)
head(adv_hat)

# Define wind speed categories and reorder levels to change legend order
wind_plot_df <- adv_hat %>%
  mutate(FF_interval = factor(cut(FF_hat_kmh,
                            breaks = c(0, 1, 1.2, 1.4, 1.6, 1.8, 2, Inf),
                            labels = c("<1", "1-1.2", "1.2-1.4", "1.4-1.6",
                                       "1.6-1.8", "1.8-2", ">2"),
                            right = FALSE),
                            levels = c(">2", "1.8-2", "1.6-1.8",
                                       "1.4-1.6", "1.2-1.4", "1-1.2", "<1")))  # Reverse order

# Define wind direction order
wind_plot_df$cardDir <- factor(wind_plot_df$cardDir,
                    levels = c("N", "NE", "E", "SE", "S", "SW", "W", "NW"))

# remove NA values
wind_plot_df <- wind_plot_df[!is.na(wind_plot_df$cardDir), ]
# Custom color palette
custom_colors <- c("#0335258e", "#1a755a8e", "#5b99868e", "#98d6c48e",
                   "#d9e4e8", "#f0f0f0", "#f7f7f7", "#ffffff")

# Plot wind rose
advplot <- ggplot(wind_plot_df, aes(x = cardDir, fill = FF_interval)) +
  geom_bar(position = "stack", width = 1, color = "black") +
  coord_polar(start = 15 * pi / 8) +
  scale_fill_manual(values = custom_colors, name = "Force (km/h)") +
  labs(x = "Direction", y = "Count", title = "") +
  theme_minimal() +
  theme(axis.text.x = element_text(size = 12, face = "bold"),
        legend.title = element_text(size = 12, face = "bold"),
        legend.text = element_text(size = 10))

file_path <- paste(im_folder, "optim/comephore/excess_episodes/estimated_adv_",
                    df_result$q*1000, "_", df_result$delta, "_",
                    df_result$min_spatial_dist, "_", df_result$tau,
                    "_", df_result$wind_dir, ".png", sep = "")
ggsave(file_path, plot = advplot, width = 20, height = 15, units = "cm")



# VARIOGRAM PLOTS ##############################################################

# compute variogram with parameters
result <- df_result
# convert last column into numeric
result[, c("beta1", "beta2", "alpha1", "alpha2", "eta1", "eta2")] <-
  lapply(result[, c("beta1", "beta2", "alpha1", "alpha2", "eta1", "eta2")],
         as.numeric)

beta1 <- result$beta1
beta2 <- result$beta2
alpha1 <- result$alpha1
alpha2 <- result$alpha2
eta1 <- result$eta1
eta2 <- result$eta2

# Compute advection for each row
adv_df <- wind_df %>%
  mutate(
    adv_x = (abs(vx)^eta1) * sign(vx) * eta2,
    adv_y = (abs(vy)^eta1) * sign(vy) * eta2
  )
  
# remove NA values
adv_df <- adv_df[!is.na(adv_df$adv_x), ]
head(adv_df)
list_lags <- lapply(seq_along(list_lags), function(i) {
  list_lags[[i]] %>%
    mutate(vx = wind_df$vx[i], vy = wind_df$vy[i])
})

params <- c(beta1, beta2, alpha1, alpha2, adv_df$adv_x, adv_df$adv_y)

list_chi <- lapply(seq_along(list_lags), function(i) {
  theoretical_chi(params, list_lags[[i]], latlon = TRUE)
})

head(list_chi[[1]])

df_chi_all <- bind_rows(list_chi, .id = "episode") 
head(df_chi_all)


library(ggplot2)

df_chi_all$hlagabs <- abs(df_chi_all$hlag)
ggplot(df_chi_all, aes(x = hlagabs, y = vario)) +
  geom_line() +
  facet_wrap(~ tau, scales = "free_x",
      labeller = labeller(tau = function(x) paste0("tau = ", x))) +
  labs(
    title = "Directional variogram",
    x = "Spatial lag",
    y = "Variogram"
  ) +
  scale_color_discrete(name = expression(tau), breaks = sort(unique(df_chi_all$tau))) +
  theme_minimal()

n_episodes <- length(unique(df_chi_all$episode))
# save the plot
end_filename <- paste(df_result$q*1000, "_", df_result$min_spatial_dist, "km_",
                     df_result$wind_dir, "_delta", delta, ".png", sep = "")
filename <- paste(im_folder, "optim/comephore/vario_estim/variogram_rpareto_all",
                    "_", end_filename, sep = "")
ggsave(filename, width = 20, height = 15, units = "cm")

library(ggplot2)


# Order tau (ascending or your custom order)
df_chi_all$tau <- factor(df_chi_all$tau, levels = rev(sort(unique(df_chi_all$tau))))

# Now plot
ggplot(df_chi_all, aes(x = hlagabs, y = vario, color = tau)) +
  geom_line() +
  labs(
    title = "Directional variogram",
    x = "Spatial lag",
    y = "Variogram",
    color = expression(tau)
  ) +
  theme_minimal()

# save the plot
end_filename <- paste(df_result$q*1000, "_", df_result$min_spatial_dist, "km_",
                     df_result$wind_dir, "_delta", delta, ".png", sep = "")
filename <- paste(im_folder, "optim/comephore/vario_estim/variogram_rpareto_all",
                    "_", end_filename, sep = "")
ggsave(filename, width = 20, height = 15, units = "cm")


# Plot only for tau = 10
df_chi_tau10 <- df_chi_all[df_chi_all$tau == 10, ]

ggplot(df_chi_tau10, aes(x = hlagabs, y = vario)) +
  geom_line() +
  labs(
    title = "",
    x = "Spatial lag",
    y = "Variogram"
  ) +
  theme_minimal()

# save the plot
filename <- paste(im_folder, "optim/comephore/vario_estim/variogram_rpareto_tau10",
                    "_", end_filename, sep = "")
ggsave(filename, width = 20, height = 15, units = "cm")


# Plot only for tau = 10
df_chi_tau10 <- df_chi_all[df_chi_all$tau == 0, ]

ggplot(df_chi_tau10, aes(x = hlagabs, y = vario)) +
  geom_line() +
  labs(
    title = "",
    x = "Spatial lag",
    y = "Variogram"
  ) +
  theme_minimal()

# save the plot
filename <- paste(im_folder, "optim/comephore/vario_estim/variogram_rpareto_tau0",
                    "_", end_filename, sep = "")
ggsave(filename, width = 20, height = 15, units = "cm")


# Revoir la pres du variogram: V=0 et une avec V qui vient du sud est et comparaison