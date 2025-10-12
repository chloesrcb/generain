
filename_rain <- paste0(data_folder, "omsev/omsev_5min/rain_mtp_5min_2019_2024_cleaned.csv")
rain_omsev <- read.csv(filename_rain, sep = ",")

filename_wind <- paste0(data_folder, "wind/data_gouv/wind_mtp.csv")
wind_mtp_h <- read.csv(filename_wind)
head(wind_mtp_h)

filename_era5 <- paste0(data_folder, "wind/ERA5/wind_era5.csv")
wind_era5 <- read.csv(filename_era5)
head(wind_era5)

wind_era5$time <- ifelse(nchar(wind_era5$time) == 10,
                            paste0(wind_era5$time, " 00:00:00"),
                            wind_era5$time)

wind_era5$time <- as.POSIXct(wind_era5$time,
                             format = "%Y-%m-%d %H:%M:%S", tz = "UTC")

head(wind_era5$time)
# When datetime do not have hours and minutes, add 00:00:00
wind_mtp_h$datetime <- ifelse(nchar(wind_mtp_h$datetime) == 10,
                            paste0(wind_mtp_h$datetime, " 00:00:00"),
                            wind_mtp_h$datetime)
# Convert datetime to POSIXct
wind_mtp_h$datetime <- as.POSIXct(wind_mtp_h$datetime,
                                format = "%Y-%m-%d %H:%M:%S", tz = "UTC")
head(wind_mtp_h$datetime)

# get only data between 2019 and 2025
wind_mtp_h <- wind_mtp_h[year(wind_mtp_h$datetime) >= 2019 & year(wind_mtp_h$datetime) <= 2025, ]
# Apply function to the wind_direction column
wind_mtp_h$cardDir <- sapply(wind_mtp_h$DD, convert_to_cardinal)
wind_mtp_h$cardDir <- as.character(wind_mtp_h$cardDir)  # Ensure it's character
wind_mtp_h$cardDir[is.na(wind_mtp_h$DD)] <- NA
wind_mtp_h <- wind_mtp_h[, c("datetime", "FF", "DD", "cardDir")]
colnames(wind_mtp_h) <- c("datetime", "speed", "direction", "cardDir")

q <- 0.95
delta <- 30
min_spatial_dist <- 7
# get omsev advection for each episode
folder_adv <- paste0(data_folder, "comephore/adv_estim/")
list_files <- list.files(folder_adv)
adv_filename <- paste(data_folder, "/comephore/adv_estim/advection_results_q",
                      q * 100, "_delta", delta, "_dmin", min_spatial_dist,
                      ".csv", sep = "")
adv_df <- read.csv(adv_filename, sep = ",")
head(adv_df)
# number of 0 advection
nrow(adv_df[adv_df$mean_dx_kmh == 0 & adv_df$mean_dy_kmh == 0, ])
nrow(adv_df)
adv_df$duration_hours
# if durations is more than 48 hours, put nan in adv
adv_df$mean_dx_kmh[adv_df$duration_hours > 24] <- NA
adv_df$mean_dy_kmh[adv_df$duration_hours > 24] <- NA

# count na 
nrow(adv_df[is.na(adv_df$mean_dx_kmh) | is.na(adv_df$mean_dy_kmh), ])

# Transform t0 to POSIXct
adv_df$t0 <- as.POSIXct(adv_df$t0, format = "%Y-%m-%d %H:%M:%S", tz = "UTC")
wind_mtp_h$datetime <- as.POSIXct(wind_mtp_h$datetime, format = "%Y-%m-%d %H:%M:%S", tz = "UTC")
head(adv_df)
head(wind_era5)

adv_df_hour <- adv_df
# round t0 to the nearest hour
adv_df_hour$t0 <- round_date(adv_df_hour$t0, unit = "hour")
head(adv_df_hour)

merged <- adv_df_hour %>%
  left_join(wind_mtp_h, by = c("t0" = "datetime"))

tail(merged)

# remove rows with NA in speed or direction
merged <- merged[!is.na(merged$speed) & !is.na(merged$direction), ]
head(merged)
# compute wind vectors
merged$wind_x <- merged$speed * sin((90 - merged$direction) * pi / 180)
merged$wind_y <- merged$speed * cos((90 - merged$direction) * pi / 180)

# from m/s to km/h
merged$wind_x <- merged$wind_x * 3.6
merged$wind_y <- merged$wind_y * 3.6

# keep only useful columns
merged <- merged[, c("t0", "mean_dx_kmh", "mean_dy_kmh", "wind_x", "wind_y", "speed", "direction")]
head(merged)


# keep only non nul mean_dx_kmh and mean_dy_kmh
merged <- merged[!(merged$mean_dx_kmh == 0 & merged$mean_dy_kmh == 0), ]
nrow(merged)
# plot adv vs wind
p1 <- ggplot(merged, aes(x = mean_dx_kmh, y = wind_x)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE, color = "blue") +
  labs(title = "Advection in x vs Wind in x",
       x = "Advection in x (km/h)",
       y = "Wind in x (km/h)") +
  theme_minimal()

p2 <- ggplot(merged, aes(x = mean_dy_kmh, y = wind_y)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE, color = "blue") +
  labs(title = "Advection in y vs Wind in y",
       x = "Advection in y (km/h)",
       y = "Wind in y (km/h)") +
  theme_minimal()

# print plots
print(p1)
print(p2)   

# # check correlation
# cor(merged$mean_dx_kmh, merged$wind_y, use="complete.obs")
# cor(merged$mean_dy_kmh, merged$wind_y, use="complete.obs")
# lm(mean_dx_kmh ~ wind_y, data = merged)

# vitesses (normes)
merged$adv_speed <- sqrt(merged$mean_dx_kmh^2 + merged$mean_dy_kmh^2)

#plot adv speed vs wind speed
p3 <- ggplot(merged, aes(x = adv_speed, y = speed)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE, color = "blue") +
  labs(title = "Advection speed vs Wind speed",
       x = "Advection speed (km/h)",
       y = "Wind speed (km/h)") +
  theme_minimal()

print(p3)
colnames(merged)
# directions
merged$adv_angle <- atan2(merged$mean_dx_kmh, merged$mean_dy_kmh) * 180/pi
merged$wind_angle <- atan2(merged$wind_x, merged$wind_y) * 180/pi

# plot adv angle vs wind angle
p4 <- ggplot(merged, aes(x = adv_angle, y = direction)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE, color = "blue") +
  labs(title = "Advection angle vs Wind angle",
       x = "Advection angle (degrees)",
       y = "Wind angle (degrees)") +
  theme_minimal()

print(p4)
library(ggplot2)
library(circular)

# --- 1. Calculer les directions en degrés --------------------------
# atan2(x, y) -> angle en radians par rapport au nord
merged$adv_angle <- atan2(merged$mean_dx_kmh, merged$mean_dy_kmh) * 180/pi
merged$wind_angle <- atan2(merged$wind_x, merged$wind_y) * 180/pi

# ramener entre 0 et 360
merged$adv_angle <- (merged$adv_angle + 360) %% 360
merged$wind_angle <- (merged$wind_angle + 360) %% 360

# --- 2. Histogrammes de directions --------------------------------
p1 <- ggplot(merged, aes(x = adv_angle)) +
  geom_histogram(binwidth = 15, fill = "steelblue", color = "white") +
  coord_polar(start = -pi/2) +
  labs(title = "Distribution des directions d'advection",
       x = "Direction (°)", y = "Fréquence") +
  theme_minimal()

p2 <- ggplot(merged, aes(x = wind_angle)) +
  geom_histogram(binwidth = 15, fill = "darkred", color = "white") +
  coord_polar(start = -pi/2) +
  labs(title = "Distribution des directions de vent (10 m)",
       x = "Direction (°)", y = "Fréquence") +
  theme_minimal()

print(p1)
print(p2)

# --- 3. Corrélation circulaire -----------------------------------
adv_circ  <- circular(merged$adv_angle, units = "degrees", template = "geographics")
wind_circ <- circular(merged$wind_angle, units = "degrees", template = "geographics")

cor_circ <- cor.circular(adv_circ, wind_circ)
print(cor_circ)



head(rain_omsev)
rain_omsev$dates <- as.POSIXct(rain_omsev$dates, format="%Y-%m-%d %H:%M:%S", tz="UTC")
rain_omsev$datetime_rounded <- ceiling_date(rain_omsev$dates, unit = "hour")
# correlation between rain and wind speed
wind_wind <- wind_era5 %>%
  left_join(wind_mtp_h,
            by = c("dates" = "datetime"))

sites <- c("archie","archiw","cefe","chu1","chu2","chu3","chu4","chu5","chu6","chu7")

for (s in sites) {
  cat("Site:", s, "\n")
  print(cor(rain_wind[[s]], rain_wind$u, use="complete.obs"))
}

wind_wind$u_kmh_mtp <- wind_wind$speed.x * sin((90 - wind_wind$direction.x) * pi / 180)
wind_wind$v_kmh_mtp <- wind_wind$speed.x * cos((90 - wind_wind$direction.x) * pi / 180)

# plot rain vs wind speed
p1 <- ggplot(wind_wind, aes(x = u_kmh_mtp, y = speed.y)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE, color = "blue") +
  labs(title = "Rain vs Wind speed",
       x = "Wind speed (km/h)",
       y = "Rain (mm/5min)") +
  theme_minimal() 

print(p1)

filename_com <- paste0(data_folder, "comephore/comephore_2008_2024_5km.csv")
# filename_com <- paste0(data_folder, "comephore/zoom_5km.csv")
comephore_raw <- read.csv(filename_com, sep = ",")
# filename_loc <- paste0(data_folder, "comephore/coords_pixels_10km.csv")
filename_loc <- paste0(data_folder, "comephore/loc_px_zoom_5km.csv")
loc_px <- read.csv(filename_loc, sep = ",")

# remove pixel in loc_px that are not in comephore_raw
loc_px <- loc_px[loc_px$pixel_name %in% colnames(comephore_raw)[-1], ]
# reindex
rownames(loc_px) <- NULL
df_comephore <- as.data.frame(comephore_raw)
df_comephore$date <- as.POSIXct(df_comephore$date,
                                format = "%Y-%m-%d %H:%M:%S", tz = "UTC")

# Take only data after 2007
df_comephore <- df_comephore[df_comephore$date >= "2008-01-01", ]
rownames(df_comephore) <- format(as.POSIXct(df_comephore$date), "%Y-%m-%d %H:%M:%S")
rain_com <- df_comephore

rain_wind <- rain_com %>%
  left_join(adv_df,
            by = c("date" = "t0"))

rain_wind$speed <- sqrt(rain_wind$mean_dx_kmh^2 + rain_wind$mean_dy_kmh^2)
rain_wind$mean_dx_kmh
# subset to dates after 2019-01-01
rain_wind <- rain_wind[rain_wind$date >= "2019-01-01", ]
rain_wind$speed
# plot rain at pixel 18 vs wind speed with lag of 15 min
ggplot(rain_wind, aes(x = speed, y = p18)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE, color = "blue") +
  labs(title = "Rain at pixel 18 vs Wind speed",
       x = "Wind speed (km/h)",
       y = "Rain (mm/5min)") +
  theme_minimal()


library(dplyr)
rain_wind <- rain_wind %>%
  arrange(date) %>%
  mutate(wind_speed_lag15 = lag(speed, 3))  # si données 5 min, lag 1 = 5 min
cor(rain_wind$p18, rain_wind$wind_speed_lag15, use="complete.obs")
