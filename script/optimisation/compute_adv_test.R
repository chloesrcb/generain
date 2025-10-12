library(dplyr)
library(geosphere)  # pour distances en km
library(tidyr)

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
comephore <- df_comephore[-1] # remove dates column

filename_era5 <- paste0(data_folder, "wind/ERA5/wind_era5.csv")
wind_era5 <- read.csv(filename_era5)
head(wind_era5)

wind_era5$time <- ifelse(nchar(wind_era5$time) == 10,
                            paste0(wind_era5$time, " 00:00:00"),
                            wind_era5$time)

wind_era5$time <- as.POSIXct(wind_era5$time,
                             format = "%Y-%m-%d %H:%M:%S", tz = "UTC")

filename_wind <- paste0(data_folder, "wind/data_gouv/wind_mtp.csv")
wind_mtp_h <- read.csv(filename_wind)
head(wind_mtp_h)

wind_mtp_h$datetime <- ifelse(nchar(wind_mtp_h$datetime) == 10,
                            paste0(wind_mtp_h$datetime, " 00:00:00"),
                            wind_mtp_h$datetime)
# Convert datetime to POSIXct
wind_mtp_h$datetime <- as.POSIXct(wind_mtp_h$datetime,
                                format = "%Y-%m-%d %H:%M:%S", tz = "UTC")
head(wind_mtp_h$datetime)


library(dplyr)
library(tidyr)
library(geosphere)

# COMEPHORE en format long : date, pixel_name, rain
come_long <- df_comephore %>%
  pivot_longer(-date, names_to = "pixel_name", values_to = "rain") %>%
  left_join(loc_px, by = "pixel_name")

head(come_long)
bary <- come_long %>%
  group_by(date) %>%
  summarise(
    rain_tot = sum(rain, na.rm = TRUE),
    lon_c = if (rain_tot > 0.5) weighted.mean(Longitude, rain, na.rm = TRUE) else NA_real_,
    lat_c = if (rain_tot > 0.5) weighted.mean(Latitude, rain, na.rm = TRUE) else NA_real_
  ) %>%
  arrange(date)


bary <- bary %>%
  mutate(
    lon_next = lead(lon_c),
    lat_next = lead(lat_c),
    dt_h = as.numeric(difftime(lead(date), date, units = "hours")),
    dist_km = distVincentyEllipsoid(cbind(lon_c, lat_c),
                                    cbind(lon_next, lat_next)) / 1000,
    speed_prec_kmh = dist_km / dt_h,
    bearing = bearing(cbind(lon_c, lat_c), cbind(lon_next, lat_next)),
    vx_prec = speed_prec_kmh * sin(bearing*pi/180),
    vy_prec = speed_prec_kmh * cos(bearing*pi/180)
  )

# convertir ERA5 en km/h
wind_era5 <- wind_era5 %>%
  mutate(u_kmh = u*3.6, v_kmh = v*3.6)

merged <- left_join(bary, wind_era5, by = c("date" = "time")) %>%
  mutate(
    speed_wind = sqrt(u_kmh^2 + v_kmh^2),
    dot = vx_prec*u_kmh + vy_prec*v_kmh,
    angle_deg = acos(pmin(1, pmax(-1, dot/(speed_prec_kmh*speed_wind))))*180/pi
  )



library(dplyr)
library(geosphere)

# 1️⃣ Calcul de l'advection réelle (barycentre)
bary <- come_long %>%
  group_by(date) %>%
  summarise(
    rain_tot = sum(rain, na.rm = TRUE),
    lon_c = if (rain_tot > 10) weighted.mean(Longitude, rain, na.rm = TRUE) else NA_real_,
    lat_c = if (rain_tot > 10) weighted.mean(Latitude, rain, na.rm = TRUE) else NA_real_
  ) %>%
  arrange(date) %>%
  mutate(
    lon_next = lead(lon_c),
    lat_next = lead(lat_c),
    dt_h = as.numeric(difftime(lead(date), date, units = "hours")),
    dist_km = distVincentyEllipsoid(cbind(lon_c, lat_c),
                                    cbind(lon_next, lat_next)) / 1000,
    speed_prec_kmh = dist_km / dt_h,
    bearing = bearing(cbind(lon_c, lat_c), cbind(lon_next, lat_next)),
    vx_prec = speed_prec_kmh * sin(bearing*pi/180),
    vy_prec = speed_prec_kmh * cos(bearing*pi/180)
  )
head(bary)


# 2️⃣ Préparer les vents ERA5 en km/h
wind_era5 <- wind_era5 %>%
  mutate(
    u_kmh = u*3.6,
    v_kmh = v*3.6
  )

wind_mtp_h <- wind_mtp_h %>%
    mutate(
        u_mf = FF * sin((90 - DD) * pi / 180) * 3.6,
        v_mf = FF * cos((90 - DD) * pi / 180) * 3.6
    )

# rename datetime to date
colnames(wind_mtp_h)[colnames(wind_mtp_h) == "datetime"] <- "date"
# 3️⃣ Ajouter le vent Météo France si disponible (10 m)
# Supposons meteoFrance contient: date, u_kmh, v_kmh
merged <- bary %>%
  left_join(wind_era5, by = c("date" = "time")) %>%
  left_join(wind_mtp_h, by = "date") %>%
  mutate(
    speed_wind = sqrt(u_kmh^2 + v_kmh^2),
    speed_wind_mf = sqrt(u_mf^2 + v_mf^2),
    dot_era = vx_prec*u_kmh + vy_prec*v_kmh,
    dot_mf = vx_prec*u_mf + vy_prec*v_mf,
    angle_mf = acos(pmin(1, pmax(-1, dot_mf/(speed_prec_kmh*speed_wind_mf))))*180/pi,
    angle_era = acos(pmin(1, pmax(-1, dot_era/(speed_prec_kmh*speed_wind))))*180/pi,
  )


# keep only rain_tot > 10
merged <- merged[merged$rain_tot > 100, ]
nrow(merged)

# 4️⃣ Modèle linéaire multivarié (vx et vy séparés)
lm_vx <- lm(vx_prec ~ u_kmh + v_kmh + u_mf + v_mf, data = merged)
lm_vy <- lm(vy_prec ~ u_kmh + v_kmh + u_mf + v_mf, data = merged)

summary(lm_vx)
summary(lm_vy)

# 5️⃣ Prédiction et évaluation
merged <- merged %>%
  mutate(
    vx_pred = predict(lm_vx, newdata = merged),
    vy_pred = predict(lm_vy, newdata = merged),
    speed_pred = sqrt(vx_pred^2 + vy_pred^2),
    angle_pred = atan2(vy_pred, vx_pred)*180/pi
  )

# Comparaison : advection réelle vs prédite
plot(merged$speed_prec_kmh, merged$speed_pred,
     xlab="vitesse barycentre (km/h)", ylab="vitesse prédite (km/h)")
abline(0,1, col="red")
