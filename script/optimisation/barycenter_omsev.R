# get rain data from omsev
filename_omsev <- paste0(data_folder,
                         "omsev/omsev_5min/rain_mtp_5min_2019_2024.RData")

load(filename_omsev)


rain_omsev <- as.data.frame(rain.all5[, c(1, 6:(ncol(rain.all5) - 1))])
tail(rain_omsev$dates)
# remove rain before september 2019 and after january 2024
# rain$dates <- as.POSIXct(rain.all5$dates, tz = "Europe/Paris")
rain_omsev$dates <- as.POSIXct(rain_omsev$dates, format = "%Y-%m-%d %H:%M:%S", tz = "UTC")
rain_omsev <- rain_omsev[rain_omsev$dates >= "2019-09-01" & rain_omsev$dates <= "2025-02-01", ] #!!!!!!!(retirer 2019 ou pas ?? TODO)
rownames(rain_omsev) <- rain_omsev$dates
head(rain_omsev)
tail(rain_omsev)

# get location of each rain gauge
filename_loc <- paste0(data_folder, "omsev/loc_rain_gauges.csv")
location_gauges <- read.csv(filename_loc)
location_gauges$Station <- c("iem", "mse", "poly", "um", "cefe", "cnrs",
                             "crbm", "archiw", "archie", "um35", "chu1",
                             "chu2", "chu3", "chu4", "chu5", "chu6", "chu7",
                             "cines", "brives", "hydro")



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
rain_com <- df_comephore[df_comephore$date >= "2008-01-01", ]
rownames(rain_com) <- format(as.POSIXct(df_comephore$date), "%Y-%m-%d %H:%M:%S")


# get wind data
filename_wind <- paste0(data_folder, "wind/RADOME/wind_data_6min.csv")

filename_wind <- paste0(data_folder, "wind/data_gouv/wind_mtp.csv")
wind_mtp_h <- read.csv(filename_wind)
head(wind_mtp_h)

# When datetime do not have hours and minutes, add 00:00:00
wind_mtp_h$datetime <- ifelse(nchar(wind_mtp_h$datetime) == 10,
                            paste0(wind_mtp_h$datetime, " 00:00:00"),
                            wind_mtp_h$datetime)
# Convert datetime to POSIXct
wind_mtp_h$datetime <- as.POSIXct(wind_mtp_h$datetime,
                                format = "%Y-%m-%d %H:%M:%S", tz = "UTC")
head(wind_mtp_h$datetime)
year(wind_mtp_h$datetime)
# get only data between 2019 and 2025
wind_mtp_h <- wind_mtp_h[year(wind_mtp_h$datetime) >= 2019 & year(wind_mtp_h$datetime) <= 2025, ]
# Apply function to the wind_direction column
wind_mtp_h$cardDir <- sapply(wind_mtp_h$DD, convert_to_cardinal)
wind_mtp_h$cardDir <- as.character(wind_mtp_h$cardDir)  # Ensure it's character
wind_mtp_h$cardDir[is.na(wind_mtp_h$DD)] <- NA
wind_mtp_h <- wind_mtp_h[, c("datetime", "FF", "DD", "cardDir")]
colnames(wind_mtp_h) <- c("datetime", "speed", "direction", "cardDir")


# Convertir OMSEV en format long
omsev_long <- rain_omsev %>%
  pivot_longer(-dates, names_to = "Station", values_to = "rain") %>%
  left_join(location_gauges, by = "Station")

# Calculer barycentre par épisode ou par fenêtre 5-10min
bary_omsev <- omsev_long %>%
  group_by(dates) %>%
  summarise(
    rain_tot = sum(rain, na.rm = TRUE),
    lon_c = ifelse(rain_tot > 0, weighted.mean(Longitude, rain, na.rm = TRUE), NA_real_),
    lat_c = ifelse(rain_tot > 0, weighted.mean(Latitude, rain, na.rm = TRUE), NA_real_)
  ) %>%
  arrange(dates) %>%
  mutate(
    lon_next = lead(lon_c),
    lat_next = lead(lat_c),
    dt_h = as.numeric(difftime(lead(dates), dates, units="hours")),
    dist_km = distVincentyEllipsoid(cbind(lon_c, lat_c),
                                    cbind(lon_next, lat_next)) / 1000,
    speed_prec_kmh = dist_km / dt_h,
    bearing = bearing(cbind(lon_c, lat_c), cbind(lon_next, lat_next)),
    vx_prec_omsev = speed_prec_kmh * sin(bearing*pi/180),
    vy_prec_omsev = speed_prec_kmh * cos(bearing*pi/180)
  )

library(ggplot2)
library(geosphere) # pour distVincentyEllipsoid et bearing si besoin
library(dplyr)
# Filtrer les données valides

# keep only positive barycentre
bary_t <- bary_omsev %>% filter(!is.na(lon_c) & !is.na(lat_c))


bary_sample <- bary_t %>%
  dplyr::slice(1:50)  # juste un sous-ensemble pour tester

# Tracer
ggplot(bary_sample) +
  # pluie : cercles proportionnels à rain_tot
  geom_point(aes(x = lon_c, y = lat_c, size = rain_tot), color = "blue", alpha = 0.5) +
  # vecteurs d'advection
  geom_segment(aes(x = lon_c, y = lat_c,
                   xend = lon_c + vx_prec_omsev/50,  # factor to scale for plotting
                   yend = lat_c + vy_prec_omsev/50),
               arrow = arrow(length = unit(0.2,"cm")),
               color = "red") +
  scale_size_continuous(name = "Rain (mm)") +
  labs(title = "OMSEV advection vectors with rainfall intensity",
       x = "Longitude", y = "Latitude") +
  theme_minimal() +
  coord_fixed()

# COMEPHORE en format long : date, pixel_name, rain
come_long <- rain_com %>%
  pivot_longer(-date, names_to = "pixel_name", values_to = "rain") %>%
  left_join(loc_px, by = "pixel_name")


bary_com <- come_long %>%
  group_by(date) %>%
  summarise(
    rain_tot = sum(rain, na.rm = TRUE),
    lon_c = if (rain_tot > 0.5) weighted.mean(Longitude, rain, na.rm = TRUE) else NA_real_,
    lat_c = if (rain_tot > 0.5) weighted.mean(Latitude, rain, na.rm = TRUE) else NA_real_
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
    vx_prec_com = speed_prec_kmh * sin(bearing*pi/180),
    vy_prec_com = speed_prec_kmh * cos(bearing*pi/180)
  )
# rename date en dates pour faire la jointure
colnames(bary_com)[1] <- "dates"
head(bary_com)
head(bary_omsev)
# Exemple de jointure par episode_id
final_data <- bary_omsev %>%
  left_join(bary_com, by="dates", suffix=c("_omsev", "_com")) %>%
  mutate(
    vx_final = ifelse(!is.na(vx_prec_com), vx_prec_com, vx_prec_omsev),
    vy_final = ifelse(!is.na(vy_prec_com), vy_prec_com, vy_prec_omsev)
  )
head(final_data)

final_data <- final_data %>%
  left_join(wind_mtp_h %>% 
              mutate(u_kmh = -speed * sin((90 - direction) * pi / 180),
                     v_kmh = -speed * cos((90 - direction) * pi / 180)),
            by = c("dates" = "datetime")) %>%
  mutate(
    # Si barycentre COMEPHORE manquant, mix OMSEV + RADOME
    vx_final = ifelse(is.na(vx_prec_com),
                      0.7 * vx_prec_omsev + 0.3 * u_kmh,
                      vx_final),
    vy_final = ifelse(is.na(vy_prec_com),
                      0.7 * vy_prec_omsev + 0.3 * v_kmh,
                      vy_final)
  )

eta1 <- 1.0  # à estimer
eta2 <- 1.0  # à estimer

final_data <- final_data %>%
  mutate(
    vx_trans = eta1 * sign(vx_final) * abs(vx_final)^eta2,
    vy_trans = eta1 * sign(vy_final) * abs(vy_final)^eta2
  )

library(ggplot2)


ggplot(final_data, aes(x = vx_final, y = vy_final)) +
  geom_point(alpha = 0.3) +
  geom_hline(yintercept = 0, linetype="dashed") +
  geom_vline(xintercept = 0, linetype="dashed") +
  labs(title="Vents finaux : vx vs vy", x="vx (km/h)", y="vy (km/h)")


# wind rose with only ggplot2
ggplot(final_data, aes(x = vx_final, y = vy_final)) +
  geom_point(alpha = 0.3) +
  geom_hline(yintercept = 0, linetype="dashed") +
  geom_vline(xintercept = 0, linetype="dashed") +
  labs(title="Vents finaux : vx vs vy", x="vx (km/h)", y="vy (km/h)") +
  coord_fixed() +
  theme_minimal() 

final_data$speed_final <- sqrt(final_data$vx_final^2 + final_data$vy_final^2)
final_data$direction_final <- (atan2(final_data$vx_final, final_data$vy_final ) * 180 / pi) %% 360

# categrorize direction_final into 8 cardinal directions
final_data$cardDir_final <- cut(final_data$direction_final,
                                breaks = c(-22.5, 22.5, 67.5, 112.5, 157.5, 202.5, 247.5, 292.5, 337.5, 382.5),
                                labels = c("N", "NE", "E", "SE", "S", "SW", "W", "NW", "N"),
                                include.lowest = TRUE)  

# categrorize speed_final into 5 categories
final_data$speed_cat <- cut(final_data$speed_final,
                            breaks = c(-Inf, 1, 2, 5, 10, 15, Inf),
                            labels = c("0-1", "1-2", "2-5", "5-10", "10-15", ">15"),
                            include.lowest = TRUE)  

colnames(final_data)
library(ggplot2)
library(dplyr)

# Résumer les occurrences par direction et vitesse
wind_rose_df <- final_data %>%
  group_by(cardDir_final, speed_cat) %>%
  summarise(count = n(), .groups = "drop") %>%
  mutate(freq = count / sum(count) * 100)  # pourcentage

# remove NA
wind_rose_df <- wind_rose_df[!is.na(wind_rose_df$cardDir_final), ]
wind_rose_df <- wind_rose_df[!is.na(wind_rose_df$speed_cat), ]
# Plot rose des vents avec ggplot
ggplot(wind_rose_df, aes(x = cardDir_final, y = freq, fill = speed_cat)) +
  geom_bar(stat = "identity", width = 1, color = "black") +
  coord_polar(start = -pi/8) +  # aligne N vers le haut
  scale_fill_brewer(palette = "YlGnBu", name = "Vitesse (km/h)") +
  labs(title = "Wind Rose (final data)", x = "", y = "Frequency (%)") +
  theme_minimal() +
  theme(axis.text.y = element_blank(),
        axis.ticks = element_blank(),
        panel.grid.major = element_line(color = "gray80"))

