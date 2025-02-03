library(generain)
library(reshape2)
library(ggplot2)
source("load_libraries.R")
library(kableExtra)
library(extRemes)
library(bbmle)
library(ismev)
library(extRemes)
library(evd)
library(latex2exp)
library(geosphere)
library(dplyr)

btf_theme <- theme_minimal() +
  theme(axis.text.x = element_text(size =  6, angle = 0),
        axis.text.y = element_text(size = 6),
        axis.title.x = element_text(size = 8),
        axis.title.y = element_text(size = 8),
        legend.title = element_text(size = 8),
        legend.text = element_text(size = 6),
        title = element_text(size = 10),
        axis.line = element_blank(),  # Remove axis lines
        panel.border = element_blank(),  # Remove plot border
        panel.background = element_rect(fill = "transparent", color = NA),
        # Remove plot background
        axis.text = element_blank(),  # Remove axis text
        axis.ticks = element_blank(),  # Remove axis ticks
        plot.background = element_rect(fill = "transparent", color = NA),
        panel.grid = element_line(color = "#5c595943"))

# my green color
btfgreen <- "#69b3a2"

# fromhttps://www.data.gouv.fr/fr/datasets/donnees-climatologiques-de-base-horaires/

wind_data_2000_2009 <- read.csv(gzfile("./data/wind/data_gouv/H_34_2000-2009.csv.gz"), sep = ";")
wind_data_2010_2019 <- read.csv(gzfile("./data/wind/data_gouv/H_34_2010-2019.csv.gz"), sep = ";")
# wind_data_2020_2023 <- read.csv(gzfile("./data/wind/data_gouv/H_34_previous-2020-2023.csv.gz"), sep = ";")
head(wind_data_2010_2019)

# FF: force du vent sur 10 minutes mesuré à 10 metres en m/s et 1/10
# DD: direction de FF en degrésen rose de 0 à 360
# FXY: valeur maximale de FF dans l’heure (en m/s et 1/10)
# DXY: direction de FXY (rose de 360)
# HXY         : heure de FXY (hhmm)

wind_df <- wind_data_2010_201
# get stations names
stations <- unique(wind_df$NOM_USUEL)
length(stations)
# get unique locations by station
stations_loc <- wind_df %>% select(NOM_USUEL, LAT, LON) %>% unique()

# plot all stations on a map around Montpellier
# get the map
library(leaflet)
# Create a map centered on Montpellier
leaflet() %>%
  addTiles() %>%
  addMarkers(lng = stations_loc$LON, lat = stations_loc$LAT, popup = stations_loc$NOM_USUEL)



# keep only "Montpellier-Aéroport" and "Grabels"
wind_df <- wind_df[wind_df$NOM_USUEL %in% c("MONTPELLIER-AEROPORT"),]

# convert AAAAMMJJHH to datetime
wind_df$AAAAMMJJHH <- as.POSIXct(strptime(as.character(wind_df$AAAAMMJJHH), format = "%Y%m%d%H"), tz = "UTC")
# convert HXY to time
# wind_df$HXY <- as.POSIXct(strptime(as.character(wind_df$HXY), format = "%H%M"), tz = "UTC")
# rename columns
colnames(wind_df)[colnames(wind_df) == "AAAAMMJJHH"] <- "datetime"
# remove AAAAMMJJHH
wind_mtp <- wind_df[,c("datetime", "FF", "DD", "FXY", "DXY", "HXY")]
head(wind_mtp)

# for one year 2019
wind_mtp_2019 <- wind_mtp[wind_mtp$datetime >= as.POSIXct("2019-01-01", tz = "UTC") & wind_mtp$datetime <= as.POSIXct("2019-12-31", tz = "UTC"),]

# plot wind speed
ggplot(wind_mtp_2019, aes(x = datetime, y = FF)) +
  geom_line(color = btfgreen) +
  labs(title = "Wind speed at Montpellier-Aéroport",
       x = "Date",
       y = "Wind speed (m/s)") +
  btf_theme

library(ggforce)
# plot wind direction
# Création de la rose des vents
ggplot(wind_mtp_2019, aes(x = DD_rad, y = FF)) +
  geom_histogram(stat = "identity", aes(fill = FF), binwidth = pi/8) +
  coord_polar(start = -pi/2) +
  scale_x_continuous(breaks = seq(0, 2*pi, pi/4), 
                     labels = c("E", "NE", "N", "NW", "W", "SW", "S", "SE")) +
  labs(title = "Rose des vents - Montpellier", x = "Direction du vent", y = "Fréquence") +
  theme_minimal()

library(dplyr)

# for one year 2018
wind_mtp_2018 <- wind_mtp[wind_mtp$datetime >= as.POSIXct("2018-01-01", tz = "UTC") & wind_mtp$datetime <= as.POSIXct("2018-12-31", tz = "UTC"),]

wind_mtp_2018 <- wind_mtp_2018 %>%
  mutate(direction_bin = cut(DD, breaks = seq(0, 360, by = 30), include.lowest = TRUE))

wind_data <- wind_mtp_2018 %>%
  group_by(direction_bin) %>%
  summarise(count = n())

angle_labels <- c("E", "NE", "N", "NW", "W", "SW", "S", "SE")
angle_breaks <- seq(0, 315, by = 45)

ggplot(wind_data, aes(x = as.numeric(direction_bin), y = count, fill = count)) +
  geom_bar(stat = "identity", width = 1, color = "black") +
  coord_polar(start = -pi/2) +
  scale_x_continuous(breaks = angle_breaks, labels = angle_labels) +
  scale_fill_viridis_c(option = "rocket") +
  labs(x = "Direction of wind", y = "Frequence") +
  theme_minimal()

# save
ggsave("wind_rose_2018.png", width = 10, height = 10, units = "cm")

wind_mtp_2019 <- wind_mtp_2019 %>%
  mutate(direction_bin = cut(DD, breaks = seq(0, 360, by = 30), include.lowest = TRUE))

wind_data <- wind_mtp_2019 %>%
  group_by(direction_bin) %>%
  summarise(count = n())

angle_labels <- c("E", "NE", "N", "NW", "W", "SW", "S", "SE")
angle_breaks <- seq(0, 315, by = 45) 

ggplot(wind_data, aes(x = as.numeric(direction_bin), y = count, fill = count)) +
  geom_bar(stat = "identity", width = 1, color = "black") +
  coord_polar(start = -pi/2) +
  scale_x_continuous(breaks = angle_breaks, labels = angle_labels) +
  scale_fill_viridis_c(option = "rocket") + 
  labs(x = "Direction of wind", y = "Frequence") +
  theme_minimal()


ggsave("wind_rose_2019.png", width = 10, height = 10, units = "cm")


library(ggplot2)
library(dplyr)

wind_mtp_2018 <- wind_mtp[wind_mtp$datetime >= as.POSIXct("2018-01-01", tz = "UTC") & wind_mtp$datetime <= as.POSIXct("2018-12-31", tz = "UTC"),]
wind_mtp_2018 <- wind_mtp_2018 %>%
  mutate(direction_bin = cut(DD, breaks = seq(0, 360, by = 30), include.lowest = TRUE),
         year = "2018")

wind_mtp_2019 <- wind_mtp[wind_mtp$datetime >= as.POSIXct("2019-01-01", tz = "UTC") & wind_mtp$datetime <= as.POSIXct("2019-12-31", tz = "UTC"),]
wind_mtp_2019 <- wind_mtp_2019 %>%
  mutate(direction_bin = cut(DD, breaks = seq(0, 360, by = 30), include.lowest = TRUE),
         year = "2019")

wind_data_combined <- bind_rows(wind_mtp_2018, wind_mtp_2019) %>%
  group_by(direction_bin, year) %>%
  summarise(count = n(), .groups = "drop")

angle_labels <- c("E", "NE", "N", "NW", "W", "SW", "S", "SE")
angle_breaks <- seq(0, 315, by = 45)

ggplot(wind_data_combined, aes(x = as.numeric(direction_bin), y = count, fill = count)) +
  geom_bar(stat = "identity", width = 1, color = "black") +
  coord_polar(start = -pi/2) +
  scale_x_continuous(breaks = angle_breaks, labels = angle_labels) +
  scale_fill_viridis_c(option = "rocket", limits = c(0, max(wind_data_combined$count))) + 
  facet_wrap(~year) +
  labs(x = "Direction of wind", y = "Frequence") +
  btf_theme + 
  theme(
    legend.text = element_text(size = 5),    # Reduce the size of the legend text
    legend.title = element_text(size = 5),   # Reduce the size of the legend title
    legend.key.size = unit(0.2, "cm"),        # Reduce the size of the legend symbols
    legend.position = "right"                # Position the legend on the right (optional)
  )

ggsave("../phd_extremes/wind_rose_comparison_2018_2019.png", width = 10, height = 5, units = "cm")
