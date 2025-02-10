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

# wind_data_1990_1999 <- read.csv(gzfile("./data/wind/data_gouv/H_34_1990-1999.csv.gz"), sep = ";")
# wind_data_2000_2009 <- read.csv(gzfile("./data/wind/data_gouv/H_34_2000-2009.csv.gz"), sep = ";")
# wind_data_2010_2019 <- read.csv(gzfile("./data/wind/data_gouv/H_34_2010-2019.csv.gz"), sep = ";")
# wind_data_2020_2023 <- read.csv(gzfile("./data/wind/data_gouv/H_34_previous-2020-2023.csv.gz"), sep = ";")
# # head(wind_data_2010_2019)
# # head(wind_data_2020_2023)

# # put all wind data together
# wind_data <- rbind(wind_data_1990_1999, wind_data_2000_2009,
#                                 wind_data_2010_2019, wind_data_2020_2023)
# head(wind_data)
# tail(wind_data)
# # FF: force du vent sur 10 minutes mesuré à 10 metres en m/s et 1/10
# # DD: direction de FF en degrésen rose de 0 à 360
# # FXY: valeur maximale de FF dans l’heure (en m/s et 1/10)
# # DXY: direction de FXY (rose de 360)
# # HXY         : heure de FXY (hhmm)

# wind_df <- wind_data
# get stations names
# stations <- unique(wind_df$NOM_USUEL)
# length(stations)
# # get unique locations by station
# stations_loc <- wind_df %>% select(NOM_USUEL, LAT, LON) %>% unique()

# # plot all stations on a map around Montpellier
# # get the map
# library(leaflet)
# # Create a map centered on Montpellier
# leaflet() %>%
#   addTiles() %>%
#   addMarkers(lng = stations_loc$LON, lat = stations_loc$LAT, popup = stations_loc$NOM_USUEL)



# # keep only "Montpellier-Aéroport" and "Grabels"
# wind_df <- wind_df[wind_df$NOM_USUEL %in% c("MONTPELLIER-AEROPORT"),]

# # convert AAAAMMJJHH to datetime
# wind_df$AAAAMMJJHH <- as.POSIXct(strptime(as.character(wind_df$AAAAMMJJHH), format = "%Y%m%d%H"), tz = "UTC")
# # convert HXY to time
# # wind_df$HXY <- as.POSIXct(strptime(as.character(wind_df$HXY), format = "%H%M"), tz = "UTC")
# # rename columns
# colnames(wind_df)[colnames(wind_df) == "AAAAMMJJHH"] <- "datetime"
# # remove AAAAMMJJHH
# wind_mtp <- wind_df[,c("datetime", "FF", "DD", "FXY", "DXY", "HXY")]
# head(wind_mtp)

# # save to csv
# write.csv(wind_mtp, file = "./data/wind/data_gouv/wind_mtp.csv", row.names = FALSE)

# get csv
wind_mtp <- read.csv("./data/wind/data_gouv/wind_mtp.csv")
head(wind_mtp)
# remove rows with NA for DD
wind_mtp <- wind_mtp[!is.na(wind_mtp$DD),]

# convert datetime to POSIXct
wind_mtp$datetime <- as.POSIXct(wind_mtp$datetime, tz = "UTC")

convert_to_cardinal <- function(degrees) {
  if (is.na(degrees)) {
    return(NA)  # Return NA if the input is NA
  }
  
  directions <- c("N", "NE", "E", "SE", "S", "SW", "W", "NW", "N")  
  breaks <- c(0, 22.5, 67.5, 112.5, 157.5, 202.5, 247.5, 292.5, 337.5, 360)  
  
  return(directions[findInterval(degrees, breaks, rightmost.closed = TRUE)])
}

# Apply function to the DD column
wind_mtp$cardDir <- sapply(wind_mtp$DD, convert_to_cardinal)
wind_mtp$cardDir <- as.character(wind_mtp$cardDir)  # Ensure it's character, not factor
wind_mtp$cardDir[is.na(wind_mtp$DD)] <- NA  

# Filter data for the years 2018 to 2023
wind_mtp_years <- wind_mtp %>%
  filter(datetime >= as.POSIXct("2018-01-01", tz = "UTC") & 
         datetime <= as.POSIXct("2023-12-31", tz = "UTC"))


library(lubridate)  # For date handling

# Extract year from datetime
wind_mtp_years$year <- year(as.POSIXct(wind_mtp_years$datetime, format="%Y-%m-%d %H:%M:%S", tz="UTC"))

# Ensure cardDir is a factor
wind_mtp_years$cardDir <- factor(wind_mtp_years$cardDir, 
                                 levels = c("N", "NE", "E", "SE", "S", "SW", "W", "NW"))

ggplot(wind_mtp_years, aes(x = cardDir)) +
  geom_bar(color = "#6965659d", fill=btfgreen, alpha= 0.5, width = 1) +  # Bar plot counting occurrences
  coord_polar(start = 15*pi/8) +  # Set North at the top
  labs(x = "Wind direction", y = "Count", title = "") +
  facet_wrap(~ year) +  # One wind rose per year
  theme_minimal() +
  theme(axis.text.x = element_text(size = 12, face = "bold"))


filename <- paste0(im_folder, 
                  "wind/datagouv/wind_rose_comparison_2018_2023.png")
ggsave(filename, width = 30, height = 15, units = "cm")
