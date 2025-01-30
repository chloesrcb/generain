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

wind_df <- wind_data_2010_2019 
# get stations names
stations <- unique(wind_df$NOM_USUEL)
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
