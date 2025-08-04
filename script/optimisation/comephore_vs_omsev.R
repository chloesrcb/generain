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
rain$dates <- as.POSIXct(rain.all5$dates, tz = "UTC")
rownames(rain) <- rain$dates

# Filter dates and remove NA rows from all columns
start_date <- as.POSIXct("2023-09-20 08:00:00", tz = "UTC")
end_date <- as.POSIXct("2023-09-21 01:00:00", tz = "UTC")
dateshours_str <- paste("from_", format(start_date, "%Y-%m-%d_%H%M"),
                   "_to_", format(end_date, "%Y-%m-%d_%H%M"), sep = "")
month_year_str <- paste(format(start_date, "%Y-%m"), collapse = "_")

subrain <- rain
subrain <- subrain %>%
  filter(dates >= start_date & dates <= end_date)

# remove brives, cines, hydro
subrain <- subrain[, !(colnames(subrain) %in% c("brives", "cines", "hydro"))]

# Get all rain gauge columns (exclude 'dates')
rain_gauges <- setdiff(names(subrain), "dates")

foldername <- paste0(im_folder, "rain/OMSEV/",
                     "episodes_", month_year_str, "/")

# put every rain gauge time series on the same plot
# melt the data frame to long format
library(reshape2)
rain_long <- melt(subrain, id.vars = "dates", variable.name = "Station",
                  value.name = "Rain")
colnames(rain_long) <- c("dates", "Station", "Rain")
# remove na
rain_long <- rain_long %>%
  filter(!is.na(Rain))
# Create the plot
p <- ggplot(rain_long, aes(x = dates, y = Rain, color = Station)) +
  geom_line() +
  labs(
    title = "",
    x = "Date",
    y = "Rain (mm)"
  ) +
  theme_minimal() +
  theme(legend.position = "none") 
p
# Save the plot
filename <- paste0(foldername, "5min_", dateshours_str, ".png")
ggsave(
  filename = filename,
  plot = p,
  width = 10,
  height = 6,
  dpi = 300
)

# remove na in subrain
subrain <- na.omit(subrain)
# Aggregate the data to 1h intervals (5 min to 1h)
rain_1h <- subrain %>%
  group_by(dates = ceiling_date(dates, "hour")) %>%
  summarise(across(everything(), ~ sum(.x, na.rm = TRUE))) %>%
  ungroup()

library(reshape2)
rain_long <- melt(as.data.frame(rain_1h), id.vars = "dates",
                    variable.name = "Station", value.name = "Rain")
colnames(rain_long) <- c("dates", "Station", "Rain")
# remove na
rain_long <- rain_long %>%
  filter(!is.na(Rain))
# Create the plot
p <- ggplot(rain_long, aes(x = dates, y = Rain, color = Station)) +
  geom_line() +
  labs(
    title = "",
    x = "Date",
    y = "Rain (mm)"
  ) +
  theme_minimal() +
  theme(legend.position = "none") 
p
# Save the plot
filename <- paste0(foldername, "1h_", dateshours_str, ".png")
ggsave(
  filename = filename,
  plot = p,
  width = 10,
  height = 6,
  dpi = 300
)


########
# Check comephore values
filename_com <- paste0(data_folder, "comephore/comephore_2008_2024.csv")
comephore_raw <- read.csv(filename_com, sep = ",")

# get comephore values for the selected episodes datetime
df_comephore <- as.data.frame(comephore_raw)
colnames(df_comephore)[1] <- "date"
# if there is a "Y-M-D" change to "Y-M-D 00:00:00"
df_comephore$date <- ifelse(nchar(df_comephore$date) == 10,
                            paste0(df_comephore$date, " 00:00:00"),
                            df_comephore$date)

df_comephore$date <- as.POSIXct(df_comephore$date,
                                format = "%Y-%m-%d %H:%M:%S", tz = "UTC")
head(df_comephore)

# add one day to start_date and end_date
start_date <- start_date
end_date <- end_date + lubridate::hours(12)
df_comephore <- df_comephore[df_comephore$date >= start_date &
                              df_comephore$date <= end_date, ]


library(reshape2)
rain_long <- melt(as.data.frame(df_comephore), id.vars = "date",
                  variable.name = "Pixel", value.name = "Rain")

colnames(rain_long) <- c("dates", "Pixel", "Rain")
# remove na
rain_long <- rain_long %>%
  filter(!is.na(Rain))
# Create the plot
p <- ggplot(rain_long, aes(x = dates, y = Rain, color = Pixel)) +
  geom_line() +
  labs(
    title = "",
    x = "Date",
    y = "Rain (mm)"
  ) +
  theme_minimal() +
  theme(legend.position = "none")
p

foldername_com <- paste0(im_folder, "rain/comephore/",
                     "episodes_", month_year_str, "/")
filename <- paste0(foldername_com, "all_pixels", "_",
                   dateshours_str, ".png")
# Save the plot
ggsave(
  filename = filename,
  plot = p,
  width = 10,
  height = 6,
  dpi = 300
)
