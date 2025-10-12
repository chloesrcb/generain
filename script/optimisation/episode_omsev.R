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
# first_date <- as.POSIXct("2019-10-01", tz = "UTC")
# last_date <- as.POSIXct("2024-01-31", tz = "UTC")
# rain <- rain[rain$dates >= first_date & rain$dates <= last_date, ]
rownames(rain) <- rain$dates
head(rain)
tail(rain)

# Filter dates and remove NA rows from all columns
first_date <- as.POSIXct("2022-09-06 23:00:00", tz = "UTC")
last_date <- as.POSIXct("2022-09-07 05:00:00", tz = "UTC")
subrain <- rain %>%
  filter(dates >= first_date & dates <= last_date)

# remove brives, cines, hydro
subrain <- subrain[, !(colnames(subrain) %in% c("brives", "cines", "hydro"))]

# Get all rain gauge columns (exclude 'dates')
rain_gauges <- setdiff(names(subrain), "dates")
foldername <- paste0(im_folder, "rain/OMSEV/episodes_6_8_sept_2022/")

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
  theme(legend.position = "right") 
p
# Save the plot
ggsave(
  filename = paste0(foldername, "all_gauges_5min_6_7_sept_2022_from1hto4h.png"),
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

plot(subrain$dates, subrain$mse, type = "l", col = btfgreen,
     xlab = "Date", ylab = "Rain (mm)", main = "Rainfall at MSE station (5 min aggregation)")
     
plot(rain_1h$dates, rain_1h$mse, type = "l", col = btfgreen,
     xlab = "Date", ylab = "Rain (mm)", main = "Rainfall at MSE station (1h aggregation)")


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
  theme(legend.position = "right") 
p
# Save the plot
ggsave(
  filename = paste0(foldername, "all_gauges_1h_6_7_sept_2022_from1hto4h.png"),
  plot = p,
  width = 10,
  height = 6,
  dpi = 300
)



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

rain1min$dates <- as.POSIXct(rain1min$dates, tz = "UTC")
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

# in rain remove when all data are NA<
rain <- rain[rowSums(is.na(rain)) < ncol(rain), ]

head(rain$mse)
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


sites_names <- colnames(rain)
# in rain remove when all data are NA<
rain <- rain[rowSums(is.na(rain)) < ncol(rain), ]

# get rain for only september 2019
# subrain <- rain[rain_clean$dates >= "2019-09-05" & rain_clean$dates <= "2019-10-01", ]
# plot(subrain$mse)

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
breaks <- seq(floor(min(df_threshold$u_s0)), ceiling(max(df_threshold$u_s0)), by = 0.1)

ggplot(df_threshold, aes(x = u_s0)) +
  geom_histogram(breaks = breaks, fill = btfgreen, alpha = 0.5, color = "#5f5d5d") +
  btf_theme +
  xlab(TeX(paste0("Threshold for quantile $q = ", q, "$"))) +
  ylab("Count")
filename <- paste(im_folder, "optim/omsev/threshold_histogram_q",
                  q * 1000, "_delta", delta, "_dmin", min_spatial_dist,
                  ".png", sep = "")
ggsave(filename, width = 20, height = 15, units = "cm")

n_episodes <- length(selected_points$s0)
print(n_episodes)

list_episodes_points <- get_extreme_episodes(selected_points, rain,
                                     episode_size = episode_size, unif = FALSE,
                                     beta = 0)

list_episodes <- list_episodes_points$episodes
episode <- list_episodes[[1]]

# check that for each episode ep, for s0=s0_list[ep] and t0=t0_list[ep] 
# we have an excess above the threshold
for (ep in 1:length(list_episodes)) {
  s0 <- selected_points$s0[ep]
  t0 <- 1
  u_s0 <- selected_points$u_s0[ep]
  episode <- list_episodes[[ep]]
  if (episode[t0, s0] <= u_s0) {
    stop(paste("Excess is not above threshold for s0 =", s0, "and t0 =", t0))
  }
}


plots_list <- list()
s0_list <- unique(selected_points$s0)

for (s0 in s0_list) {
  episode_indices <- which(selected_points$s0 == s0)
  episodes_for_site <- list_episodes[episode_indices]
  u_site <- unique(selected_points$u_s0[episode_indices])  # seuil pour ce site

  # Concaténer les épisodes en format long
  episodes_combined <- rbindlist(
    lapply(seq_along(episodes_for_site), function(i) {
      df <- as.data.frame(episodes_for_site[[i]])
      data.frame(
        Time = 0:(nrow(df) - 1),
        Value = df[[s0]],
        Episode = paste0("Episode ", i)
      )
    })
  )

  # Créer le plot
  p <- ggplot(episodes_combined, aes(x = Time, y = Value, color = Episode)) +
    geom_line(size = 1, show.legend = FALSE, alpha = 0.7) +  # Cache la légende
    geom_hline(yintercept = u_site, color = "red", linetype = "dashed", size = 1.2) +
    labs(title = paste("Site", s0),
        x = "Relative time", y = "Rainfall (mm)") +
    theme_minimal()

  # Stocker le plot
  plots_list[[s0]] <- p
}


print(plots_list[["cnrs"]])  # Remplace "mse" par le nom du site souhaité


folder_episode <- paste0(im_folder, "optim/omsev/episodes_plot/q",
                         q * 100, "_delta", delta, "_dmin", min_spatial_dist, "/")
if (!dir.exists(folder_episode)) {
  dir.create(folder_episode, recursive = TRUE)
}

for (s0 in names(plots_list)) {
  filename <- paste0(folder_episode, "episodes_", s0, ".png")
  ggsave(filename = filename, plot = plots_list[[s0]],
         width = 8, height = 4)
}


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

unique(time_lookup$t0_date) # check unique t0 values
sort(unique(time_lookup$month))


# distribution of months
month_counts <- table(time_lookup$month)
month_df <- as.data.frame(month_counts) 
colnames(month_df) <- c("Month", "Count")
month_df$Month <- factor(month_df$Month, levels = 1:12,
                         labels = c("Jan", "Feb", "Mar", "Apr", "May",
                                    "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"))
ggplot(month_df, aes(x = Month, y = Count)) +
  geom_bar(stat = "identity", fill = btfgreen, alpha = 0.5, color = "#5f5d5d") +
  btf_theme +
  xlab("Month") +
  ylab("Count") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# save histo 
filename <- paste(im_folder, "optim/omsev/months_histogram_q",
                  q * 100,
                  ".png", sep = "")
ggsave(filename, width = 20, height = 15, units = "cm")

head(time_lookup)
# selected_points$t0_date <- time_lookup$t0_date[match(selected_points$t0, time_lookup$t0)]
selected_points$t0_date <- as.POSIXct(selected_points$t0_date, format = "%Y-%m-%d %H:%M:%S", tz = "UTC")
# get month of t0_date inside selected_points
selected_points$month <- month(selected_points$t0_date)

# plot selected points t0 dates
ggplot(selected_points, aes(x = t0_date)) +
  geom_histogram(bins = 50, fill = btfgreen, alpha = 0.5, color = "#5f5d5d") +
  btf_theme +
  xlab("t0 Date") +
  ylab("Count") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# plot months of t0 dates
ggplot(selected_points, aes(x = factor(month))) +
  geom_bar(fill = btfgreen, alpha = 0.5, color = "#5f5d5d") +
  btf_theme +
  xlab("Month") +
  ylab("Count") +
  scale_x_discrete(labels = c("Jan", "Feb", "Mar", "Apr", "May",
                               "Jun", "Jul", "Aug", "Sep", "Oct",
                               "Nov", "Dec")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# save plot
filename <- paste(im_folder, "optim/omsev/selected_episodes_months_histogram_q",
                  q * 100, "_delta", delta, "_dmin", min_spatial_dist,
                  ".png", sep = "")

ggsave(filename, width = 20, height = 15, units = "cm")

# if there is "Y-M-D" change to "Y-M-D 00:00:00"
selected_points$t0_date <- ifelse(
  nchar(format(selected_points$t0_date, "%Y-%m-%d %H:%M:%S")) == 10,
  paste0(format(selected_points$t0_date, "%Y-%m-%d"), " 00:00:00"),
  format(selected_points$t0_date, "%Y-%m-%d %H:%M:%S")
)

datetimes <- selected_points$t0_date



########
# Check comephore values
filename_com <- paste0(data_folder, "comephore/comephore_2008_2024.csv")
# filename_com <- paste0(data_folder, "comephore/zoom_5km.csv")
comephore_raw <- read.csv(filename_com, sep = ",")
filename_loc <- paste0(data_folder, "comephore/loc_px.csv")
loc_px <- read.csv(filename_loc, sep = ",")


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

# Take only data from 6 to 7 September 2022
start_date <- as.POSIXct("2022-09-06 00:00:00", tz = "UTC")
end_date <- as.POSIXct("2022-09-10 20:00:00", tz = "UTC")
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

foldername <- paste0(im_folder, "rain/comephore/episodes_6_8_sept_2022/")
# Save the plot
ggsave(
  filename = paste0(foldername, "all_pixels_6_9_sept_2022_from0to20h.png"),
  plot = p,
  width = 10,
  height = 6,
  dpi = 300
)

# zoom comephore values from 6 September 2022 18:00 to 7 September 2022 1:00
start_time <- as.POSIXct("2022-09-06 18:00:00", tz = "UTC")
end_time <- as.POSIXct("2022-09-07 01:00:00", tz = "UTC")
df_comephore <- df_comephore[df_comephore$date >= start_time & 
                              df_comephore$date <= end_time, ]

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

foldername <- paste0(im_folder, "rain/comephore/episodes_6_8_sept_2022/")
# Save the plot
ggsave(
  filename = paste0(foldername, "all_pixels_6_7_sept_2022_timezoom.png"),
  plot = p,
  width = 10,
  height = 6,
  dpi = 300
)





# Take only data after 2019 october 1st
df_comephore <- df_comephore[df_comephore$date >= "2019-10-01" & 
                              df_comephore$date <= "2024-01-31", ]
# # put date in index
# rownames(df_comephore) <- df_comephore$date
tail(df_comephore)
# comephore <- df_comephore[-1] # remove dates column


library(lubridate)

# For all episodes, check that comephore_episode has at least one value > 0 (excluding the date column)
for (episode_id in seq_along(list_episodes)) {
    episode <- list_episodes[[episode_id]]
    dates_episode <- rownames(episode)
    dates_episode <- as.POSIXct(dates_episode, format = "%Y-%m-%d %H:%M:%S", tz = "UTC")
    # Get the range from the previous day to the next day, by hour
    min_hour <- floor_date(min(dates_episode), unit = "day") - hours(1)
    max_hour <- ceiling_date(max(dates_episode), unit = "day") + hours(23)
    hours_episode <- seq(from = min_hour, to = max_hour, by = "hour")
    
    comephore_episode <- df_comephore[df_comephore$date %in% hours_episode, ]
    
    # Exclude the date column and check if any value > 0
    comephore_values <- comephore_episode[, !(colnames(comephore_episode) %in% "date"), drop = FALSE]
    if (!any(comephore_values > 0, na.rm = TRUE)) {
        warning(paste("No comephore value > 0 for episode", episode_id))
    }
}


episode <- list_episodes[[216]]


