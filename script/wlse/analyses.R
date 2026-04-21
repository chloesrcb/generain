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


rain1min$dates <- as.POSIXct(rain1min$dates, tz = "Europe/Paris")
rain1min$dates <- with_tz(rain1min$dates, tzone = "UTC")

# remove brives, hydro, cines
rain1min <- rain1min[, !(colnames(rain1min) %in% c("brives", "hydro", "cines", "dates"))]
# count number of NAs in each column

na_counts <- sapply(rain1min, function(col) {
  sum(is.na(col))
})
print(na_counts)

# plot by site the na counts with frequency of NA values
na_freq <- sapply(rain1min, function(col) {
  mean(is.na(col))
})
print(na_freq)


na_freq_df <- data.frame(Site = names(na_freq), NA_Frequency = na_freq)
ggplot(na_freq_df, aes(x = Site, y = NA_Frequency)) +
  geom_bar(stat = "identity", fill = btfgreen) +
  theme_minimal() +
  xlab("Site") + ylab("Frequency of missing values") +
  btf_theme +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

folder <- paste0(im_folder, "rain/OMSEV/")
if (!dir.exists(folder)) {
  dir.create(folder, recursive = TRUE)
}
ggsave(paste0(folder, "na_frequency_by_site_1min.png"), width = 8, height = 6, dpi = 400)

# total number of NA values
total_na <- sum(is.na(rain1min))
print(paste("Total number of NA values:", total_na))

total_na_freq <- mean(is.na(rain1min))
print(paste("Overall frequency of NA values:", total_na_freq))


# count na within first and last non NA value for each column
na_within_range <- sapply(rain1min, function(col) {
  non_na_indices <- which(!is.na(col))
  if (length(non_na_indices) > 0) {
    first_non_na <- non_na_indices[1]
    last_non_na <- non_na_indices[length(non_na_indices)]
    sum(is.na(col[first_non_na:last_non_na]))
  } else {
    NA
  }
})
print(na_within_range)

# proportion of NA within first and last non NA value for each column
na_within_range_freq <- sapply(rain1min, function(col) {
  non_na_indices <- which(!is.na(col))
  if (length(non_na_indices) > 0) {
    first_non_na <- non_na_indices[1]
    last_non_na <- non_na_indices[length(non_na_indices)]
    sum(is.na(col[first_non_na:last_non_na])) / (last_non_na - first_non_na + 1)
  } else {NA
  }
})

print(na_within_range_freq)


# plot by site the proportion of NA values within first and last non NA value
na_within_range_freq_df <- data.frame(Site = names(na_within_range_freq), NA_Within_Range_Frequency = na_within_range_freq)
ggplot(na_within_range_freq_df, aes(x = Site, y = NA_Within_Range_Frequency)) +
  geom_bar(stat = "identity", fill = btfgreen) +
  theme_minimal() +
  xlab("Site") + ylab("Proportion of missing values") +
  btf_theme +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(paste0(folder, "na_within_range_frequency_by_site_1min.png"), width = 8, height = 6, dpi = 400)


# count number of NA within range of first and last non NA value for each column for all sites
total_na_within_range <- sum(na_within_range, na.rm = TRUE)
# proportion of NA within range of first and last non NA value for all sites
total_na_within_range_freq <- sum(na_within_range, na.rm = TRUE)/ sum(sapply(rain1min, function(col) {
  non_na_indices <- which(!is.na(col))
  if (length(non_na_indices) > 0) {
    first_non_na <- non_na_indices[1]
    last_non_na <- non_na_indices[length(non_na_indices)]
    return(last_non_na - first_non_na + 1)
  } else {
    return(0)
  }
}))

# # Loop through all columns except "Dates"
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
# rain$mse
# Get distances matrix
dist_mat <- get_dist_mat(location_gauges)
df_dist <- reshape_distances(dist_mat)
# x11()
plot(df_dist$value)
max(df_dist$value)
