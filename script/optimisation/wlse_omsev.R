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
# rain$dates <- as.POSIXct(rain.all5$dates, tz = "Europe/Paris")
rain$dates <- with_tz(rain$dates, tzone = "UTC")
rain <- rain[rain$dates >= "2019-01-01" & rain$dates <= "2024-02-01", ]
rownames(rain) <- rain$dates
head(rain)
tail(rain)
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
# rain$mse
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

# Spatial chi WLSE #############################################################
library(knitr)
library(kableExtra)
library(classInt)

make_distance_classes <- function(distances,
                                  method = "quantile",
                                  n = 10,
                                  pairs_per_class = NULL) {
  distances <- sort(distances)
  if (method == "quantile") {
    breaks <- quantile(distances, probs = seq(0, 1, length.out = n + 1))
  } else if (method == "equal") {
    breaks <- seq(min(distances), max(distances), length.out = n + 1)
  } else if (method == "jenks") {
    breaks <- classIntervals(distances, n = n, style = "jenks")$brks
  } else if (method == "fixed_pairs") {
    if (is.null(pairs_per_class)) stop("Need pairs_per_class")
    n_classes <- floor(length(distances) / pairs_per_class)
    indices <- seq(1, length(distances), by = pairs_per_class)
    breaks <- distances[indices]
    breaks <- c(breaks, max(distances) + 10)
  } else if (method == "log") {
    breaks <- exp(seq(log(min(distances)), log(max(distances)),
                                                          length.out = n + 1))
  } else {
    stop("Unknown method")
  }
  return(unique(breaks))
}

df_dist_order <- df_dist %>%
  filter(value > 0) %>%
  arrange(value)
# Supposons que df_dist_order a déjà été filtré avec df_dist_order$value > 0
distances <- df_dist_order$value

# Choisir méthode ici :
radius <- make_distance_classes(distances,
                                method = "quantile",
                                n = 12,
                                pairs_per_class = 20)
dist_counts <- table(cut(distances, breaks = radius))
df_hist <- data.frame(Interval = names(dist_counts), Count = as.vector(dist_counts))
df_hist$Breaks <- gsub("e\\+0.", "0", df_hist$Interval)
df_hist$Breaks <- gsub("\\.", "", df_hist$Breaks)

histradius <- ggplot(df_hist, aes(x = Interval, y = Count)) +
  btf_theme +
  geom_bar(stat = "identity", fill = btfgreen, alpha = 0.8) +
  xlab("Spatial lag") +
  ylab("Pair count") +
  theme(axis.text.x = element_text(angle = 45)) +
  scale_x_discrete(labels = df_hist$Breaks)

# save histogram
filename <- paste(im_folder, "WLSE/omsev/spatial/histogram_spatial_lag_",
                 length(radius) - 1, ".pdf", sep = "")
ggsave(filename, plot = histradius, width = 20, height = 15, units = "cm",
       dpi = 600, device = "pdf")

# Create matrix of radius
rad_mat <- dist_mat

for (i in 2:length(radius)) {
  curr_radius <- radius[i]
  prev_radius <- radius[i - 1]
  rad_mat[dist_mat >= prev_radius & dist_mat < curr_radius] <- curr_radius
}

# rad_mat <- rad_mat / 1000 # convert to km
# Assigner NA aux distances non classées (>= dernier rayon)
rad_mat[dist_mat >= max(radius)] <- NA
rad_mat[dist_mat == 0] <- NA

# Triangulaire
rad_mat[lower.tri(rad_mat)] <- NA
# rad_mat_sym <- rad_mat
# rad_mat_sym[lower.tri(rad_mat)] <- t(rad_mat)[lower.tri(rad_mat)]
hmax <- max(radius) # distance max in absolute value...
nb_col <- length(unique(radius)) # we exclude 0 ?

# quantiles
q_spa_vals <- c(0.90, 0.91, 0.92, 0.93, 0.94, 0.95, 0.96, 0.97, 0.98, 0.99)
q_temp_vals <- c(0.90, 0.91, 0.92, 0.93, 0.94, 0.95)
tmax <- 10

foldername <- paste0(data_folder, "/omsev/WLSE/")
dir.create(foldername, recursive = TRUE, showWarnings = FALSE)

spa_results <- list()


# Spatial loop =================================================================
for (q_no0_spa in q_spa_vals) {

  chispa_df <- spatial_chi(rad_mat, rain,
                         quantile = q_no0_spa, zeros = FALSE, mid = TRUE)


  etachispa_df <- data.frame(
    chi = eta(chispa_df$chi),
    lagspa = log(chispa_df$lagspa)
  )

  # Plot chi(h, 0)
  chispa_plot <- ggplot(chispa_df, aes(lagspa, chi)) +
        btf_theme +
        geom_point(col = btfgreen, size = 4) +
        xlab(TeX(paste0("$||", expression(bold("h")), "||$", " (m)"))) +
        ylab(TeX(paste0("$\\widehat{\\chi}(", expression(bold("h")), ", 0)$"))) +
        theme(axis.ticks = element_blank(),
              axis.line = element_blank(),
              panel.border = element_blank(),
              legend.text = element_text(size = 15),
              legend.title = element_text(size = 15),
              legend.position = "right",
              legend.margin = margin(0.5, 0.5, 0.5, 0, "cm"),
              panel.grid = element_line(color =  "#5c595943")) +
        ylim(0, 1) +
        theme(panel.grid.major = element_line(color = "white"),
              panel.grid.minor = element_line(color = "white"))
  
  ggsave(paste0(im_folder, "WLSE/comephore/full_spatial_chi_", q_no0_spa, ".pdf"),
         plot = chispa_plot, width = 20, height = 15, units = "cm")
  
  # Estimation WLSE
  spa_result <- get_estimate_variospa(chispa_df, weights = "exp", summary = TRUE)
  c1 <- as.numeric(spa_result[[1]])
  beta1 <- as.numeric(spa_result[[2]])
  alpha1 <- as.numeric(spa_result[[3]])

  # Plot zeta(chi)
  chispa_eta_estim <- ggplot(etachispa_df, aes(lagspa, chi)) +
    btf_theme +
    geom_point(col = btfgreen, size = 4) +
    xlab(TeX(r"($\log(h)$)")) +
    ylab(TeX(r"($\zeta(\widehat{\chi}(h, 0))$)")) +
    theme(axis.ticks = element_blank(),
          axis.line = element_blank(),
          panel.border = element_blank(),
          legend.text = element_text(size = 15),
          legend.title = element_text(size = 15),
          legend.position = "right",
          legend.margin = margin(0.5, 0.5, 0.5, 0, "cm"),
          panel.grid = element_line(color = "#5c595943")) +
    geom_line(aes(x = lagspa, y = alpha1 * lagspa + c1), alpha = 0.6,
              color = "darkred", linewidth = 1.5)

  ggsave(paste0(im_folder, "WLSE/omsev/spatial/full_spatial_chi_zeta_estim_exp_", q_no0_spa, ".pdf"),
         plot = chispa_eta_estim, width = 20, height = 15, units = "cm")
  
  # Save results
  spa_results[[length(spa_results) + 1]] <- data.frame(
    q_spa = q_no0_spa,
    beta1 = beta1,
    alpha1 = alpha1,
    signif_beta1 = spa_result[4],
    signif_alpha1 = spa_result[5]
  )
}

df_spa_result <- bind_rows(spa_results)


# Temporal chi WLSE ============================================================
temp_results <- list()
for (q_no0_temp in q_temp_vals) {
  print(q_no0_temp)
  # Chi(0, tau) boxplot
  chimat_dtlag <- temporal_chi(rain, quantile = q_no0_temp,
                               tmax = tmax, mean = FALSE, zeros = FALSE)
  
  chi_df_dt <- as.data.frame(chimat_dtlag)[, -1]
  colnames(chi_df_dt) <- (1:tmax) * 5 # convert to minutes
  
  df_gathered <- chi_df_dt %>%
    pivot_longer(everything(), names_to = "variable", values_to = "value") %>%
    mutate(group = factor(variable, levels = (1:tmax) * 5))
  
  chitemp <- ggplot(df_gathered, aes(x = group, y = value)) +
    geom_boxplot(fill = btfgreen, alpha = 0.4, outlier.color = btfgreen) +
    btf_boxplot_theme +
    xlab(TeX(r"($\tau$ (minutes))")) +
    ylab(TeX(r"($\widehat{\chi}(0,\tau)$)"))
  
  ggsave(paste0(im_folder, "WLSE/omsev/temporal/full_temporal_chi_boxplot_", q_no0_temp, ".pdf"),
         plot = chitemp, width = 20, height = 15, units = "cm")
  
  # Chi(0, tau) mean
  chimat_dt_mean <- temporal_chi(rain, tmax, quantile = q_no0_temp,
                                 mean = TRUE, zeros = FALSE)
  
  df_chi <- data.frame(lag = (0:tmax), chi = chimat_dt_mean)
  df_chi_not0 <- df_chi[df_chi$lag > 0, ]
  
  temp_result <- get_estimate_variotemp(df_chi_not0, weights = "exp", summary = TRUE)
  
  dftemp <- data.frame(
    lag = log(df_chi_not0$lag),
    chi = eta(df_chi_not0$chi)
  )

  c2 <- as.numeric(temp_result[[1]])
  beta2 <- as.numeric(temp_result[[2]])
  alpha2 <- as.numeric(temp_result[[3]])

  chitemp_eta_estim <- ggplot(dftemp, aes(x = lag, y = chi)) +
    geom_point(color = btfgreen) +
    btf_theme +
    xlab(TeX(r"($\log(\tau)$)")) +
    ylab(TeX(r"($\zeta(\widehat{\chi}(0,\tau))$)")) +
    geom_line(aes(x = lag, y = alpha2 * lag + c2),
              alpha = 0.5, color = "darkred", linewidth = 0.5)
  
  ggsave(paste0(im_folder, "WLSE/omsev/temporal/full_temporal_chi_zeta_estim_", q_no0_temp, ".pdf"),
         plot = chitemp_eta_estim, width = 20, height = 15, units = "cm")
  
  temp_results[[length(temp_results) + 1]] <- data.frame(
    q_temp = q_no0_temp,
    beta2 = temp_result[2],
    alpha2 = temp_result[3],
    signif_beta2 = temp_result[4],
    signif_alpha2 = temp_result[5]
  )
}

df_temp_result <- bind_rows(temp_results)


# Save results to CSV
write.csv(df_spa_result, file = paste0(foldername, "wlse_spa_summary.csv"), row.names = FALSE)
write.csv(df_temp_result, file = paste0(foldername, "wlse_temp_summary.csv"), row.names = FALSE)

# get csv
foldername <- paste0(data_folder, "omsev/WLSE/")
df_spa_result <- read.csv(paste0(foldername, "wlse_spa_summary.csv"))
df_temp_result <- read.csv(paste0(foldername, "wlse_temp_summary.csv"))
df_result_all <- merge(df_spa_result, df_temp_result, by = NULL)

df_result_all <- df_result_all[, c("q_spa", "q_temp",
                                   "beta1", "alpha1", "beta2", "alpha2")]
# put q_spa and q_temp in the first two columns
colnames(df_result_all)

# round values to 4 decimal places
df_result_all$beta1 <- round(df_result_all$beta1, 4)
df_result_all$alpha1 <- round(df_result_all$alpha1, 4)
df_result_all$beta2 <- round(as.numeric(df_result_all$beta2), 4)
df_result_all$alpha2 <- round(as.numeric(df_result_all$alpha2), 4)

df_result_all

# LaTeX
df_temp_result <- df_temp_result[, c("q_temp", "beta2", "alpha2")]
# round values to 4 decimal places
df_temp_result$beta2 <- round(as.numeric(df_temp_result$beta2), 4)
df_temp_result$alpha2 <- round(as.numeric(df_temp_result$alpha2), 4)
kable(df_temp_result, format = "latex", digits = 3) %>%
  kable_styling(latex_options = c("HOLD_position"))


# LaTeX
df_spa_result <- df_spa_result[, c("q_spa", "beta1", "alpha1")]
# round values to 4 decimal places
df_spa_result$beta1 <- round(as.numeric(df_spa_result$beta1), 4)
df_spa_result$alpha1 <- round(as.numeric(df_spa_result$alpha1), 4)
kable(df_spa_result, format = "latex", digits = 3) %>%
  kable_styling(latex_options = c("HOLD_position"))

# choice of one row
df_result <- df_result_all[df_result_all$q_spa == 0.95 & df_result_all$q_temp == 0.95, ]
