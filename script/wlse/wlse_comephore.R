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

library(sf)
library(dplyr)
library(ggplot2)
library(gridExtra)

# load("workspace.RData")
# library(generain)

# LOAD DATA ####################################################################
filename_com <- paste0(data_folder, "comephore/comephore_2008_2024_10km.csv")
comephore_raw <- read.csv(filename_com, sep = ",")
# filename_loc <- paste0(data_folder, "comephore/coords_pixels_10km.csv")
filename_loc <- paste0(data_folder, "comephore/loc_px_zoom_10km.csv")
loc_px <- read.csv(filename_loc, sep = ",")
# colnames(comephore_raw)
head(comephore_raw)
# colnames(comephore_raw)[1] <- "date" # rename date column

# remove pixel in loc_px that are not in comephore_raw
loc_px <- loc_px[loc_px$pixel_name %in% colnames(comephore_raw)[-1], ]
nrow(loc_px) # number of pixels
# reindex
rownames(loc_px) <- NULL
# length(colnames(comephore_raw)) - 1
# length(loc_px$pixel_name)
# remove columns date 
df_comephore <- as.data.frame(comephore_raw)
df_comephore$date <- as.POSIXct(df_comephore$date,
                                format = "%Y-%m-%d %H:%M:%S", tz = "UTC")
# when date is just a date, put 00:00:00 as time
# if (nchar(df_comephore$date[1]) == 10) {
#   df_comephore$date <- paste(df_comephore$date, "00:00:00")
# }


head(df_comephore)

length(unique(df_comephore$date))


# Take only data after 2007
# colnames(df_comephore)[1] <- "date"
df_comephore <- df_comephore[df_comephore$date >= "2008-01-01", ]
head(df_comephore)

# show duplicates dates
# duplicated_dates <- df_comephore[duplicated(df_comephore$date), ]
# put date in index
rownames(df_comephore) <- format(as.POSIXct(df_comephore$date), "%Y-%m-%d %H:%M:%S")
comephore <- df_comephore[-1] # remove dates column

# DISTANCE AND COORDS ##########################################################

# Get distances matrix
# dist_mat <- get_dist_mat(loc_px)

# Get number of sites
nsites <- nrow(loc_px)

# Get coords
sites_coords <- loc_px[, c("Longitude", "Latitude")]
rownames(sites_coords) <- loc_px$pixel_name
sites_coords_sf <- st_as_sf(sites_coords, coords = c("Longitude", "Latitude"),
                            crs = 4326)

sites_coords_sf <- st_transform(sites_coords_sf, crs = 2154)

coords_m <- st_coordinates(sites_coords_sf)

grid_coords_km <- sites_coords
grid_coords_m <- sites_coords
grid_coords_m$x_m <- (coords_m[, "X"] - min(coords_m[, "X"]))
grid_coords_m$y_m <- (coords_m[, "Y"] - min(coords_m[, "Y"]))
grid_coords_km$x_km <- (coords_m[, "X"] - min(coords_m[, "X"])) / 1000
grid_coords_km$y_km <- (coords_m[, "Y"] - min(coords_m[, "Y"]))  / 1000

sites_coords_sf <- st_transform(sites_coords_sf, crs = 2154)
coords_m <- st_coordinates(sites_coords_sf)
grid_coords_km <- as.data.frame(coords_m / 1000)
colnames(grid_coords_km) <- c("Longitude", "Latitude")
rownames(grid_coords_km) <- rownames(sites_coords)

dist_mat <- get_dist_mat(grid_coords_km, latlon = FALSE)

# Spatial chi
df_dist_km <- reshape_distances(dist_mat)
# df_dist$value <- round(df_dist$value / 1000, 1) * 1000  # / 1000 in km
# df_dist_km <- df_dist
df_dist_km$value <- round(df_dist_km$value, 1)


# Spatial chi WLSE #############################################################

library(knitr)
library(kableExtra)

# quantiles
q_spa_vals <- c(0.90, 0.91, 0.92, 0.93, 0.94, 0.95, 0.96, 0.97, 0.98, 0.99)
q_temp_vals <- c(0.90, 0.91, 0.92, 0.93, 0.94, 0.95, 0.96, 0.97)
hmax <- 10
tmax <- 10

foldername <- paste0(data_folder, "/comephore/WLSE/")
dir.create(foldername, recursive = TRUE, showWarnings = FALSE)

spa_results <- list()
temp_results <- list()

# Spatial loop =================================================================
for (q_no0_spa in q_spa_vals) {

  chispa_df <- spatial_chi_alldist(df_dist_km, data_rain = comephore,
                                   quantile = q_no0_spa, hmax = hmax, 
                                   zeros = FALSE)

  etachispa_df <- data.frame(
    chi = zeta(chispa_df$chi),
    lagspa = log(chispa_df$lagspa)
  )

  # Plot chi(h, 0)
  chispa_plot <- ggplot(chispa_df, aes(lagspa, chi)) +
    btf_theme +
    geom_point(col = btfgreen) +
    xlab(TeX(r"($h$)")) +
    ylab(TeX(r"($\widehat{\chi}(h, 0)$)")) +
    ylim(0, 1)
  
  ggsave(paste0(im_folder, "WLSE/comephore/full_spatial_chi_", q_no0_spa, ".pdf"),
         plot = chispa_plot, width = 20, height = 15, units = "cm")

  # Estimation WLSE
  spa_result <- get_estimate_variospa(chispa_df, weights = "exp", summary = TRUE)
  alpha <- as.numeric(spa_result[[3]])
  beta <- as.numeric(spa_result[[2]])
  c <- as.numeric(spa_result[[1]])
  # Plot zeta(chi)
  chispa_eta_estim <- ggplot(etachispa_df, aes(lagspa, chi)) +
    btf_theme +
    geom_point(col = btfgreen) +
    xlab(TeX(r"($\log(h)$)")) +
    ylab(TeX(r"($\zeta(\widehat{\chi}(h, 0))$)")) +
    geom_line(aes(x = lagspa, y = alpha * lagspa + c),
              alpha = 0.6, color = "darkred", linewidth = 0.5)
  
  ggsave(paste0(im_folder, "WLSE/comephore/full_spatial_chi_zeta_estim_", q_no0_spa, ".pdf"),
         plot = chispa_eta_estim, width = 20, height = 15, units = "cm")
  
  # Save results
  spa_results[[length(spa_results) + 1]] <- data.frame(
    q_spa = q_no0_spa,
    beta1 = beta,
    alpha1 = alpha,
    signif_beta1 = spa_result[[4]],
    signif_alpha1 = spa_result[[5]]
  )
}

df_spa_result <- bind_rows(spa_results)

# beta2 <- 0.02717288
# alpha2 <- 1.065663

# Temporal chi WLSE ============================================================
for (q_no0_temp in q_temp_vals) {
  print(q_no0_temp)
  # Chi(0, tau) boxplot
  chimat_dtlag <- temporal_chi(comephore, quantile = q_no0_temp,
                               tmax = tmax, mean = FALSE, zeros = FALSE)
  
  chi_df_dt <- as.data.frame(chimat_dtlag)[, -1]
  colnames(chi_df_dt) <- 1:tmax
  
  df_gathered <- chi_df_dt %>%
    pivot_longer(everything(), names_to = "variable", values_to = "value") %>%
    mutate(group = factor(variable, levels = 1:tmax))
  
  chitemp <- ggplot(df_gathered, aes(x = group, y = value)) +
    geom_boxplot(fill = btfgreen, alpha = 0.4, outlier.color = btfgreen) +
    btf_boxplot_theme +
    xlab(TeX(r"($\tau$ (hours))")) +
    ylab(TeX(r"($\widehat{\chi}(0,\tau)$)"))
  
  ggsave(paste0(im_folder, "WLSE/comephore/full_temporal_chi_boxplot_", q_no0_temp, ".pdf"),
         plot = chitemp, width = 20, height = 15, units = "cm")
  
  # Chi(0, tau) mean
  chimat_dt_mean <- temporal_chi(comephore, tmax, quantile = q_no0_temp,
                                 mean = TRUE, zeros = FALSE)
  
  df_chi <- data.frame(lag = 0:tmax, chi = chimat_dt_mean)
  df_chi_not0 <- df_chi[df_chi$lag > 0, ]
  
  temp_result <- get_estimate_variotemp(df_chi_not0, weights = "exp", summary = TRUE)
  
  dftemp <- data.frame(
    lag = log(df_chi_not0$lag),
    chi = zeta(df_chi_not0$chi)
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
  
  ggsave(paste0(im_folder, "WLSE/comephore/full_temporal_chi_zeta_estim_", q_no0_temp, ".pdf"),
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

# choice of one row
df_result <- df_result_all[df_result_all$q_spa == 0.95 & df_result_all$q_temp == 0.95, ]
