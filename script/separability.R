# Check separability assumption

# Load libraries
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

################################################################################
# LOCATION ---------------------------------------------------------------------
################################################################################
# get location of each rain gauge
location_gauges <- read.csv("../data/PluvioMontpellier_1min/pluvio_mtp_loc.csv")
location_gauges$codestation <- c("iem", "mse", "poly", "um", "cefe", "cnrs",
                                 "crbm", "archiw", "archie", "um35", "chu1",
                                 "chu2", "chu3", "chu4", "chu5", "chu6", "chu7")

# Get distances matrix
dist_mat <- get_dist_mat(location_gauges)
df_dist <- reshape_distances(dist_mat)

################################################################################
# DATA -------------------------------------------------------------------------
################################################################################
# get rain measurements
# load data
load("../data/PluvioMontpellier_1min/rain_mtp_5min_2019_2022.RData")
rain <- rain.all5[c(1, 6:ncol(rain.all5))]
# Remove the Time column to focus on site data
rain <- rain[, -1]

# Compute covariance matrix
calculate_spatio_temporal_covariance <- function(data) {
  n_sites <- ncol(data)
#   n_times <- nrow(data)

  # Initialize a covariance matrix
  cov_matrix <- matrix(0, n_sites, n_sites)

  # Calculate covariance between each pair of sites across all time points
  for (i in 1:n_sites) {
    for (j in 1:n_sites) {
      # remove NA for the pair of sites
      data_pair <- na.omit(data.frame(data[[i]], data[[j]]))
      cov_matrix[i, j] <- cov(data_pair[[1]], data_pair[[2]])
    }
  }
  return(cov_matrix)
}

cov_matrix <- calculate_spatio_temporal_covariance(rain)

plot(cov_matrix)

temporal_covariances <- sapply(rain, function(site) {
  site <- na.omit(site)
  cov(site, site)
})

print(temporal_covariances)
plot(temporal_covariances)

transposed_data <- t(rain)

spatial_covariances <- apply(transposed_data, 2, function(instant) {
  instant <- na.omit(instant)
  cov(instant, instant)
})

plot(spatial_covariances)

# Separability asumption
C_h_t <- spatial_covariances + temporal_covariances


library(geosphere)

# Obtenir les distances géodésiques
dist_matrix_geo <- distm(
  location_gauges[, c("Longitude", "Latitude")],
  fun = distHaversine
)

colnames(dist_matrix_geo) <- location_gauges$codestation
rownames(dist_matrix_geo) <- location_gauges$codestation

# Convertir en dataframe
dist_df <- as.data.frame(as.table(dist_matrix_geo))
names(dist_df) <- c("Site1", "Site2", "Distance")

dist_df <- dist_df[as.character(dist_df$Site1) != as.character(dist_df$Site2), ]

breaks <- seq(min(dist_df$Distance), max(dist_df$Distance), length.out = 11)
dist_df$Class <- cut(dist_df$Distance, breaks = breaks, labels = 1:10, include.lowest = TRUE)
table(dist_df$Class)

dist_df$Class <- cut(dist_df$Distance, breaks = quantile(dist_df$Distance, probs = seq(0, 1, length.out = 11)), 
                     labels = 1:10, include.lowest = TRUE)

table(dist_df$Class)

mean_dist <- aggregate(Distance ~ Class, data = dist_df, FUN = mean)

labels <- 1:10

df_lags$class <- cut(lags_h, breaks = breaks, labels = labels, include.lowest = TRUE)

table(df_lags$class)

kmeans_result <- kmeans(h_vect, centers = 10)

df_lags$class <- factor(kmeans_result$cluster)


lags_tau <- 0:10 

results <- expand.grid(h = lags_h, tau = lags_tau)
results$C_h_tau <- mapply(function(h, tau) spatio_temporal_covariance(site_data, h, tau),
                          results$h, results$tau)
results$C_S_h <- sapply(results$h, function(h) spatial_covariance(site_data, h))
results$C_T_tau <- sapply(results$tau, function(tau) temporal_covariance(site_data, tau))
results$Sum_C_S_T <- results$C_S_h + results$C_T_tau

# C(h, tau) ≈ C_S(h) + C_T(tau) ???
# results$Equality <- abs(results$C_h_tau - results$Sum_C_S_T) < 1e-6

print(results)

library(geosphere)

dist_matrix <- distm(
location_gauges[, c("Longitude", "Latitude")],
fun = distHaversine
)

colnames(dist_matrix) <- location_gauges$codestation
rownames(dist_matrix) <- location_gauges$codestation

dist_df <- as.data.frame(as.table(dist_matrix))
names(dist_df) <- c("Site1", "Site2", "Distance")

dist_df <- dist_df[as.character(dist_df$Site1) != as.character(dist_df$Site2), ]

autocorr_pairs <- data.frame(
Site1 = location_gauges$codestation,
Site2 = location_gauges$codestation,
Distance = 0
)
dist_df <- rbind(dist_df, autocorr_pairs)
breaks <- quantile(dist_df$Distance[dist_df$Distance > 0],
                     probs = seq(0, 1, length.out = 11))
dist_df$Class <- cut(dist_df$Distance, breaks = c(0, breaks),
                       labels = 0:10, include.lowest = TRUE)


spatio_temporal_covariance <- function(data, dist_df, tau_values,
                                        num_classes = 10) {
  num_classes <- length(unique(dist_df$Class)) - 1

  class_means <- aggregate(Distance ~ Class, data = dist_df, mean)

  results <- expand.grid(Class = 0:num_classes, Tau = tau_values)
  results$Covariance <- NA
  results$h_mean <- NA

  for (tau in tau_values) {
    for (class in 0:num_classes) {
      pairs <- dist_df[dist_df$Class == class, ]

      if (class == 0) {
        covariances <- sapply(1:nrow(pairs), function(i) {
          site <- as.integer(pairs$Site1[i])

          if (tau >= nrow(data)) return(NA)

          original <- data[1:(nrow(data) - tau), site, drop = TRUE]
          shifted <- data[(1 + tau):nrow(data), site, drop = TRUE]

          safe_cov <- function(x, y) cov(x, y, use = "complete.obs")
          safe_cov(original, shifted)
        })
      } else {
        covariances <- sapply(1:nrow(pairs), function(i) {
          site1 <- as.integer(pairs$Site1[i])
          site2 <- as.integer(pairs$Site2[i])

          if (tau >= nrow(data)) return(NA)

          original <- data[1:(nrow(data) - tau), site1, drop = TRUE]
          shifted <- data[(1 + tau):nrow(data), site2, drop = TRUE]

          safe_cov <- function(x, y) cov(x, y, use = "complete.obs")
          safe_cov(original, shifted)
        })
      }

      mean_cov <- mean(covariances, na.rm = TRUE)
      results$Covariance[results$Class == class & results$Tau == tau] <- mean_cov

      if (class > 0) {
        results$h_mean[results$Class == class & results$Tau == tau] <- class_means$Distance[class_means$Class == class]
      } else {
        results$h_mean[results$Class == class & results$Tau == tau] <- 0
      }
    }
  }
  
  return(results)
}
