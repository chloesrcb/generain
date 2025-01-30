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
location_gauges <- read.csv("./data/PluvioMontpellier_1min/pluvio_mtp_loc.csv")
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
load("./data/PluvioMontpellier_1min/rain_mtp_5min_2019_2022.RData")
rain <- rain.all5[c(1, 6:ncol(rain.all5))]

# put dates as index 
rownames(rain) <- rain$dates
# Remove the Time column to focus on site data
rain <- rain[, -1]

# Remove rows where all values are NA
rain <- rain[rowSums(is.na(rain)) != ncol(rain), ]
head(rain)

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

# from breaks get intervals for each class
class_intervals <- data.frame(
  Class = 0:10,
  Interval = paste0("(", round(breaks, 0), ", ", round(c(breaks[-1], Inf), 0), "]")
  )


# Calculate the frequency of distances in each class
class_frequency <- table(dist_df$Class)

# Create a data frame with class intervals, frequency, and the width for each class
class_interval_summary <- data.frame(
  Class = as.integer(names(class_frequency)),
  Frequency = as.integer(class_frequency),
  IntervalStart = breaks[class_intervals$Class + 1],  # Start of each interval
  IntervalEnd = c(breaks[-1], Inf),  # End of each interval
  IntervalWidth = c(diff(breaks), Inf)  # Width of each interval
)


# Now create the data for geom_rect
df_for_rect <- data.frame(
  xmin = class_interval_summary$IntervalStart,  # Starting position (x-axis)
  xmax = class_interval_summary$IntervalEnd,    # Ending position (x-axis)
  ymin = 0,                                     # Bottom of the bars (y-axis)
  ymax = class_interval_summary$Frequency,       # Height of the bars (y-axis)
  Class = class_interval_summary$Class,          # Class for coloring
  IntervalWidth = class_interval_summary$IntervalWidth  # For sizing
)

# Load ggplot2 for visualization
library(ggplot2)

# Plot the histogram using geom_rect for variable bar widths
ggplot(df_for_rect) +
  geom_rect(aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax), colour = "#f3eeee", fill=btfgreen, alpha=0.5) +
  theme_minimal() +
  labs(
    x = "Distance Class Interval",
    y = "Frequency (number of pairs)"
  ) +
  scale_x_continuous(breaks = breaks, labels = round(breaks, 0)) +  # Custom x-axis breaks
  btf_theme

# save plot
ggsave("./images/separability/histogram_distance_class_intervals.png", width = 6, height = 4)


# Plot
ggplot(data, aes(ymin = 0)) + 
    geom_rect(aes(xmin = left, xmax = right, ymax = value, colour = group, fill = group)) +
    xlab("number of obs") + 
    ylab("value") +
    theme_ipsum() +
    theme(legend.position="none") 

spatio_temporal_covariance <- function(data, dist_df, tau_values, num_classes = 10) {
  num_classes <- length(unique(dist_df$Class)) - 1

  class_means <- aggregate(Distance ~ Class, data = dist_df, mean)

  results <- expand.grid(Class = 0:num_classes, Tau = tau_values)
  results$Covariance <- NA
  results$h_mean <- NA

  for (tau in tau_values) {
    for (class in 0:num_classes) {
      pairs <- dist_df[dist_df$Class == class, ]
      
      covariances <- sapply(1:nrow(pairs), function(i) {
        if (class == 0) {
          # For autocovariance (same site)
          site <- as.integer(pairs$Site1[i])

          if (tau >= nrow(data)) return(NA)

          original <- data[1:(nrow(data) - tau), site, drop = TRUE]
          shifted <- data[(1 + tau):nrow(data), site, drop = TRUE]

        } else {
          # For cross-covariance (different sites)
          site1 <- as.integer(pairs$Site1[i])
          site2 <- as.integer(pairs$Site2[i])

          if (tau >= nrow(data)) return(NA)

          original <- data[1:(nrow(data) - tau), site1, drop = TRUE]
          shifted <- data[(1 + tau):nrow(data), site2, drop = TRUE]
        }

        # Only calculate covariance for non-zero pairs
        valid_indices <- !(original == 0 | shifted == 0)
        safe_cov <- function(x, y) cov(x[valid_indices], y[valid_indices], use = "complete.obs")
        safe_cov(original, shifted)
      })

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


spatial_covariance <- function(data, dist_df) {
  num_classes <- length(unique(dist_df$Class)) - 1
  
  class_means <- aggregate(Distance ~ Class, data = dist_df, mean)
  
  results <- data.frame(Class = 0:num_classes, covspa = NA, h_mean = NA)
  
  for (class in 0:num_classes) {
    pairs <- dist_df[dist_df$Class == class, ]
    
    covariances <- sapply(1:nrow(pairs), function(i) {
      if (class == 0) {
        # Compute covariance for the same site (i.e., autocovariance)
        site <- as.integer(pairs$Site1[i])
        original <- data[, site]
        shifted <- data[, site]
        
      } else {
        # Compute covariance between different sites
        site1 <- as.integer(pairs$Site1[i])
        site2 <- as.integer(pairs$Site2[i])
        original <- data[, site1]
        shifted <- data[, site2]
      }
      
      # Only calculate covariance for non-zero pairs
      valid_indices <- !(original == 0 | shifted == 0)
      safe_cov <- function(x, y) cov(x[valid_indices], y[valid_indices], use = "complete.obs")
      safe_cov(original, shifted)
    })
    
    mean_cov <- mean(covariances, na.rm = TRUE)
    results$covspa[results$Class == class] <- mean_cov
    
    if (class > 0) {
      results$h_mean[results$Class == class] <- class_means$Distance[class_means$Class == class]
    } else {
      results$h_mean[results$Class == class] <- 0
    }
  }
  
  return(results)
}

temporal_covariance <- function(data, tau_values) {
  results <- data.frame(Tau = tau_values, covtemp = NA)
  
  for (tau in tau_values) {
    covariances <- sapply(1:ncol(data), function(site) {
      if (tau >= nrow(data)) return(NA)
      
      original <- data[1:(nrow(data) - tau), site, drop = TRUE]
      shifted <- data[(1 + tau):nrow(data), site, drop = TRUE]
      
      # Only calculate covariance for non-zero pairs
      valid_indices <- !(original == 0 | shifted == 0)
      safe_cov <- function(x, y) cov(x[valid_indices], y[valid_indices], use = "complete.obs")
      safe_cov(original, shifted)
    })
    
    # Compute the mean covariance for the given tau
    mean_cov <- mean(covariances, na.rm = TRUE)
    results$covtemp[results$Tau == tau] <- mean_cov
  }
  
  return(results)
}


spatial_results <- spatial_covariance(rain, dist_df)
ggplot(spatial_results, aes(x = h_mean, y = covspa)) +
  geom_point() + geom_line() +
  ggtitle("Estimated Spatial Covariance C_S(h)")


tau_values <- 0:10  # Define the range of temporal lags to compute covariance
temporal_results <- temporal_covariance(rain, tau_values)

# Plot temporal covariance function C_T(τ)
ggplot(temporal_results, aes(x = Tau, y = covtemp)) +
  geom_point() + geom_line() +
  ggtitle("Estimated Temporal Covariance C_T(τ)")


C_plus <- spatial_results$covspa + temporal_results$covtemp
# spatial_results$C_plus <- C_plus

# # Plot the sum of spatial and temporal covariances
# ggplot(spatial_results, aes(x = h_mean, y = C_plus)) +
#   geom_point() + geom_line() +
#   ggtitle("Estimated Sum of Spatial and Temporal Covariances C_S(h) + C_T(τ)")

tau_values <- 0:10
results <- spatio_temporal_covariance(rain, dist_df, tau_values)

results


df_res <- results

# New dataframe with C(h, tau), C_S + C_T, and the class intervals and tau
df_res <- merge(df_res, class_intervals, by = "Class")
df_res <- merge(df_res, temporal_results, by = "Tau")
df_res <- merge(df_res, spatial_results, by = "Class")

# remove h_mean.y column
df_res <- df_res[, -which(names(df_res) == "h_mean.y")]

# rename h_mean.x to h_mean
colnames(df_res)[colnames(df_res) == "h_mean.x"] <- "h_mean"

df_res$covsep <- df_res$covspa + df_res$covtemp

# reorder columns by class and tau
df_res <- df_res[order(df_res$Tau), ]



# same against tau for fixed h
class <- 8
df_res_h <- df_res[df_res$Class == class, ]
interval <- unique(df_res_h$Interval)
ggplot(df_res_h, aes(x = Tau)) +
  geom_line(aes(y = Covariance, color = "C(h,tau)"), size = 1) +  # Plot C(h,tau)
  geom_point(aes(y = Covariance, color = "C(h,tau)"), size = 2) +  # Plot C(h,tau)
  geom_line(aes(y = covsep, color = "C_s(h) + C_t(tau)"), size = 1, linetype = "dashed") +  # Plot covsep (C_s(h) + C_t(tau))
  geom_point(aes(y = covsep, color = "C_s(h) + C_t(tau)"), size = 2) +  # Plot covsep points
  labs(x = expression("Time lag" ~ tau ~ "(5 minutes)"),  # LaTeX for x-axis
       y = expression("Covariance"),  # LaTeX for y-axis
       color = "Covariance Type") +
  scale_color_manual(values = c("C(h,tau)" = btfgreen, "C_s(h) + C_t(tau)" = "#b68080"), labels = c(TeX(r"($C(h, \tau)$)"), TeX(r"($C(h,0) + C(0,\tau)$)"))) +
  theme_minimal() +
  ggtitle(paste("Covariance for Class = ", interval)) # LaTeX for title

# save plot
ggsave(paste0("./images/separability/covariances_against_tau_class_", class, ".png"), width = 6, height = 4)

# for a fixed tau, plot C(h, tau) against h
tau <- 6
df_res_tau <- df_res[df_res$Tau == tau, ]

ggplot(df_res_tau, aes(x = h_mean)) +
  geom_line(aes(y = Covariance, color = "C(h,tau)"), size = 1) +  # Plot C(h,tau)
  geom_point(aes(y = Covariance, color = "C(h,tau)"), size = 2) +  # Plot C(h,tau)
  geom_line(aes(y = covsep, color = "C_s(h) + C_t(tau)"), size = 1, linetype = "dashed") +  # Plot covsep (C_s(h) + C_t(tau))
  geom_point(aes(y = covsep, color = "C_s(h) + C_t(tau)"), size = 2) +  # Plot covsep points
  labs(x = expression("Distance" ~ h ~ "(meters)"),  # LaTeX for x-axis
       y = expression("Covariance"),  # LaTeX for y-axis
       color = "Covariance Type") +
  scale_color_manual(values = c("C(h,tau)" = btfgreen, "C_s(h) + C_t(tau)" = "#b68080"), labels = c(TeX(r"($C(h, \tau)$)"), TeX(r"($C(h,0) + C(0,\tau)$)"))) +
  theme_minimal() +
  ggtitle(expression(paste("Covariance for Time Lag = ", tau, " minutes")))  # LaTeX for title

# save plot
ggsave(paste0("./images/separability/covariances_against_h_tau_", tau, ".png"), width = 6, height = 4)

# Rename the class labels with class intervals
results <- merge(results, class_intervals, by = "Class")

# Plot the results C(h, tau) vs. time lag
ggplot(results, aes(x = Tau, y = Covariance, color = factor(Interval))) +
  geom_line() +
  geom_point() +
  labs(x = "Time lag (minutes)",
       y = "Covariance",
       color = "Distance class") +
  theme_minimal()

# Ensure Interval is a factor with the correct order
results$Interval <- factor(results$Interval, levels = unique(results$Interval), ordered = TRUE)

# Plot with correct legend order
ggplot(results, aes(x = Tau, y = Covariance, color = Interval)) +
  geom_line() +
  geom_point() +
  labs(x = "Time lag (minutes)",
       y = "Covariance",
       color = "Distance class") +
  scale_color_manual(values = brewer.pal(n = 11, name = "Paired")) +
  btf_theme

# Save plot
ggsave("./images/separability/covariances_against_tau.png", width = 6, height = 4)

library(RColorBrewer)
# Plot the results C(h, tau) vs. distance
ggplot(results, aes(x = h_mean, y = Covariance, color = factor(Tau))) +
  geom_line() +
  geom_point() +
  labs(x = "Distance (meters)",
       y = "Covariance",
       color = "Time lag (5 minutes)") +
  scale_color_manual(values = brewer.pal(n = 11, name = "Paired")) +
  btf_theme

ggsave("./images/separability/covariances_against_h.png", width = 6, height = 4)


C_h <- aggregate(Covariance ~ Class, data = results, FUN = mean)
colnames(C_h) <- c("Class", "C_h")

C_tau <- aggregate(Covariance ~ Tau, data = results, FUN = mean)
colnames(C_tau) <- c("Tau", "C_tau")

# Plot the results C(h) vs. distance
ggplot(C_h, aes(x = Class, y = C_h)) +
  geom_line() +
  geom_point() +
  labs(x = "Distance class",
       y = "Covariance") +
  btf_theme

# Plot the results C(tau) vs. time lag
ggplot(C_tau, aes(x = Tau, y = C_tau)) +
  geom_line() +
  geom_point() +
  labs(x = "Time lag (minutes)",
       y = "Covariance") +
  btf_theme


# Merge C_h and C_tau into results
results <- merge(results, C_h, by = "Class")
results <- merge(results, C_tau, by = "Tau")

# Compute R(h, τ)
results$R_h_tau <- results$Covariance / (results$C_h  + results$C_tau)

ggplot(results, aes(x = h_mean, y = R_h_tau, color = factor(Tau))) +
  geom_line() +
  geom_point() +
  labs(x = "Distances (meters)",
       y = expression(R(h, tau)),
       color = "Time lag (5 minutes)" ) +
  theme_minimal()
