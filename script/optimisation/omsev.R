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
rain$dates <- as.POSIXct(rain.all5$dates, tz = "Europe/Paris")
rain <- rain[rain$dates >= "2019-09-06" & rain$dates <= "2024-01-31", ]
rain$dates <- with_tz(rain$dates, tzone = "UTC")
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
rain$mse
# Get distances matrix
dist_mat <- get_dist_mat(location_gauges)
df_dist <- reshape_distances(dist_mat)
x11()
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

################################################################################
# QUANTILE ---------------------------------------------------------------------
################################################################################

sites_names <- colnames(rain)
# site_combinations <- combn(sites_names, 2, simplify = FALSE)

# for (site in sites_names) {
#   filename <- paste0(im_folder, "mrlplot/omsev/", site, ".png")
#   png(filename, width = 10, height = 5, units = "in", res = 300)
#   par(mfrow = c(1, 1))
#   rain_pair <- rain[, c(site)]
#   # remove na
#   rain_pair <- rain_pair[complete.cases(rain_pair), ]
#   rain_pair_no0 <- rain_pair[rowSums(rain_pair) > 0, ]
#   mrlplot(rain_pair)

#   # get corresponding quantile for threshold u

#   dev.off()
# }

# library(ggplot2)


# sites_names <- colnames(rain)
# site_combinations <- combn(sites_names, 2, simplify = FALSE)
# q <- 0.955
# excess_counts_spa <- data.frame(site1 = character(), site2 = character(), 
#                       n_excesses = integer(), stringsAsFactors = FALSE)
# for (pair in site_combinations) {
#   site1 <- pair[1]
#   site2 <- pair[2]
#   filename <- paste0(im_folder, "chiplot/omsev/", site1, "_", site2, ".png")
#   png(filename, width = 10, height = 5, units = "in", res = 300)

#   par(mfrow = c(1, 1))
#   rain_pair <- rain[, c(site1, site2)]
#   # remove na
#   rain_pair <- rain_pair[complete.cases(rain_pair), ]
#   colnames(rain_pair)
#   rain_pair_no0 <- rain_pair[rowSums(rain_pair) > 0, ]
#   chiplot(rain_pair_no0, xlim = c(0.9, 1), ylim1 = c(0, 1), which = 1,
#           qlim = c(0.9, 0.995), main1 = "Without zeros")
#   abline(v = q, col = "red", lty = 2)
#   dev.off()
#   n <- nrow(rain_pair_no0)
#   data_unif <- cbind(rank(rain_pair_no0[, 1]) / (n + 1),
#                     rank(rain_pair_no0[, 2]) / (n + 1))

#   count_excesses <- sum(data_unif[, 1] > q & data_unif[, 2] > q)
#   excess_counts_spa <- rbind(excess_counts_spa, data.frame(site1 = site1,
#                                                    site2 = site2,
#                                                    n_excesses = count_excesses,
#                                                    stringsAsFactors = FALSE))
 
# }

# excess_save <- excess_counts_spa
# # get latex table for excess counts
# excess_counts_spa <- excess_counts_spa %>%
#   mutate(site1 = factor(site1, levels = sites_names),
#          site2 = factor(site2, levels = sites_names)) %>%
#   arrange(site1, site2)


# ggplot(excess_counts_spa %>% filter(n_excesses > 0),
#        aes(x = site2, y = site1, size = n_excesses)) +
#   geom_point(alpha = 0.7, color = btfgreen) +
#   geom_text(aes(label = ifelse(n_excesses < 40, n_excesses, "")),
#             size = 3, vjust = -1) +
#   scale_size_continuous(range = c(1, 10)) +
#   theme_minimal() +
#   labs(title = "", size = "Number of excesses") +
#   theme(axis.text.x = element_text(angle = 45, hjust = 1))
# # save plot
# filename <- paste(im_folder, "WLSE/omsev/spatial/excess_counts_spatial_", 
#                     q, ".png", sep = "")
# ggsave(filename, width = 20, height = 15, units = "cm")




# templag <- 0:10
# q <- 0.94
# excess_counts_temp <- data.frame(site1 = character(), tau = integer(), 
#                       n_excesses = integer(), stringsAsFactors = FALSE)
# for (site in sites_names) {
#   for (tau in templag) {
#     site1 <- site
#     site2 <- paste0(site, "_lag", tau)
#     # create a lagged version of the site
#     rain_lagged <- cbind(rain[[site]][1:(nrow(rain) - tau)],
#                          rain[[site]][(tau + 1):nrow(rain)])
#     colnames(rain_lagged) <- c(site1, site2)
#     rain_pair <- as.data.frame(rain_lagged)
#     # remove na
#     rain_pair <- rain_pair[complete.cases(rain_pair), ]
#     rain_pair_no0 <- rain_pair[rowSums(rain_pair) > 0, ]
#     filename <- paste0(im_folder, "chiplot/omsev/", site1, "_", site2, ".png")
#     png(filename, width = 10, height = 5, units = "in", res = 300)

#     par(mfrow = c(1, 1))
#     rain_nolag <- rain[[site1]][1:(length(rain[[site1]]) - tau)]
#     rain_lag <- rain[[site1]][(tau + 1):length(rain[[site1]])]
#     rain_pair <- cbind(rain_nolag, rain_lag)
#     # remove na
#     rain_pair <- rain_pair[complete.cases(rain_pair), ]
#     colnames(rain_pair)
#     rain_pair_no0 <- rain_pair[rowSums(rain_pair) > 0, ]
#     chiplot(rain_pair_no0, xlim = c(0.9, 1), ylim1 = c(0, 1), which = 1,
#             qlim = c(0.9, 0.98), main1 = "Without zeros")
#     abline(v = q, col = "red", lty = 2)
#     dev.off()
#     n <- nrow(rain_pair_no0)
#     data_unif <- cbind(rank(rain_pair_no0[, 1]) / (n + 1),
#                       rank(rain_pair_no0[, 2]) / (n + 1))

#     count_excesses <- sum(data_unif[, 1] > q & data_unif[, 2] > q)
#     excess_counts_temp <- rbind(excess_counts_temp, data.frame(site = site,
#                                                    tau = tau,
#                                                    n_excesses = count_excesses,
#                                                    stringsAsFactors = FALSE))
#   }
 
# }

# library(ggplot2)
# library(scales)

# ymin <- min(excess_counts_temp$n_excesses, na.rm = TRUE)

# ggplot(excess_counts_temp, aes(x = tau, y = n_excesses, color = site)) +
#   geom_line() +
#   geom_point() +
#   theme_minimal() +
#   scale_y_continuous(
#     limits = c(0, NA),
#     breaks = pretty(c(0, excess_counts_temp$n_excesses), n = 15)
#   ) +
#   scale_x_continuous(
#     breaks = templag
#   ) +
#   labs(
#     title = "",
#     x = "Temporal lag",
#     y = "Number of joint excesses",
#     color = "Site"
#   )


# # save plot
# filename <- paste(im_folder, "WLSE/omsev/temporal/excess_counts_temporal_",
#                     q, ".png", sep = "")
# ggsave(filename, width = 20, height = 15, units = "cm")

################################################################################
# EXTREMOGRAM ------------------------------------------------------------------
################################################################################
# Spatial and temporal CHIPLOT and CHI value for all pairs

# TEMPORAL CHI -----------------------------------------------------------------

# remove stations  cines, hydro, brives
# rain <- rain[, -c(18, 19, 20)] # remove cines, hydro, brives
tmax <- 10
nsites <- ncol(rain) # number of sites
q_temp <- 0.93 # quantile for temporal chi
# Temporal chi with spatial lag fixed at 0
# compute chiplot of every site with itself but lagged in time
start_time <- Sys.time()
chimat_dtlag <- temporal_chi(rain, tmax, quantile = q_temp, mean = FALSE, zeros = FALSE)
end_time <- Sys.time()
elapsed_time <- end_time - start_time
print(elapsed_time)


chi_df_dt <- data.frame(chimat_dtlag)
colnames(chi_df_dt) <- c(0:tmax) * 5 # temporal lags from 1 to tmax
rownames(chi_df_dt) <- c(1:nsites) # stations
# remove lag 0
chi_df_dt <- chi_df_dt[, -1] # remove lag 0

# boxplot all pixels values for chi temp
# Reshape data using gather function
df_gathered <- chi_df_dt %>% gather(key = "variable", value = "value")
df_gathered$group <- factor(as.integer(df_gathered$variable),
                            levels = seq(0, tmax) * 5)

# Plot boxplots
chitemp <- ggplot(df_gathered, aes(x = group, y = value)) +
  geom_boxplot(fill = btfgreen, alpha = 0.4, outlier.color = btfgreen) +
  btf_boxplot_theme +
  xlab(TeX(r"($\tau$ (minutes))")) +
  ylab(TeX(r"($\widehat{\chi}(0,\tau)$)"))

chitemp

filename <- paste(im_folder, "WLSE/omsev/temporal/full_temporal_chi_boxplot_", q_temp,
                 ".pdf", sep = "")

# save plot
ggsave(filename, plot = chitemp, width = 20, height = 15, units = "cm",
       dpi = 600, device = "pdf")


chimat_dt_mean <- temporal_chi(rain, tmax, quantile = q_temp, mean = TRUE,
                               zeros = FALSE)
# get h axis in minutes ie x5 minutes
df <- data.frame(lag = c(0:tmax) * 5, chi = chimat_dt_mean)

# remove lag 0
df <- df[-1, ] # remove lag 0
ggplot(df, aes(x = lag, y = chi)) +
  geom_point(color = btfgreen) +
  btf_theme +
  xlab(TeX(r"($\tau$ (minutes))")) +
  ylab(TeX(r"($\hat{\chi}$)"))


dftemp <- data.frame(lag = log(df$lag), chi = eta(df$chi))

wlse_temp <- get_estimate_variotemp(df,
                                    weights = "residuals", summary = TRUE)
c2 <- as.numeric(wlse_temp[[1]])
beta2 <- as.numeric(wlse_temp[[2]])
alpha2 <- as.numeric(wlse_temp[[3]])
y <- alpha2 * df$lag + c2  #idem

chitemp_eta_estim <- ggplot(dftemp, aes(x = lag, y = chi)) +
  geom_point(color = btfgreen, size = 4) +
  btf_theme +
  xlab(TeX(r"($\log(\tau)$)")) +
  ylab(TeX(r"($\zeta(\widehat{\chi}(0,\tau))$)")) +
  geom_line(aes(x = lag, y = alpha2 * lag + c2),
            alpha = 0.6, color = "darkred", linewidth = 1.5)

chitemp_eta_estim

# save plot
filename <- paste(im_folder, "WLSE/omsev/temporal/full_temporal_chi_eta_estim_residuals_", q_temp,
                 ".pdf", sep = "")

chitemp_eta_estim + theme(
  axis.title = element_text(size = 10),
  axis.text = element_text(size = 8)
)

ggsave(filename, width = 20, height = 15, units = "cm", dpi = 600, device = "pdf")


# SPATIAL CHI ------------------------------------------------------------------

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

# Assigner NA aux distances non classées (>= dernier rayon)
rad_mat[dist_mat >= max(radius)] <- NA
rad_mat[dist_mat == 0] <- NA

# Triangulaire
rad_mat[lower.tri(rad_mat)] <- NA
rad_mat_sym <- rad_mat
rad_mat_sym[lower.tri(rad_mat)] <- t(rad_mat)[lower.tri(rad_mat)]
hmax <- max(radius) # distance max in absolute value...
nb_col <- length(unique(radius)) # we exclude 0 ?

q_spa <- 0.95
# plot chi for each distance
chispa_df <- spatial_chi(rad_mat, rain,
                         quantile = q_spa, zeros = FALSE, mid = TRUE)

etachispa_df <- data.frame(chi = eta(chispa_df$chi),
                           lagspa = log(chispa_df$lagspa))

par(mfrow = c(1, 1))

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

chispa_plot

ymin <- min(etachispa_df$chi)

chispa_eta <- ggplot(etachispa_df, aes(lagspa, chi)) +
  btf_theme +
  geom_point(col = "darkred", size = 4) +
  xlab(TeX(paste0("$ log ||", expression(bold("h")), "||$"))) +
  ylab(TeX(paste0("$\\widehat{\\chi}(", expression(bold("h")), ", 0)$"))) +
  theme(axis.ticks = element_blank(),
        axis.line = element_blank(),
        panel.border = element_blank(),
        legend.text = element_text(size = 15),
        legend.title = element_text(size = 15),
        legend.position = "right",
        legend.margin = margin(0.5, 0.5, 0.5, 0, "cm"),
        panel.grid = element_line(color = "#5c595943"))

chispa_eta


chispa_df$lagspa <- chispa_df$lagspa # convert to km
wlse_spa <- get_estimate_variospa(chispa_df, weights = "residuals", summary = TRUE, bw = 1.5)
c1 <- wlse_spa[[1]]
beta1 <- wlse_spa[[2]]
alpha1 <- wlse_spa[[3]]
y <- alpha1 * etachispa_df$lagspa + c1  #idem
print(paste("c1:", c1, "beta1:", beta1, "alpha1:", alpha1))

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
            color = "darkred", linewidth = 1.5) +
        ylim(-2, 0) # set y limits to ymin and max chi

chispa_eta_estim

# save plot
filename <- paste(im_folder, "WLSE/omsev/spatial/spatial_chi_eta_estim_residuals_",
                 q_spa, "_radius_", length(radius) - 1, ".pdf", sep = "")
ggsave(filename, width = 20, height = 15, units = "cm")


################################################################################
# VARIOGRAM --------------------------------------------------------------------
################################################################################

# estimates
param_omsev <- c(beta1, beta2, alpha1, alpha2)

# in rain remove when all data are NA<
rain <- rain[rowSums(is.na(rain)) < ncol(rain), ]
q <- 0.95 # quantile
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
filename <- paste(im_folder, "optim/omsev/300m_threshold_histogram_q",
                  q * 1000, "_delta", delta, "_dmin", min_spatial_dist,
                  ".png", sep = "")
ggsave(filename, width = 20, height = 15, units = "cm")

n_episodes <- length(selected_points$s0)
print(n_episodes)
length(unique(selected_points$s0)) # can be same s0
length(unique(selected_points$t0)) # never same t0?
print(min(selected_points$u_s0)) # min threshold
t0_list <- selected_points$t0
s0_list <- selected_points$s0

list_episodes_points <- get_extreme_episodes(selected_points, rain,
                                     episode_size = episode_size, unif = FALSE,
                                     beta = 0)

list_episodes <- list_episodes_points$episodes
episode <- list_episodes[[1]]


library(dplyr)
library(lubridate)
library(readr)
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


date <- datetimes[1]
rain_date <- rain[date, ]


library(lubridate)

# Convert to POSIXct if needed
selected_points$t0_date <- as.POSIXct(selected_points$t0_date, format = "%Y-%m-%d %H:%M:%S", tz = "UTC")

# Round to the nearest hour
selected_points$t0_date_rounded <- round_date(selected_points$t0_date, unit = "hour")

datetimes <- unique(selected_points$t0_date_rounded)


# save datetimes of episode in a csv file
# datetimes_df <- data.frame( datetime = selected_points$t0_date,
#                             day = selected_points$day,
#                             month = selected_points$month,
#                             year = selected_points$year,
#                             hour = selected_points$hour,
#                             minute = selected_points$minute,
#                             second = selected_points$second)

datetimes_df <- data.frame( datetime =datetimes)
datetimes_df$datetime <- format(datetimes_df$datetime, "%Y-%m-%d %H:%M:%S")

filename <- paste0(data_folder, "omsev/t0_episodes",
           "_q", q * 100, "_delta", delta,
           "_dmin", min_spatial_dist, ".csv")

write.csv(datetimes_df, file = filename, row.names = FALSE)

datetimes_df <- read_csv(filename, col_types = cols(
  datetime = col_datetime(format = "")
))


# check the episode
head(episode)
library(ggplot2)
library(reshape2)  # for melting wide data to long format

# Convert matrix to data frame
index <- 105
sort(t0_list)
which(t0_list == t0_list[index])
episode_test <- list_episodes[[index]]
df_episode <- as.data.frame(episode_test)
df_episode$Time <- 0:(nrow(df_episode)-1)  # Add a time column
s0_list[index]
u_episode <- selected_points$u_s0[index]
# t0_episode <- t0_list[index]
# Convert from wide to long format
df_long <- melt(df_episode, id.vars = "Time")
head(df_long)
colnames(df_long) <- c("Time", "Pixel", "Value")
ggplot(df_long, aes(x = Time, y = Value, group = Pixel)) +
  geom_line(color = btfgreen) +
  geom_hline(yintercept = u_episode, color = "red", linetype = "dashed") +
  theme_minimal()

filename <- paste(im_folder, "optim/omsev/extreme_episode", index, "_min", min_spatial_dist,
                  "m_max", tmax, "_5min_delta_", delta, ".png", sep = "")
# filename <- "test.png"
ggsave(filename, width = 20, height = 15, units = "cm")


list_episodes_unif_points <- get_extreme_episodes(selected_points, rain,
                                      episode_size = episode_size, unif = TRUE)

list_episodes_unif <- list_episodes_unif_points$episodes

s0_list <- selected_points$s0
t0_list <- selected_points$t0
u_list <- selected_points$u_s0

library(parallel)

tau_vect <- 0:10
thresholds_by_site <- apply(rain, 2, function(col) {
  col <- col[!is.na(col) & col > 0]
  if (length(col) < 30) return(NA)
  quantile(col, probs = q)
})

tmax <- max(tau_vect)
df_coords <- as.data.frame(sites_coords)
# Compute the lags and excesses for each conditional point
list_results <- mclapply(1:length(s0_list), function(i) {
  s0 <- s0_list[i]
  col_s0 <- which(colnames(rain) == s0)
  s0_coords <- df_coords[col_s0, ]
  # t0 <- t0_list[i]
  episode <- list_episodes[[i]]
  ind_t0_ep <- 0 # index of t0 in the episode
  lags <- get_conditional_lag_vectors(df_coords, s0_coords, ind_t0_ep,
                                tau_vect, latlon = TRUE)
  # lags$hnorm <- lags$hnorm / 1000 # convert to km
  # excesses <- empirical_excesses_rpar(episode, q, lags, t0 = ind_t0_ep)
  # Seuils calculés sans les zéros

  excesses <- empirical_excesses_rpar(episode, thresholds = thresholds_by_site,
                                      lags, t0 = ind_t0_ep)

  lags$tau <- lags$tau
  list(lags = lags, excesses = excesses)
}, mc.cores = detectCores() - 1)

list_lags <- lapply(list_results, `[[`, "lags")
list_excesses <- lapply(list_results, `[[`, "excesses")

s0 <- s0_list[1]
col_s0 <- which(colnames(rain) == s0)
s0_coords <- df_coords[col_s0, ]
excesses <- list_excesses[[20]]
sum(excesses$kij)
df_lags <- list_lags[[1]]
tail(df_lags)


# ADD ADV DATA ################################################################

adv_filename <- paste0(data_folder, "omsev/adv_estim/advection_results.csv")
adv_df <- read.csv(adv_filename, sep = ",")
head(adv_df)
nrow(adv_df) # number of advection estimates
# convert adv_df$t0 to POSIXct
adv_df$t0 <- as.POSIXct(adv_df$t0, format = "%Y-%m-%d %H:%M:%S", tz = "UTC")
# convert adv_df
length(adv_df$t0) # should be same as length(t0_list)
length(selected_points$t0_date) # should be same as length(t0_list)
selected_points$adv_x <- adv_df$mean_dx_mps[match(selected_points$t0_date_rounded, adv_df$t0)]
selected_points$adv_y <- adv_df$mean_dy_mps[match(selected_points$t0_date_rounded, adv_df$t0)]


selected_points_nona <- selected_points[!is.na(selected_points$adv_x) & !is.na(selected_points$adv_y), ]

# remove those with 0 advection
# selected_points_adv <- selected_points_nona[
  # selected_points_nona$adv_x != 0 | selected_points_nona$adv_y != 0, ]
selected_points_adv <- selected_points_nona
s0_list <- selected_points_nona$s0
s0_list <- selected_points_adv$s0
list_episodes_points <- get_extreme_episodes(selected_points_adv, rain,
                              episode_size = episode_size, unif = FALSE)

list_episodes <- list_episodes_points$episodes

tmax <- max(tau_vect)
df_coords <- as.data.frame(sites_coords)
# Compute the lags and excesses for each conditional point
list_results <- mclapply(1:length(s0_list), function(i) {
  s0 <- s0_list[i]
  col_s0 <- which(colnames(rain) == s0)
  s0_coords <- df_coords[col_s0, ]
  # t0 <- t0_list[i]
  episode <- list_episodes[[i]]
  ind_t0_ep <- 0 # index of t0 in the episode
  lags <- get_conditional_lag_vectors(df_coords, s0_coords, ind_t0_ep,
                                tau_vect, latlon = TRUE)
  # lags$hnorm <- lags$hnorm / 1000 # convert to km
  # excesses <- empirical_excesses_rpar(episode, q, lags, t0 = ind_t0_ep)
  # Seuils calculés sans les zéros

  excesses <- empirical_excesses_rpar(episode, thresholds = thresholds_by_site,
                                      lags, t0 = ind_t0_ep)

  lags$tau <- lags$tau
  list(lags = lags, excesses = excesses)
}, mc.cores = detectCores() - 1)

list_lags <- lapply(list_results, `[[`, "lags")
list_excesses <- lapply(list_results, `[[`, "excesses")

s0 <- s0_list[1]
col_s0 <- which(colnames(rain) == s0)
s0_coords <- df_coords[col_s0, ]
excesses <- list_excesses[[20]]
sum(excesses$kij)
df_lags <- list_lags[[1]]
tail(df_lags)


# wind_df test
# adv_df <- data.frame(vx = 1, vy = 1)
adv_df <- selected_points_nona[, c("adv_x", "adv_y")]
colnames(adv_df) <- c("vx", "vy")
init_param <- c(beta1, beta2, alpha1, alpha2, 1, 1)
# init_param <- c(0.02, 0.5, alpha1, alpha2, 1, 1)

result <- optim(par = init_param, fn = neg_ll_composite,
        list_lags = list_lags, list_episodes = list_episodes,
        list_excesses = list_excesses, hmax = hmax,
        wind_df = adv_df,
        latlon = TRUE,
        directional = TRUE,
        fixed_eta1 = FALSE,
        fixed_eta2 = TRUE,
        convert_in_hours = FALSE,
        method = "L-BFGS-B",
        lower = c(1e-08, 1e-08, 1e-08, 1e-08, 1e-08, 1e-08),
        upper = c(10, 10, 1.999, 1.999, 10, 10),
        control = list(maxit = 10000,
                      trace = 1,
                      parscale = c(1, 1, 1, 1, 1, 1)),
        hessian = F)

result

#0.64170004 0.23508995 0.07947333 1.06460933 1.01845062 1.00000000


init_param <- c(beta1, beta2, alpha1, alpha2)

result <- optim(par = init_param, fn = neg_ll_composite,
        list_lags = list_lags, list_episodes = list_episodes,
        list_excesses = list_excesses, hmax = sqrt(17),
        latlon = TRUE,
        directional = FALSE,
        convert_in_hours = FALSE,
        method = "L-BFGS-B",
        lower = c(1e-08, 1e-08, 1e-08, 1e-08),
        upper = c(10, 10, 1.999, 1.999),
        control = list(maxit = 10000,
                      trace = 1,
                      parscale = c(1, 1, 1, 1)),
        hessian = F)

result

# [1] 0.01392341 0.79141409 0.47458862 0.59647469
#  Init WLSE: [1] 0.01392341 0.16940296 0.47458862 0.61530967



init_param <- c(beta1, beta2, alpha1, alpha2, 0.1, 0.2)

result <- optim(par = init_param, fn = neg_ll_composite,
        list_lags = list_lags, list_episodes = list_episodes,
        list_excesses = list_excesses, hmax = sqrt(17),
        latlon = TRUE,
        directional = FALSE,
        convert_in_hours = FALSE,
        method = "L-BFGS-B",
        lower = c(1e-08, 1e-08, 1e-08, 1e-08, 1e-08, 1e-08),
        upper = c(10, 10, 1.999, 1.999, 10, 10),
        control = list(maxit = 10000,
                      trace = 1,
                      parscale = c(1, 1, 1, 1, 1, 1)),
        hessian = F)

result

# [1] 0.005393757 0.791221062 0.475214154 0.597500562 0.650280008 0.649669269
#  Init WLSE: [1] 0.01392341 0.16940296 0.47458862 0.61530967 0 0

# [1] 0.6307329 0.7782762 0.4765049 0.5992543 0.1339897 0.2673149
#  Init WLSE: [1] 0.01392341 0.16940296 0.47458862 0.61530967 0.1 0.2


################################################################################


