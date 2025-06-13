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
rain <- rain.all5[, c(1, 6:(ncol(rain.all5) - 1))]


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

head(rain_clean$mse)

rain <- rain_clean
colnames(rain)
ncol(rain)

rain <- rain[rain$dates >= "2019-12-01",]
# put dates as rownames
rownames(rain) <- rain$dates
rain <- rain[-1] # remove dates column

# Get distances matrix
dist_mat <- get_dist_mat(location_gauges)
df_dist <- reshape_distances(dist_mat)
plot(df_dist$value)
max(df_dist$value)


library(dplyr)

# remove all rows with all NAs
rain <- rain %>%
  filter(!if_all(everything(), is.na))


# remove cines, hydro, brives
rain <- rain[, -c(18, 19, 20)] # remove cines, hydro, brives
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


library(ggplot2)

site1 <- "cines"
site2 <- "hydro"
filename <- paste0(im_folder, "chiplot/omsev/quantile_omsev_", site1, "_", site2, ".png")
png(filename, width = 10, height = 5, units = "in", res = 300)

par(mfrow = c(1, 2))
rain_pair <- rain[, c(site1, site2)]
# remove na
rain_pair <- rain_pair[complete.cases(rain_pair), ]
colnames(rain_pair)
rain_pair_no0 <- rain_pair[rowSums(rain_pair) > 0, ]
chiplot(rain_pair_no0, xlim = c(0.9, 1), ylim1 = c(-0.1, 1), which = 1,
        qlim = c(0.8, 0.997), main1 = "Without zeros")
# abline(v = 0.93, col = "red", lty = 2)

# Same with all zeros
chiplot(rain_pair, xlim = c(0.9995, 1), ylim1 = c(0, 1), which = 1,
        qlim = c(0.9995, 0.99995), main1 = "With zeros")
# abline(v = 0.999875, col = "red", lty = 2)

dev.off()
nrow(rain_pair_no0)
# count conjoint excesses
q <- 0.98
# uniformize the data
n <- nrow(rain_pair_no0)
data_unif <- cbind(rank(rain_pair_no0[, 1]) / (n + 1),
                   rank(rain_pair_no0[, 2]) / (n + 1))

count_excesses <- sum(data_unif[, 1] > q & data_unif[, 2] > q)
print(count_excesses)


q_no0_spa <- 0.94
# get the corresponding quantile in the data without 0 and with 0
q_value_no0 <- quantile(rain_pair_no0[,1], probs = q_no0_spa, na.rm = TRUE)
q_value_no0 <- quantile(rain_pair_no0[,2], probs = q_no0_spa, na.rm = TRUE)


filename <- paste0(im_folder, "chiplot/omsev/quantile_omsev_", site1, "_", site1, ".png")
png(filename, width = 10, height = 5, units = "in", res = 300)

# # Temporal chi
# par(mfrow = c(1,2))
# rain_nolag <- rain[[site1]][1:(length(rain[[site1]]) - 5)]
# rain_lag <- rain[[site1]][6:length(rain[[site1]])]
# rain_pair <- cbind(rain_nolag, rain_lag)
# rain_pair_no0 <- rain_pair[rowSums(rain_pair) > 0, ]
# chiplot(rain_pair_no0, xlim = c(0.9, 1), ylim1 = c(0, 1), which = 1,
#         qlim = c(0.9, 0.995), main1 = "Without zeros")
# abline(v = 0.98, col = "red", lty = 2)

# chiplot(rain_pair, xlim = c(0.99975, 1), ylim1 = c(0, 1), which = 1,
#         qlim = c(0.99975, 0.99998), main1 = "With zeros")
# abline(v = 0.99995, col = "red", lty = 2)

# dev.off()



# q_no0_temp <- 0.98


################################################################################
# EXTREMOGRAM ------------------------------------------------------------------
################################################################################
# Spatial and temporal CHIPLOT and CHI value for all pairs

# TEMPORAL CHI -----------------------------------------------------------------

# remove stations  cines, hydro, brives
# rain <- rain[, -c(18, 19, 20)] # remove cines, hydro, brives
tmax <- 10
nsites <- ncol(rain) # number of sites
q_temp <- 0.998 # quantile for temporal chi
# Temporal chi with spatial lag fixed at 0
# compute chiplot of every site with itself but lagged in time
start_time <- Sys.time()
chimat_dtlag <- temporal_chi(rain, tmax, quantile = q_temp, mean = FALSE, zeros = TRUE)
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

# save plot
filename <- paste(im_folder, "WLSE/omsev/full_temporal_chi_boxplot_", q_temp,
                 ".png", sep = "")
ggsave(filename, width = 20, height = 15, units = "cm")


chimat_dt_mean <- temporal_chi(rain, tmax, quantile = q_temp, mean = TRUE,
                               zeros = TRUE)
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
                                    weights = "exp", summary = TRUE)
c2 <- wlse_temp[[1]]
beta2 <- wlse_temp[[2]]
alpha2 <- wlse_temp[[3]]
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
filename <- paste(im_folder, "WLSE/omsev/full_temporal_chi_eta_estim_", q_temp,
                 ".png", sep = "")
ggsave(filename, width = 20, height = 15, units = "cm")

# SPATIAL CHI ------------------------------------------------------------------

# Trier les distances
# Construire des intervalles fixes
step_km <- 150
max_dist <- max(df_dist$value, na.rm = TRUE)
radius <- seq(0, max_dist + 10, by = step_km)

# Histogramme des paires
dist_counts <- table(cut(df_dist$value, breaks = radius, include.lowest = TRUE))
df_hist <- as.data.frame(dist_counts)
colnames(df_hist) <- c("Interval", "Count")

# Mettre à jour rad_mat
rad_mat <- matrix(NA, nrow = nrow(dist_mat), ncol = ncol(dist_mat))
for (i in seq_len(length(radius) - 1)) {
  lower <- radius[i]
  upper <- radius[i + 1]
  rad_mat[dist_mat > lower & dist_mat <= upper] <- upper
}
rad_mat[dist_mat == 0] <- 0
rad_mat[lower.tri(rad_mat)] <- NA

# SPATIAL with temporal lag fixed at 0 (ie 1 in the precedent temporal case)

# spatial structure with an almost constant amount of pairs in each intervals
df_dist_order <- df_dist[order(df_dist$value), ]
# remove when distance is 0
df_dist_order <- df_dist_order[df_dist_order$value > 0, ]
num_intervals <- 10
quantiles_rad <- quantile(df_dist_order$value,
                            probs = seq(0, 1, length.out = num_intervals + 1))
radius_intervals <- unique(quantiles_rad)
radius <- as.integer(radius_intervals)
radius[length(radius)] <- max(df_dist_order$value) + 10 # last radius is the max distance
dist_counts <- table(cut(df_dist$value, breaks = radius))

# Get dataframe for the histogram plot
df_hist <- data.frame(dist_counts)

colnames(df_hist) <- c("Interval", "Count")

df_hist$Breaks <- gsub("e\\+0.", "0", df_hist$Interval)
df_hist$Breaks <- gsub("\\.", "", df_hist$Breaks)

# # Histogram
histradius <- ggplot(df_hist, aes(x = Interval, y = Count)) +
  btf_theme +
  geom_bar(stat = "identity", fill = btfgreen, alpha = 0.8) +
  xlab("Spatial lag") +
  ylab("Pair count") +
  theme(axis.text.x = element_text(angle = 45)) +
  scale_x_discrete(labels = df_hist$Breaks) +
  scale_y_continuous(breaks = c(0, 4, 6, 8, 10, 12))

histradius

# Create matrix of radius
rad_mat <- dist_mat
# Loop through radius and set distances in matrix
for (i in 2:length(radius)) {
  curr_radius <- radius[i]
  prev_radius <- radius[i - 1]
  rad_mat[dist_mat >= prev_radius & dist_mat < curr_radius] <- curr_radius
  rad_mat[dist_mat > curr_radius] <- Inf
}

rad_mat[dist_mat == 0] <- NA
# Make a triangle
rad_mat[lower.tri(rad_mat)] <- NA
rad_mat_sym <- rad_mat
rad_mat_sym[lower.tri(rad_mat)] <- t(rad_mat)[lower.tri(rad_mat)]

hmax <- max(radius) # distance max in absolute value...
nb_col <- length(unique(radius)) # we exclude 0 ?

q_spa <- 0.998
# plot chi for each distance
chispa_df <- spatial_chi(rad_mat, rain,
                         quantile = q_spa, zeros = TRUE, mid = TRUE)

etachispa_df <- data.frame(chi = eta(chispa_df$chi),
                           lagspa = log(chispa_df$lagspa))

# Exemple :
breaks <- seq(0, 1600, by = 100)
site_names <- colnames(rain)

site_pairs <- cbind(
  site_names[df_dist_order$X],
  site_names[df_dist_order$Y]
)

# distances: vecteur des distances entre les paires
# site_pairs: matrice à deux colonnes (par exemple obtenue par combn(1:S, 2))

chispa_df <- spatial_chi_vector(
  distances = df_dist_order$value,
  site_pairs = site_pairs,
  data_rain = rain,
  quantile = q_spa,
  breaks = radius,     # les bornes de tes intervalles
  zeros = FALSE
)

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
wlse_spa <- get_estimate_variospa(chispa_df, weights = "exp", summary = TRUE)
c1 <- wlse_spa[[1]]
beta1 <- wlse_spa[[2]]
alpha1 <- wlse_spa[[3]]
y <- alpha1 * etachispa_df$lagspa + c1  #idem


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
            color = "darkred", size = 1.5) +
        ylim(-2, 2) # set y limits to ymin and max chi

chispa_eta_estim



################################################################################
# VARIOGRAM --------------------------------------------------------------------
################################################################################

# estimates
param_omsev <- c(beta1, beta2, alpha1, alpha2)



q <- 0.999 # quantile
# get central site from sites_coords

set_st_excess <- get_spatiotemp_excess(rain, q)

list_s <- set_st_excess$list_s
unique(unlist(list_s))
list_t <- set_st_excess$list_t

rain[list_t[[1]], list_s[[1]]]

list_u <- set_st_excess$list_u
list_u[[1]]

# Spatio-temporal neighborhood parameters
min_spatial_dist <- 300 # m
delta <- 10 # in * 5 min
episode_size <- delta # size of the episode
sites_coords <- location_gauges[, c("Longitude", "Latitude")]

s0t0_set <- get_s0t0_pairs(sites_coords, rain,
                            min_spatial_dist = min_spatial_dist,
                            episode_size = episode_size,
                            set_st_excess = set_st_excess,
                            n_max_episodes = 10000,
                            latlon = TRUE)

selected_points <- s0t0_set

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
# check the episode
head(episode)
library(ggplot2)
library(reshape2)  # for melting wide data to long format

# Convert matrix to data frame
index <- 100
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

# filename <- paste(im_folder, "optim/comephore/extreme_episode", index, "_min", min_spatial_dist,
#                   "km_max", tmax, "h_delta_", delta, ".png", sep = "")
# # filename <- "test.png"
# ggsave(filename, width = 20, height = 15, units = "cm")


list_episodes_unif_points <- get_extreme_episodes(selected_points, rain,
                                      episode_size = episode_size, unif = TRUE)

list_episodes_unif <- list_episodes_unif_points$episodes

s0_list <- selected_points$s0
t0_list <- selected_points$t0
u_list <- selected_points$u_s0

library(parallel)

tau_vect <- 0:9
tmax <- max(tau_vect)
df_coords <- as.data.frame(sites_coords)
# Compute the lags and excesses for each conditional point
list_results <- mclapply(1:length(s0_list), function(i) {
  s0 <- s0_list[i]
  col_s0 <- which(colnames(rain) == s0)
  s0_coords <- df_coords[col_s0, ]
  # t0 <- t0_list[i]
  episode <- list_episodes_unif[[i]]
  u <- u_list[i]
  ind_t0_ep <- 0 # index of t0 in the episode
  lags <- get_conditional_lag_vectors(df_coords, s0_coords, ind_t0_ep,
                                tau_vect, latlon = TRUE)
  # lags$hnorm <- lags$hnorm / 1000 # convert to km
  excesses <- empirical_excesses_rpar(episode, q, lags, t0 = ind_t0_ep)
  lags$tau <- lags$tau
  list(lags = lags, excesses = excesses)
}, mc.cores = detectCores() - 1)

list_lags <- lapply(list_results, `[[`, "lags")
list_excesses <- lapply(list_results, `[[`, "excesses")

s0 <- s0_list[index]
col_s0 <- which(colnames(rain) == s0)
s0_coords <- df_coords[col_s0, ]
excesses <- list_excesses[[1]]
excesses$kij
df_lags <- list_lags[[1]]
tail(df_lags)

# ADD WIND DATA ################################################################

list_episodes_points <- get_extreme_episodes(selected_points, rain,
                                episode_size = episode_size, unif = FALSE)

list_episodes <- list_episodes_points$episodes

# # get wind data
filename_wind <- paste0(data_folder, "wind/data_gouv/wind_mtp.csv")
wind_mtp <- read.csv(filename_wind)

# Convert datetime to POSIXct
wind_mtp$datetime <- as.POSIXct(wind_mtp$datetime,
                                format = "%Y-%m-%d %H:%M:%S", tz = "UTC")

# Apply function to the DD column
wind_mtp$cardDir <- sapply(wind_mtp$DD, convert_to_cardinal)
wind_mtp$cardDir <- as.character(wind_mtp$cardDir)  # Ensure it's character
wind_mtp$cardDir[is.na(wind_mtp$DD)] <- NA
# Check if NA values are properly handled
summary(wind_mtp)

head(wind_mtp$cardDir)

# wind_df test
adv_df <- data.frame(vx = 0.1, vy = 0.1)
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