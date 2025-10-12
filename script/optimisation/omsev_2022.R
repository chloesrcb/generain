# remove objects from the environment
rm(list = ls())

# clear console
cat("\014")

# Load libraries and set theme
source("./script/load_libraries.R")

# get rain data from omsev
filename_omsev <- paste0(data_folder,
                         "omsev/omsev_5min/rain_mtp_5min_2019_2022.RData")

load(filename_omsev)

# get location of each rain gauge
filename_loc <- paste0(data_folder, "omsev/pluvio_mtp_loc_till_2022.csv")
location_gauges <- read.csv(filename_loc)
location_gauges$Station <- c("iem", "mse", "poly", "um", "cefe", "cnrs",
                             "crbm", "archiw", "archie", "um35", "chu1",
                             "chu2", "chu3", "chu4", "chu5", "chu6", "chu7")


# rain.all5 is the data frame name for 5 min data
# get only stations records and dates
rain <- rain.all5[, c(1, 6:(ncol(rain.all5)))]
colnames(rain)
ncol(rain)
# put dates as rownames
rownames(rain) <- rain$dates
rain <- rain[-1] # remove dates column

# Get distances matrix
dist_mat <- get_dist_mat(location_gauges)
df_dist <- reshape_distances(dist_mat)
plot(df_dist$value)
max(df_dist$value)



################################################################################
# QUANTILE ---------------------------------------------------------------------
################################################################################


library(ggplot2)

sites_names <- colnames(rain)
site_combinations <- combn(sites_names, 2, simplify = FALSE)

for (pair in site_combinations) {
  site1 <- pair[1]
  site2 <- pair[2]
  filename <- paste0(im_folder, "mrlplot/omsev/", site1, "_", site2, ".png")
  png(filename, width = 10, height = 5, units = "in", res = 300)

  par(mfrow = c(1, 1))
  rain_pair <- rain[, c(site1, site2)]
  # remove na
  rain_pair <- rain_pair[complete.cases(rain_pair), ]
  rain_pair_no0 <- rain_pair[rowSums(rain_pair) > 0, ]
  mrlplot(rain_pair_no0)

  # get corresponding quantile for threshold u

  dev.off()
}

u <- 6
abline(v = u, col = "red", lty = 2)
quantile_value <- ecdf(rain_pair_no0[[site1]])(u)
quantile_value <- mean(rain_pair_no0[[site1]] <= u)
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

# Temporal chi
par(mfrow = c(1,2))
rain_nolag <- rain[[site1]][1:(length(rain[[site1]]) - 5)]
rain_lag <- rain[[site1]][6:length(rain[[site1]])]
rain_pair <- cbind(rain_nolag, rain_lag)
rain_pair_no0 <- rain_pair[rowSums(rain_pair) > 0, ]
chiplot(rain_pair_no0, xlim = c(0.9, 1), ylim1 = c(0, 1), which = 1,
        qlim = c(0.9, 0.995), main1 = "Without zeros")
abline(v = 0.98, col = "red", lty = 2)

chiplot(rain_pair, xlim = c(0.99975, 1), ylim1 = c(0, 1), which = 1,
        qlim = c(0.99975, 0.99998), main1 = "With zeros")
abline(v = 0.99995, col = "red", lty = 2)

dev.off()



q_temp <- 0.98




# get a matrix of high quantiles for all pair
# q <- 0.96 # quantile
# list_count_quant <- quantile_matrix(q, rain, qlim = TRUE, remove_zeros = TRUE,
#                                     count_min = 30) # with removing zeros
# quant_mat <- list_count_quant[1][[1]]
# count_mat <- list_count_quant[2]

################################################################################
# EXTREMOGRAM ------------------------------------------------------------------
################################################################################
# Spatial and temporal CHIPLOT and CHI value for all pairs

# TEMPORAL CHI -----------------------------------------------------------------

# remove stations  cines, hydro, brives
# rain <- rain[, -c(18, 19, 20)] # remove cines, hydro, brives
tmax <- 10
nsites <- ncol(rain) # number of sites
# remove when there is na every column
# rain 
q_temp <- 0.95 # quantile for temporal chi
# Temporal chi with spatial lag fixed at 0
# compute chiplot of every site with itself but lagged in time
start_time <- Sys.time()
chimat_dtlag <- temporal_chi(rain, tmax, quantile = q_temp, mean = FALSE, zeros = FALSE)
end_time <- Sys.time()
elapsed_time <- end_time - start_time
print(elapsed_time)


chi_df_dt <- data.frame(chimat_dtlag)
colnames(chi_df_dt) <- c(0:tmax) # temporal lags from 1 to tmax
rownames(chi_df_dt) <- c(1:nsites) # stations
# remove lag 0
chi_df_dt <- chi_df_dt[, -1] # remove lag 0

# boxplot all pixels values for chi temp
# Reshape data using gather function
df_gathered <- chi_df_dt %>% gather(key = "variable", value = "value")
df_gathered$group <- factor(as.integer(df_gathered$variable),
                            levels = seq(0, tmax))

# Plot boxplots
chitemp <- ggplot(df_gathered, aes(x = group, y = value)) +
  geom_boxplot(fill = btfgreen, alpha = 0.4, outlier.color = btfgreen) +
  btf_boxplot_theme +
  xlab(TeX(r"($\tau$ (5 minutes))")) +
  ylab(TeX(r"($\widehat{\chi}(0,\tau)$)"))

chitemp

# save plot
filename <- paste(im_folder, "WLSE/omsev/full_temporal_chi_boxplot_", q_temp,
                 ".png", sep = "")
ggsave(filename, width = 20, height = 15, units = "cm")


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



# remove first lag
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

# SPATIAL with temporal lag fixed at 0 (ie 1 in the precedent temporal case)

# spatial structure with an almost constant amount of pairs in each intervals
df_dist_order <- df_dist[order(df_dist$value), ]
num_intervals <- 12
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

rad_mat[dist_mat == 0] <- 0
# Make a triangle
rad_mat[lower.tri(rad_mat)] <- NA

hmax <- max(radius) # distance max in absolute value...
nb_col <- length(unique(radius)) # we exclude 0 ?
chimat_dslag <- matrix(1, nrow = n - 1, ncol = nb_col)

q_spa <- 0.998
# plot chi for each distance
chispa_df <- spatial_chi(rad_mat, rain,
                         quantile = q_spa, zeros = T)



par(mfrow = c(1, 1))

chispa_plot <- ggplot(chispa_df, aes(lagspa, chi)) +
  btf_theme +
  geom_point(col = btfgreen, size = 4) +
  xlab(TeX(paste0("$||", expression(bold("h")), "||$", " (km)"))) +
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

chispa_df$lagspa <- chispa_df$lagspa / 1000  # in km
etachispa_df <- data.frame(chi = eta(chispa_df$chi),
                           lagspa = log(chispa_df$lagspa))

wlse_spa <- get_estimate_variospa(chispa_df, weights = "exp", summary = TRUE)
c1 <- wlse_spa[[1]]
beta1 <- wlse_spa[[2]]
alpha1 <- wlse_spa[[3]]

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
            color = "darkred", size = 1.5)

chispa_eta_estim

################################################################################
# VARIOGRAM --------------------------------------------------------------------
################################################################################

# estimates
param_ohsm <- c(beta1, beta2, alpha1, alpha2)