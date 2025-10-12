# remove objects from the environment
rm(list = ls())

# clear console
cat("\014")

# Load libraries and set theme
source("./script/load_libraries.R")


################################################################################
# LOCATION ---------------------------------------------------------------------
################################################################################
# get location of each rain gauge
filename_loc <- paste0(data_folder, "omsev/pluvio_mtp_loc_till_2022.csv")
location_gauges <- read.csv(filename_loc)
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
# get rain data from omsev
filename_omsev <- paste0(data_folder,
                         "omsev/omsev_5min/rain_mtp_5min_2019_2022.RData")
# filename_omsev <- paste0(data_folder,
#                          "omsev/omsev_5min/rain_mtp_5min_2019_2024_cleaned.csv")

load(filename_omsev)
rain <- rain.all5[c(1, 6:ncol(rain.all5))]
typeof(rain)
class(rain)
colnames(rain)

################################################################################
# QUANTILE ---------------------------------------------------------------------
################################################################################

# get a matrix of high quantiles for all pair
q <- 0.99 # quantile
rownames(rain) <- rain$dates
rain_new <- rain[-1] # remove dates column
# list_count_quant <- quantile_matrix(q, rain_new, qlim = TRUE, remove_zeros = FALSE,
#                                     count_min = 80) # with removing zeros
# quant_mat <- list_count_quant[1][[1]]
# count_mat <- list_count_quant[2]
# quant_mat
# 0.998 equiv to 0.96 without 0

################################################################################
# EXTREMOGRAM ------------------------------------------------------------------
################################################################################
# Spatial and temporal CHIPLOT and CHI value for all pairs

# TEMPORAL CHI -----------------------------------------------------------------

q <- 0.9 # quantile
tmax <- 10
# rain_new <- rain
nsites <- ncol(rain_new)
# Temporal chi with spatial lag fixed at 0
# compute chiplot of every site with itself but lagged in time
start_time <- Sys.time()
chimat_dtlag <- temporal_chi(rain_new, tmax, quantile = q, mean = FALSE, zeros=F)
end_time <- Sys.time()
elapsed_time <- end_time - start_time
print(elapsed_time)

# every chi lagged mean
par(mfrow = c(1, 1))
chi_df_dt <- data.frame(chimat_dtlag)
colnames(chi_df_dt) <- c(0:tmax) # temporal lags from 1 to tmax
rownames(chi_df_dt) <- c(1:nsites) # stations
# remove lag 0 ie column 1
chi_df_dt <- chi_df_dt[-1]
# boxplot all stations values for chi temp
# Reshape data using gather function
df_gathered <- chi_df_dt %>% gather(key = "variable", value = "value")
df_gathered$group <- factor(as.integer(df_gathered$variable),
                            levels = seq(0, tmax))


# Plot boxplots
chitemp <- ggplot(df_gathered, aes(x = group, y = value)) +
  geom_boxplot(fill = btfgreen, alpha = 0.4, outlier.color = btfgreen) +
  btf_boxplot_theme +
  xlab(TeX(r"($\tau$ (minutes))")) +
  ylab(TeX(r"($\widehat{\chi}(0,\tau)$)")) +
  scale_x_discrete(breaks = c(1, 5, 10, 15, 20),
                   labels = c("5", "25", "50", "75", "100")) +
  scale_y_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1),
                   labels = c("0", "0.25", "0.5", "0.75", "1"),
                   limits = c(0, 1)) +
  theme(panel.grid.major = element_line(color = "white"),
        panel.grid.minor = element_line(color = "white"))

chitemp


chimat_dt_mean <- temporal_chi(rain_new, tmax, quantile = q, mean = TRUE, zeros=F)
# remove lag 0
chimat_dt_mean <- chimat_dt_mean[-1]


# plot mean chi
df_chi <- data.frame(lag = c(1:tmax), chi = chimat_dt_mean)
wlse_temp <- get_estimate_variotemp(df_chi, weights = "exp", summary = TRUE)

c2 <- as.numeric(wlse_temp[[1]])
beta2 <- as.numeric(wlse_temp[[2]])
alpha2 <- as.numeric(wlse_temp[[3]])

y <- alpha2 * df$lag + c2  #idem

df <- df_chi
dftemp <- data.frame(lag = log(df$lag), chi = eta(df$chi))

chitemp_eta_estim <- ggplot(dftemp, aes(x = lag, y = chi)) +
  geom_point(color = btfgreen, size = 4) +
  btf_theme +
  xlab(TeX(r"($\log(\tau)$)")) +
  ylab(TeX(r"($\eta(\widehat{\chi}(0,\tau))$)")) +
  geom_line(aes(x = lag, y = alpha2 * lag + c2),
            alpha = 0.6, color = "darkred", linewidth = 1.5)

chitemp_eta_estim


# SPATIAL CHI ------------------------------------------------------------------

# SPATIAL with temporal lag fixed at 0 (ie 1 in the precedent temporal case)

# spatial structure with an almost constant amount of pairs in each intervals
df_dist_order <- df_dist[order(df_dist$value), ]
num_intervals <- 12
quantiles_rad <- quantile(df_dist_order$value,
                            probs = seq(0, 1, length.out = num_intervals + 1))
radius_intervals <- unique(quantiles_rad)
radius <- as.integer(radius_intervals)
radius[length(radius)] <- 1550
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

# plot chi for each distance
chispa_df <- spatial_chi(rad_mat, rain_new,
                         quantile = 0.95, zeros = F)

etachispa_df <- data.frame(chi = eta(chispa_df$chi),
                           lagspa = log(chispa_df$lagspa))


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


wlse_spa <- get_estimate_variospa(chispa_df, weights = "exp", summary = TRUE)
alpha1 <- wlse_spa[[2]]
beta1 <- wlse_spa[[1]]
c1 <- log(beta1)
y <- alpha1 * etachispa_df$lagspa + c1  #idem


chispa_eta_estim <- ggplot(etachispa_df, aes(lagspa, chi)) +
  btf_theme +
  geom_point(col = btfgreen, size = 4) +
  xlab(TeX(r"($\log(h)$)")) +
  ylab(TeX(r"($\eta(\widehat{\chi}(h, 0))$)")) +
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
