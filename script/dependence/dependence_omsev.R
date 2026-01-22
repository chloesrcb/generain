# remove objects from the environment
rm(list = ls())

# clear console
cat("\014")

# Load libraries and set theme
source("./script/load_libraries.R")
functions_folder <- "./R"
files <- list.files(functions_folder, full.names = TRUE)
# load all functions in files
invisible(lapply(files, function(f) source(f, echo = FALSE)))
library(latex2exp)

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
# filename_omsev <- paste0(data_folder,
#                          "omsev/omsev_5min/rain_mtp_5min_2019_2022.RData")
# load(filename_omsev)
# rain <- rain.all5[c(1, 6:ncol(rain.all5))]
filename_omsev <- paste0(data_folder,
                         "omsev/omsev_5min/rain_mtp_5min_2019_2024_cleaned.csv")

rain <- read.csv(filename_omsev)
head(rain)

typeof(rain)
class(rain)
colnames(rain)
rownames(rain) <- rain$dates
rain_new <- rain[-1] # remove dates column

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
dftemp <- data.frame(lag = log(df$lag), chi = zeta(df$chi))

chitemp_eta_estim <- ggplot(dftemp, aes(x = lag, y = chi)) +
  geom_point(color = btfgreen, size = 4) +
  btf_theme +
  xlab(TeX(r"($\log(\tau)$)")) +
  ylab(TeX(r"($\zeta(\widehat{\chi}(0,\tau))$)")) +
  geom_line(aes(x = lag, y = alpha2 * lag + c2),
            alpha = 0.6, color = "darkred", linewidth = 1.5)

chitemp_eta_estim


# SPATIAL CHI ------------------------------------------------------------------

# SPATIAL with temporal lag fixed at 0 (ie 1 in the precedent temporal case)

# spatial structure with an almost constant amount of pairs in each intervals
df_dist_order <- df_dist[order(df_dist$value), ]
num_intervals <- 11
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

# plot chi for each distance
chispa_df <- spatial_chi(rad_mat, rain_new,
                         quantile = 0.95, zeros = F)

etachispa_df <- data.frame(chi = zeta(chispa_df$chi),
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

# remove lagspa 0
chispa_df <- chispa_df[chispa_df$lagspa!=0,]
wlse_spa <- get_estimate_variospa(chispa_df, weights = "exp", summary = TRUE)
alpha1 <- wlse_spa[[3]]
beta1 <- wlse_spa[[2]]
c1 <- wlse_spa[[1]]
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
            color = "darkred", size = 1.5)

chispa_eta_estim

#####################################################
# Spatio-temporal chi and variogram


# Ex : on veut tous les couples (s1<s2) et taus de 0 à tmax
build_df_lags <- function(dist_mat, rad_mat, tmax) {
  n <- ncol(dist_mat)
  pairs <- which(upper.tri(dist_mat), arr.ind = TRUE)
  df_pairs <- data.frame(s1 = pairs[,1], s2 = pairs[,2],
                         dist = dist_mat[upper.tri(dist_mat)],
                         lagspa = rad_mat[upper.tri(dist_mat)])
  taus <- 0:tmax
  # expand pour tous les taus
  df_lags <- tidyr::crossing(df_pairs, tau = taus)
  # ordre facultatif
  df_lags <- df_lags[order(df_lags$lagspa, df_lags$tau), ]
  rownames(df_lags) <- NULL
  return(df_lags)
}

                      
tmax <- 10    # ou ton choix
df_lags <- build_df_lags(dist_mat, rad_mat, tmax)

# df_lags <- get_lag_vectors(location_gauges, tau_vect = 0:10, latlon=TRUE)
# df_lags$hnorm <- df_lags$hnorm * 1000
q <- 0.95
chi_st <- chispatemp_empirical(data_rain = rain_new, df_lags = df_lags,
                               quantile = q, remove_zeros = TRUE)

head(chi_st)

# applique à toute la table
chi_st$vario_emp2 <- vario_spatemp(chi_st$chiemp2)


library(dplyr)

# Agréger par lag spatial et tau
vario_by_h_tau <- chi_st %>%
  group_by(lagspa, tau) %>%
  summarise(vario = mean(vario_emp2, na.rm = TRUE),
            n_pairs = sum(!is.na(vario_emp2)),
            .groups = "drop")

# Voir
head(vario_by_h_tau)


library(ggplot2)

taus_plot <- 0:10
df_plot <- vario_by_h_tau %>% filter(tau %in% taus_plot)

# Tracer tous les variogrammes, un par tau
ggplot(df_plot, aes(x = lagspa, y = vario, color = factor(tau))) +
  geom_point() +
  geom_line() +
  labs(x = "Distance",
       y = "Variogram",
       color = "Temporal lag") +
  theme_minimal()


library(ggplot2)
library(dplyr)

df_plot <- vario_by_h_tau %>% filter(tau %in% taus_plot)

ggplot(df_plot, aes(x = lagspa, y = vario, color = factor(tau))) +
  geom_point() +
  geom_smooth(se = FALSE, method = "loess", span = 1) +   # lissage local
  labs(x = "Distance",
       y = "Variogram",
       color = "Temporal lag") +
  theme_minimal()


foldername_fig <- paste0(im_folder, "variogram/empirical/omsev/variogram_separability/")
ggsave(paste0(foldername_fig, "variogram_spatemp_omsev_loess_q95_alltaus.png"),
       width = 10, height = 6)


ggplot(df_plot, aes(x = lagspa, y = vario, color = factor(tau))) +
  geom_point() +
  geom_line() +   # lignes droites entre points
  labs(x = "Distance",
       y = "Variogram",
       color = "Temporal lag") +
  theme_minimal()


foldername_fig <- paste0(im_folder, "variogram/empirical/omsev/variogram_separability/")
ggsave(paste0(foldername_fig, "variogram_spatemp_omsev_raw_q95.png"),
       width = 10, height = 6)

