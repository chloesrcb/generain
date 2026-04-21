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
filename_omsev <- paste0(data_folder,
                         "omsev/omsev_5min/rain_mtp_5min_2019_2022.RData")
load(filename_omsev)
rain <- rain.all5[c(1, 6:ncol(rain.all5))]
# filename_omsev <- paste0(data_folder,
#                          "omsev/omsev_5min/rain_mtp_5min_2019_2025.csv")
# rain <- rain.all_save
# rain <- read.csv(filename_omsev)
head(rain)

typeof(rain)
class(rain)
colnames(rain)
rownames(rain) <- rain$dates
rain_new <- rain[-1] # remove dates column

# remove brives, cines and hydro columns
# rain_new <- rain_new[, -which(colnames(rain_new) %in% c("brives", "cines", "hydro"))]

################################################################################
# EXTREMOGRAM ------------------------------------------------------------------
################################################################################
# Spatial and temporal CHIPLOT and CHI value for all pairs

# TEMPORAL CHI -----------------------------------------------------------------

q <- 0.95 # quantile
tmax <- 10
tmax_min <- tmax * 5 # in minutes
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
                   limits = c(0, 1))
                

chitemp


foldername_fig <- paste0(im_folder, "WLSE/omsev/temporal/new/")
if (!dir.exists(foldername_fig)) {
  dir.create(foldername_fig, recursive = TRUE)
}
# save
ggsave(paste0(foldername_fig, "chitemp_boxplot_omsev_q", q*100, ".pdf"),
       width = 8, height = 6)

chimat_dt_mean <- temporal_chi(rain_new, tmax, quantile = q, mean = TRUE, zeros=F)
# remove lag 0
chimat_dt_mean <- chimat_dt_mean[-1]


# plot mean chi
tmax_hour <- tmax # convert to hours
df_chi <- data.frame(lag = c(1:tmax)*5/60, chi = chimat_dt_mean)
wlse_temp <- get_estimate_variotemp(df_chi, weights = "residuals")

c2 <- as.numeric(wlse_temp$estimate[[1]])
beta2 <- as.numeric(wlse_temp$estimate[[2]])
alpha2 <- as.numeric(wlse_temp$estimate[[3]])

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

ggsave(paste0(foldername_fig, "chitemp_zeta_estim_omsev_q", q*100, ".pdf"),
       width = 8, height = 6)

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
  xlab("Spatial lag (m)") +
  ylab("Pair count") +
  theme(axis.text.x = element_text(angle = 45)) +
  scale_x_discrete(labels = df_hist$Breaks) +
  scale_y_continuous(breaks = c(0, 4, 6, 8, 10, 12))

histradius

# SAVE 
foldername_fig <- paste0(im_folder, "WLSE/omsev/spatial/new/")
if (!dir.exists(foldername_fig)) {
  dir.create(foldername_fig, recursive = TRUE)
}
ggsave(paste0(foldername_fig, "hist_radius_omsev_q", q*100, ".pdf"),
       width = 8, height = 6)
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

rad_mat <- rad_mat  / 1000 # convert to km
# plot chi for each distance
chispa_df <- spatial_chi(rad_mat, rain_new,
                         quantile = 0.95, zeros = F)

etachispa_df <- data.frame(chi = zeta(chispa_df$chi),
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

# save
foldername_fig <- paste0(im_folder, "WLSE/omsev/spatial/new/")
if (!dir.exists(foldername_fig)) {
  dir.create(foldername_fig, recursive = TRUE)
}
ggsave(paste0(foldername_fig, "chispa_omsev_q", q*100, ".pdf"),
       width = 8, height = 6)


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
wlse_spa <- get_estimate_variospa(chispa_df, weights = "exp")
alpha1 <- wlse_spa$estimate[[3]]
beta1 <- wlse_spa$estimate[[2]]
c1 <- wlse_spa$estimate[[1]]


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

# save 
ggsave(paste0(foldername_fig, "chispa_eta_estim_omsev_q", q*100, ".pdf"),
       width = 8, height = 6)


# result spatial and temporal in a table
results_table <- data.frame(
  Parameter = c("beta1", "alpha1", "beta2", "alpha2"),
  Estimate = c(beta1, alpha1, beta2, alpha2)
)

#####################################################
# Spatio-temporal chi and variogram

build_df_lags <- function(dist_mat, rad_mat, tmax) {
  n <- ncol(dist_mat)
  pairs <- which(upper.tri(dist_mat), arr.ind = TRUE)
  df_pairs <- data.frame(s1 = pairs[,1], s2 = pairs[,2],
                         dist = dist_mat[upper.tri(dist_mat)],
                         lagspa = rad_mat[upper.tri(dist_mat)])
  taus <- 0:tmax
  df_lags <- tidyr::crossing(df_pairs, tau = taus)
  df_lags <- df_lags[order(df_lags$lagspa, df_lags$tau), ]
  rownames(df_lags) <- NULL
  return(df_lags)
}

tmax <- 10
df_lags <- build_df_lags(dist_mat, rad_mat, tmax)

# convert tau in hours
# df_lags$tau <- df_lags$tau * 5 / 60
q <- 0.9
chi_st <- chispatemp_empirical(data_rain = rain_new, df_lags = df_lags,
                               quantile = q, remove_zeros = TRUE)

head(chi_st)

chi_st$vario_emp2 <- vario_spatemp(chi_st$chiemp2)


library(dplyr)

vario_by_h_tau <- chi_st %>%
  group_by(lagspa, tau) %>%
  summarise(vario = mean(vario_emp2, na.rm = TRUE),
            n_pairs = sum(!is.na(vario_emp2)),
            .groups = "drop")

head(vario_by_h_tau)

taus_plot <- 0:10
df_plot <- vario_by_h_tau %>% filter(tau %in% taus_plot)

ggplot(df_plot, aes(x = lagspa, y = vario, color = factor(tau))) +
  geom_point() +
  geom_smooth(se = FALSE, method = "loess", span = 2) +
  labs(x = "Distance (m)",
       y = "Variogram",
       color = "Temporal lag (5 min)") +
  theme_minimal() +
  btf_theme


foldername_fig <- paste0(im_folder, "variogram/empirical/omsev/variogram_separability/")
ggsave(paste0(foldername_fig, "variogram_spatemp_omsev_loess_q95_alltaus.png"),
       width = 10, height = 6)


ggplot(df_plot, aes(x = lagspa, y = vario, color = factor(tau))) +
  geom_point() +
  geom_line() +   # lignes droites entre points
  labs(x = "Distance",
       y = "Variogram",
       color = "Temporal lag") +
  theme_minimal() +
  btf_theme


foldername_fig <- paste0(im_folder, "variogram/empirical/omsev/variogram_separability/")
ggsave(paste0(foldername_fig, "variogram_spatemp_omsev_raw_q95.png"),
       width = 10, height = 6)





library(dplyr)
library(tidyr)
library(ggplot2)
library(tibble)
library(dplyr)
library(tidyr)
library(ggplot2)
library(tibble)

# ============================================================================
# Données
# ============================================================================
df <- tribble(
  ~dataset, ~unit,      ~weights,       ~beta1, ~beta1_low, ~beta1_high, ~alpha1, ~alpha1_low, ~alpha1_high, ~beta2, ~beta2_low, ~beta2_high, ~alpha2, ~alpha2_low, ~alpha2_high,
  "OMSEV",  "m/5 min",  "Residuals",    0.006,  0.001,      0.053,       0.501,   0.161,       0.840,        0.545,  0.415,      0.716,       0.576,   0.425,       0.726,
  "OMSEV",  "m/5 min",  "Exponential",  0.006,  0.001,      0.053,       0.501,   0.161,       0.840,        0.423,  0.391,      0.457,       0.846,   0.749,       0.942,
  "OMSEV",  "km/h",     "Residuals",    0.191,  0.152,      0.241,       0.501,   0.161,       0.840,        2.280,  1.969,      2.640,       0.576,   0.425,       0.726,
  "OMSEV",  "km/h",     "Exponential",  0.191,  0.152,      0.241,       0.501,   0.161,       0.840,        2.022,  1.878,      2.178,       0.353,   0.232,       0.474
)

# ============================================================================
# Mise en forme longue
# ============================================================================
df_long <- bind_rows(
  df %>% transmute(unit, weights, parameter = "hat(beta)[1]",  estimate = beta1,  low = beta1_low,  high = beta1_high),
  df %>% transmute(unit, weights, parameter = "hat(alpha)[1]", estimate = alpha1, low = alpha1_low, high = alpha1_high),
  df %>% transmute(unit, weights, parameter = "hat(beta)[2]",  estimate = beta2,  low = beta2_low,  high = beta2_high),
  df %>% transmute(unit, weights, parameter = "hat(alpha)[2]", estimate = alpha2, low = alpha2_low, high = alpha2_high)
) %>%
  mutate(
    unit = factor(unit, levels = c("m/5 min", "km/h")),
    weights = factor(weights, levels = c("Residuals", "Exponential")),
    parameter = factor(
      parameter,
      levels = c("hat(beta)[1]", "hat(alpha)[1]", "hat(beta)[2]", "hat(alpha)[2]")
    )
  )

# ============================================================================
# Fonction de plot
# ============================================================================
plot_wlse_by_weight <- function(data, chosen_weight) {

  data_sub <- data %>% filter(weights == chosen_weight)

  ggplot(data_sub, aes(x = unit, y = estimate, group = 1)) +
    geom_line(linewidth = 0.6, colour = "grey60") +
    geom_errorbar(
      aes(ymin = low, ymax = high),
      width = 0.06,
      linewidth = 0.7
    ) +
    geom_point(size = 2.8) +
    facet_wrap(~ parameter, scales = "free_y", nrow = 1, labeller = label_parsed) +
    labs(
      x = NULL,
      y = "Estimate",
      title = paste("OMSEV –", chosen_weight, "weights")
    ) +
    theme_bw(base_size = 12) +
    theme(
      legend.position = "none",
      strip.background = element_rect(fill = "grey95", colour = "grey70"),
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_blank(),
      plot.title = element_text(face = "bold")
    ) +
    btf_theme
}

# ============================================================================
# Figures
# ============================================================================
p_resid <- plot_wlse_by_weight(df_long, "Residuals")
p_exp   <- plot_wlse_by_weight(df_long, "Exponential")

p_resid
p_exp


library(dplyr)
library(ggplot2)
library(tibble)

df <- tribble(
  ~unit,      ~weights,       ~beta1, ~beta1_low, ~beta1_high, ~alpha1, ~alpha1_low, ~alpha1_high, ~beta2, ~beta2_low, ~beta2_high, ~alpha2, ~alpha2_low, ~alpha2_high,
  "m/5 min",  "Residuals",    0.006,  0.001,      0.053,       0.501,   0.161,       0.840,        0.545,  0.415,      0.716,       0.576,   0.425,       0.726,
  "m/5 min",  "Exponential",  0.006,  0.001,      0.053,       0.501,   0.161,       0.840,        0.423,  0.391,      0.457,       0.846,   0.749,       0.942,
  "km/h",     "Residuals",    0.191,  0.152,      0.241,       0.501,   0.161,       0.840,        2.280,  1.969,      2.640,       0.576,   0.425,       0.726,
  "km/h",     "Exponential",  0.191,  0.152,      0.241,       0.501,   0.161,       0.840,        2.022,  1.878,      2.178,       0.353,   0.232,       0.474
)

# ============================================================================
# Spatial
# ============================================================================
df_spatial <- bind_rows(
  df %>% transmute(unit, weights, parameter = "hat(beta)[1]",  estimate = beta1,  low = beta1_low,  high = beta1_high),
  df %>% transmute(unit, weights, parameter = "hat(alpha)[1]", estimate = alpha1, low = alpha1_low, high = alpha1_high)
) %>%
  mutate(
    unit = factor(unit, levels = c("m/5 min", "km/h")),
    weights = factor(weights, levels = c("Residuals", "Exponential")),
    parameter = factor(parameter, levels = c("hat(beta)[1]", "hat(alpha)[1]"))
  )

p_spatial <- ggplot(df_spatial, aes(x = estimate, y = parameter, colour = unit, shape = unit)) +
  geom_errorbar(
    aes(xmin = low, xmax = high),
    orientation = "y",
    width = 0.10,
    linewidth = 0.7,
    position = position_dodge(width = 0.35)
  ) +
   geom_text(
    aes(label = sprintf("%.3f", estimate)),
    position = position_dodge(width = 0.35),
    vjust = -1.3,
    size = 5
  ) +
  geom_point(
    size = 3,
    position = position_dodge(width = 0.35)
  ) +
  facet_wrap(~ weights, scales = "free_x", nrow = 1) +
  scale_y_discrete(labels = c(
    "hat(beta)[1]"  = expression(hat(beta)[1]),
    "hat(alpha)[1]" = expression(hat(alpha)[1])
  )) +
  scale_colour_manual(values = c("m/5 min" = "#1b9e77", "km/h" = "#d95f02")) +
  labs(
    x = "Estimate with 95% confidence interval",
    y = NULL,
    colour = "Space-time unit",
    shape = "Space-time unit",
    title = ""
  ) +
  theme_bw(base_size = 15) +
  theme(
    legend.position = "top",
    strip.background = element_blank(),
    panel.grid = element_line(color = "#5c595943"),
    panel.border = element_rect(color = "grey70")
  ) +
  theme(
    strip.text = element_text(size = 14, face = "bold")
  )

p_spatial

# save 
foldername_fig <- paste0(im_folder, "WLSE/omsev/spatial/new/")
if (!dir.exists(foldername_fig)) {
  dir.create(foldername_fig, recursive = TRUE)
}
ggsave(paste0(foldername_fig, "wlse_spatial_omsev.pdf"),
       width = 8, height = 6)
# ============================================================================
# Temporal
# ============================================================================
df_temporal <- bind_rows(
  df %>% transmute(unit, weights, parameter = "hat(beta)[2]",  estimate = beta2,  low = beta2_low,  high = beta2_high),
  df %>% transmute(unit, weights, parameter = "hat(alpha)[2]", estimate = alpha2, low = alpha2_low, high = alpha2_high)
) %>%
  mutate(
    unit = factor(unit, levels = c("m/5 min", "km/h")),
    weights = factor(weights, levels = c("Residuals", "Exponential")),
    parameter = factor(parameter, levels = c("hat(beta)[2]", "hat(alpha)[2]"))
  )

p_temporal <- ggplot(df_temporal, aes(x = estimate, y = parameter, colour = unit, shape = unit)) +
  geom_errorbar(
    aes(xmin = low, xmax = high),
    orientation = "y",
    width = 0.10,
    linewidth = 0.7,
    position = position_dodge(width = 0.35)
  ) +
  geom_point(
    size = 3,
    position = position_dodge(width = 0.35)
  ) +
  geom_text(
    aes(label = sprintf("%.3f", estimate)),
    position = position_dodge(width = 0.35),
    vjust = -1.3,
    size = 5
  ) +
  facet_wrap(~ weights, scales = "free_x", nrow = 1) +
  scale_y_discrete(labels = c(
    "hat(beta)[2]"  = expression(hat(beta)[2]),
    "hat(alpha)[2]" = expression(hat(alpha)[2])
  )) +
  scale_colour_manual(values = c("m/5 min" = "#1b9e77", "km/h" = "#d95f02")) +
  labs(
    x = "Estimate with 95% confidence interval",
    y = NULL,
    colour = "Space-time unit",
    shape = "Space-time unit",
    title = ""
  ) +
  theme_bw(base_size = 15) +
  theme(
    legend.position = "top",
    strip.background = element_blank(),
    panel.grid = element_line(color = "#5c595943"),
    panel.border = element_rect(color = "grey70")
  ) +
  theme(
    strip.text = element_text(size = 14, face = "bold")
  )

p_temporal

# save
foldername_fig <- paste0(im_folder, "WLSE/omsev/temporal/new/")
if (!dir.exists(foldername_fig)) {
  dir.create(foldername_fig, recursive = TRUE)
}

ggsave(paste0(foldername_fig, "wlse_temporal_omsev.pdf"),
       width = 8, height = 6)







library(ggplot2)
library(dplyr)
library(tibble)

df_comephore <- tribble(
  ~parameter,        ~estimate, ~low,   ~high,
  "hat(beta)[1]",    0.016,     0.016,  0.017,
  "hat(alpha)[1]",   1.429,     1.405,  1.453,
  "hat(beta)[2]",    1.553,     1.507,  1.600,
  "hat(alpha)[2]",   0.283,     0.256,  0.310
)

df_comephore <- df_comephore %>%
  mutate(
    parameter = factor(
      parameter,
      levels = c("hat(beta)[1]", "hat(alpha)[1]", "hat(beta)[2]", "hat(alpha)[2]")
    ),
    label_val = sprintf("%.3f", estimate)
  )

p <- ggplot(df_comephore, aes(x = parameter, y = estimate)) +
  geom_errorbar(
    aes(ymin = low, ymax = high),
    width = 0.12,
    linewidth = 0.8
  ) +
  geom_point(size = 3) +
  geom_text(
    aes(label = label_val),
    vjust = -1.1,
    hjust = 1.5,
    size = 5
  ) +
  scale_x_discrete(labels = function(x) parse(text = x)) +
  scale_y_continuous(
    expand = expansion(mult = c(0.05, 0.12))
  ) +
  labs(
    x = NULL,
    y = "Estimate (km/h)",
    title = ""
  ) +
  theme_bw(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold"),
    panel.grid.minor = element_blank()
  ) +
  btf_theme

p

# save
foldername_fig <- paste0(im_folder, "WLSE/comephore/")
if (!dir.exists(foldername_fig)) {
  dir.create(foldername_fig, recursive = TRUE)
}
ggsave(paste0(foldername_fig, "wlse_comephore_parameters.pdf"),
       width = 8, height = 6)



#####
#' Calcul du Chi Spatio-Temporel Empirique
#' @param data_rain Matrix ou DF (lignes = temps, colonnes = stations)
#' @param dist_mat Matrice des distances entre stations
#' @param tau_max Lag temporel maximum (ex: 10 pour 10 * 5min)
#' @param q Quantile de seuil (ex: 0.95)
#' 
calc_chi_spatemp <- function(data_rain, dist_mat, tau_max, q) {
  library(dplyr)
  library(tidyr)

  n_sites <- ncol(data_rain)
  n_time <- nrow(data_rain)
  
  # 1. Transformation en marges uniformes (Rangs) par station
  # On traite les NA en gardant leur position
  data_unif <- apply(data_rain, 2, function(x) {
    r <- rank(x, na.last = "keep", ties.method = "random")
    return(r / (sum(!is.na(x)) + 1))
  })

  # 2. Préparation du stockage des résultats
  results <- expand.grid(
    s1 = 1:n_sites,
    s2 = 1:n_sites,
    tau = 0:tau_max
  )
  
  # Ajout de la distance physique pour analyse ultérieure
  results$distance <- apply(results, 1, function(row) dist_mat[row["s1"], row["s2"]])
  
  # 3. Calcul du Chi empirique pour chaque combinaison
  # Chi(h, tau) = P( U(s2, t+tau) > q | U(s1, t) > q )
  results$chi <- apply(results, 1, function(row) {
    s1 <- row["s1"]
    s2 <- row["s2"]
    tau <- row["tau"]
    
    # Séries temporelles décalées
    # U1 : station s1 de t = 1 à T - tau
    # U2 : station s2 de t = 1 + tau à T
    u1 <- data_unif[1:(n_time - tau), s1]
    u2 <- data_unif[(1 + tau):n_time, s2]
    
    # Suppression des couples contenant des NA
    mask <- !is.na(u1) & !is.na(u2)
    u1 <- u1[mask]
    u2 <- u2[mask]
    
    if(length(u1) == 0) return(NA)
    
    # Nombre d'excès marginaux (Dénominateur)
    n_exc_s1 <- sum(u1 > q)
    
    if(n_exc_s1 == 0) return(0)
    
    # Nombre d'excès joints (Numérateur)
    n_joint_exc <- sum(u1 > q & u2 > q)
    
    return(n_joint_exc / n_exc_s1)
  })
  
  return(results)
}


# --- Exemple d'application ---
# 1. Calculer le chi brut
df_chi <- calc_chi_spatemp(rain_new, dist_mat, tau_max = 10, q = 0.95)

# 2. Créer des classes de distance (ex: tous les 500m)
df_chi$dist_class <- cut(df_chi$distance, breaks = seq(0, max(dist_mat), by = 500))

# 3. Agréger (Moyenne du Chi par distance et par Lag temporel)
chi_summary <- df_chi %>%
  group_by(dist_class, tau) %>%
  summarise(
    mean_dist = mean(distance, na.rm = TRUE),
    mean_chi = mean(chi, na.rm = TRUE),
    n_pairs = n(),
    .groups = "drop"
  )

# 4. Visualisation
ggplot(chi_summary, aes(x = mean_dist, y = mean_chi, color = as.factor(tau))) +
  geom_line() +
  geom_point() +
  labs(title = "Chi Spatio-Temporel Empirique",
       x = "Distance (m)", y = expression(chi(h, tau)),
       color = "Lag temporel (5 min)") +
  theme_minimal()