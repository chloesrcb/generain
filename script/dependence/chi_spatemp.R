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
                         "omsev/omsev_5min/rain_mtp_5min_2019_2024.RData")
load(filename_omsev)

rain_raw <- rain.all_save
# keep only CUM_5 == TRUE
# rain_raw <- rain_raw[rain_raw$CUM_5 == TRUE, ]
# rain <- rain_raw[c(1, 6:ncol(rain_raw))]
# remove 
filename_omsev <- paste0(data_folder,
                         "omsev/omsev_5min/rain_mtp_5min_2019_2025.csv")
rain <- read.csv(filename_omsev)
head(rain)

typeof(rain)
class(rain)
colnames(rain)
rownames(rain) <- rain$dates
rain_new <- rain[-1] # remove dates column

# remove brives, cines and hydro sites 
stations_to_remove <- c("brives", "cines", "hydro")
rain_new <- rain_new %>% select(-all_of(stations_to_remove))

calc_chi_spatemp <- function(data_rain, dist_mat, tau_max, q, remove_zeros = TRUE) {
  library(dplyr)
  library(tidyr)

  data_rain <- as.data.frame(data_rain)
  n_sites <- ncol(data_rain)
  n_time  <- nrow(data_rain)

  # keep only unique station pairs
  pairs <- which(upper.tri(dist_mat), arr.ind = TRUE)

  results <- tidyr::crossing(
    data.frame(s1 = pairs[,1], s2 = pairs[,2]),
    tau = 0:tau_max
  )

  results$distance <- mapply(function(i, j) dist_mat[i, j], results$s1, results$s2)
  results$chi <- NA_real_
  results$n_used <- NA_integer_

  for (k in seq_len(nrow(results))) {
    s1  <- results$s1[k]
    s2  <- results$s2[k]
    tau <- results$tau[k]

    # raw aligned series first
    x1 <- data_rain[[s1]][1:(n_time - tau)]
    x2 <- data_rain[[s2]][(1 + tau):n_time]

    # remove NA after alignment
    keep <- !is.na(x1) & !is.na(x2)

    # optionally remove dry-dry pairs
    if (remove_zeros) {
      keep <- keep & (x1 != 0 | x2 != 0)
    }

    x1 <- x1[keep]
    x2 <- x2[keep]

    results$n_used[k] <- length(x1)

    if (length(x1) < 2) {
      results$chi[k] <- NA_real_
      next
    }

    # rank transform on aligned filtered sample
    u1 <- rank(x1, ties.method = "average") / (length(x1) + 1)
    u2 <- rank(x2, ties.method = "average") / (length(x2) + 1)

    n_exc_s1 <- sum(u1 > q)

    if (n_exc_s1 == 0) {
      results$chi[k] <- NA_real_
      next
    }

    n_joint_exc <- sum(u1 > q & u2 > q)
    results$chi[k] <- n_joint_exc / n_exc_s1
  }

  results
}

q <- 0.9995

df_chi <- calc_chi_spatemp(rain_new, dist_mat, tau_max = 10, q = q, remove_zeros = FALSE)

make_radius_quantile <- function(dist_mat, n_bins = 10) {
  
  d <- dist_mat[upper.tri(dist_mat)]
  
  # quantile-based breaks
  breaks <- quantile(d, probs = seq(0, 1, length.out = n_bins + 1), na.rm = TRUE)
  
  # remove duplicates (important)
  breaks <- unique(breaks)
  
  return(breaks)
}

unique_h_values <- sort(unique(df_chi$distance))
radius <- make_radius_quantile(dist_mat, n_bins = 12)
df_chi$dist_class <- cut(
  df_chi$distance,
  breaks = radius,
  include.lowest = TRUE
)


get_midpoints <- function(breaks) {
  (breaks[-1] + breaks[-length(breaks)]) / 2
}

midpoints <- get_midpoints(radius)

chi_summary <- df_chi %>%
  group_by(dist_class, tau) %>%
  summarise(
    mean_chi = mean(chi, na.rm = TRUE),
    n_pairs = sum(!is.na(chi)),
    .groups = "drop"
  ) %>%
  mutate(
    h = midpoints[as.numeric(dist_class)]
  )


# Visualisation
ggplot(chi_summary, aes(x = h, y = mean_chi, color = as.factor(tau))) +
  geom_line() +
  geom_point() +
  labs(title = "Chi Spatio-Temporel Empirique",
       x = "Distance (m)", y = expression(chi(h, tau)),
       color = "Temporal lag (5 min)") +
  theme_minimal()


# from chi compute variogram
vario_spatemp <- function(chi) {
  chi <- pmin(pmax(chi, 1e-6), 1 - 1e-6)
  invphi <- qnorm(1 - 0.5 * chi)
  2 * invphi^2
}
df_vario <- chi_summary %>%
  mutate(vario = vario_spatemp(mean_chi))

ggplot(df_vario, aes(x = h, y = vario, color = as.factor(tau))) +
  geom_line() +
  geom_point() +
  labs(
       x = "Distance (m)", y = expression(gamma(h, tau)),
       color = "Temporal lag (5 min)") +
  theme_minimal()




# loess smoothing
ggplot(df_vario, aes(x = h, y = vario, color = as.factor(tau))) +
  geom_point() +
  geom_smooth(method = "loess", se = F, span =2) +
  labs(title = "",
       x = "Distance (m)", y = expression(gamma(h, tau)),
       color = "Temporal lag (5 min)") +
  btf_theme

# save 
ggsave(paste0(im_folder, "variogram/empirical/omsev/variogram_separability/chi_spatemp_q", q * 100, ".png"), width = 8, height = 7)


library(dplyr)

# référence tau = 0
chi_ref_df <- chi_summary %>%
  filter(tau == 0) %>%
  select(h, chi_ref = mean_chi)

# calcul du ratio
Ratio <- chi_summary %>%
  filter(tau > 0) %>%
  left_join(chi_ref_df, by = "h") %>%
  mutate(ratio = mean_chi / chi_ref)

ggplot(Ratio, aes(x = h, y = ratio, color = as.factor(tau))) +
    geom_line() +
    geom_point() +
    labs(title = "Ratio of Chi to Tau=0",
         x = "Distance (m)", y = expression(chi(h, tau) / chi(h, 0)),
         color = "Temporal lag (5 min)") +
    theme_minimal()

# do it with gamma(h, tau) / gamma(h, 0)
vario_ref_df <- df_vario %>%
  filter(tau == 0) %>%
  select(h, vario_ref = vario)
Ratio_vario <- df_vario %>%
  filter(tau > 0) %>%
  left_join(vario_ref_df, by = "h") %>%
  mutate(ratio_vario = vario / vario_ref)

ggplot(Ratio_vario, aes(x = h, y = ratio_vario, color = as.factor(tau))) +
    geom_line() +
    geom_point() +
    labs(title = "Ratio of Variogram to Tau=0",
         x = "Distance (m)", y = expression(gamma(h, tau) / gamma(h, 0)),
         color = "Temporal lag (5 min)") +
    theme_minimal()


# same gamma(h,tau) - gamma(h, 0)
Ratio_vario_diff <- df_vario %>%
  filter(tau > 0) %>%
  left_join(vario_ref_df, by = "h") %>%
  mutate(diff_vario = vario - vario_ref)

ggplot(Ratio_vario_diff, aes(x = h, y = diff_vario, color = as.factor(tau))) +
    geom_smooth(method = "loess", se = F, span = 2) +
    labs(title = "Difference of Variogram to Tau=0",
         x = "Distance (m)", y = expression(gamma(h, tau) - gamma(h, 0)),
         color = "Temporal lag (5 min)") +
    theme_minimal()


# with a spatial ref h0
h0 <- min(df_vario$h)
Ratio_vario_diff_h0 <- df_vario %>%
  left_join(df_vario %>% filter(h == h0) %>% select(tau, vario_h0 = vario), by = "tau") %>%
  mutate(diff_vario_h0 = vario - vario_h0)
ggplot(Ratio_vario_diff_h0, aes(x = h, y = diff_vario_h0, color = as.factor(tau))) +
    geom_smooth(method = "loess", se = F, span = 2) +
    labs(title = "",
         x = "Distance (m)", y = expression(gamma(h, tau) - gamma(0, tau)),
         color = "Temporal lag (5 min)") +
    btf_theme

# save
ggsave(paste0(im_folder, "variogram/empirical/omsev/2025/variogram_separability/vario_diff_h0_q", q * 100, ".pdf"), width = 10, height = 5)

