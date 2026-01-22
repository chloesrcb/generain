# remove objects from the environment
rm(list = ls())

# clear console
cat("\014")

# Load libraries and set theme
source("./script/load_libraries.R")

functions_folder <- "./R"
files <- list.files(functions_folder, full.names = TRUE)
# # load all functions in files
invisible(lapply(files, function(f) source(f, echo = FALSE)))
library(latex2exp)

################################################################################
# DATA 2022 --------------------------------------------------------------------
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

# get rain measurements
# load data
# get rain data from omsev
filename_omsev <- paste0(data_folder,
                         "omsev/omsev_5min/rain_mtp_5min_2019_2022.RData")

load(filename_omsev)
rain <- rain.all5[c(1, 6:ncol(rain.all5))]
typeof(rain)
class(rain)
colnames(rain)
rownames(rain) <- rain$dates
rain_new <- rain[-1] # remove dates column

# TEMPORAL CHI -----------------------------------------------------------------

q <- 0.9 # quantile
tmax <- 10
nsites <- ncol(rain_new)
chimat_dtlag <- temporal_chi(rain_new, tmax, quantile = q, mean = FALSE, zeros=F)


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

# Save the plot
filename <- paste(im_folder, "WLSE/omsev/2022/bp_temporal_chi_", q,
                  ".pdf", sep = "")
ggsave(filename, width = 20, height = 15, units = "cm")


# for multiple quantiles
qs <- seq(0.9, 0.99, by = 0.01)
tmax <- 10  # par exemple

# Initialiser une liste pour stocker les résultats
chi_list <- lapply(qs, function(q) {
  chi_vals <- temporal_chi(rain_new, tmax, quantile = q, mean = TRUE, zeros = F)
  chi_vals <- chi_vals[-1]  # enlever lag 0
  data.frame(
    lag = (1:tmax) * 5,
    chi = chi_vals,
    q = q
  )
})

# Combiner tous les data.frames en un seul
chi_df <- bind_rows(chi_list)

# Tracer
chitemp_multi <- ggplot(chi_df, aes(x = lag, y = chi, color = factor(q))) +
  geom_line(alpha = 0.8, linewidth = 1.2) +
  xlab(TeX(r"($\tau$ (minutes))")) +
  ylab(TeX(r"($\widehat{\chi}(0,\tau)$)")) +
  ylim(0, 1) +
  scale_color_brewer(palette = "BuPu", name = "Quantile") +
  theme_minimal() +
  theme(panel.grid.major = element_line(color = "white"),
        panel.grid.minor = element_line(color = "white")) +
  btf_theme

chitemp_multi

# Save the plot
filename <- paste(im_folder, "WLSE/omsev/2022/temporal_chi_multi_q_zeros.pdf", sep = "")
ggsave(filename, width = 20, height = 15, units = "cm")


# WLSE transformation ----------------------------------------------------------
quantiles <- c(0.9, 0.95, 0.99, 0.992, 0.995)  # quantiles to consider
tmax <- 10  # adapte si besoin

results <- data.frame()
points_df <- data.frame()
lines_df <- data.frame()


for (q in quantiles) {
  chimat_dt_mean <- temporal_chi(rain_new, tmax, quantile = q, mean = TRUE, zeros = FALSE)
  df_chi <- data.frame(lag = c(0:tmax) * 5, chi = chimat_dt_mean)

  wlse_temp <- get_estimate_variotemp(df_chi, weights = "exp", summary = F, lag_unit = 5)
  
  c2 <- as.numeric(wlse_temp[[1]])
  beta2 <- as.numeric(wlse_temp[[2]])
  alpha2 <- as.numeric(wlse_temp[[3]])
  signif_beta <- wlse_temp[[4]]
  signif_alpha <- wlse_temp[[5]]

  coef_table <- data.frame(
    quantile = q,
    parameter = c("beta", "alpha"),
    value = c(beta2, alpha2),
    significance = c(signif_beta, signif_alpha)
  )
  results <- rbind(results, coef_table)

  df_chi_nz <- df_chi[-1, ]
  dftemp <- data.frame(
    lag = log(df_chi_nz$lag),
    chi = zeta(df_chi_nz$chi),
    q = factor(q)
  )
  points_df <- rbind(points_df, dftemp)

  lag_seq <- seq(0, tmax * 5, by = 1)
  lag_seq <- log(lag_seq[lag_seq > 0])  # log transformation
  fit_line <- data.frame(
    lag = lag_seq,
    chi = alpha2 * lag_seq + c2,
    q = factor(q)
  )
  lines_df <- rbind(lines_df, fit_line)
}

results

# get table with quantiles, beta and alpha as columns
results <- results %>%
  pivot_wider(names_from = parameter, values_from = value) %>%
  mutate(
    quantile = as.character(quantile),
    beta = round(beta, 3),
    alpha = round(alpha, 3)
  ) %>%
  select(quantile, beta, alpha, significance)


library(ggplot2)
library(latex2exp)  # pour TeX() si utilisé

ggplot() +
  geom_point(data = points_df, aes(x = lag, y = chi, color = q), size = 3, alpha = 0.7) +
  geom_line(data = lines_df, aes(x = lag, y = chi, color = q), size = 1.2, alpha = 0.8) +
  btf_theme +  # ton thème perso
  xlab(TeX(r"($\log(\tau)$)")) +
  ylab(TeX(r"($\zeta(\widehat{\chi}(0,\tau))$)")) +
  scale_color_brewer(palette = "Purples") +  # ou autre palette adaptée
  theme(
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10),
    panel.grid.major = element_line(color = "white"),
    panel.grid.minor = element_line(color = "white")
  )

# Save the plot
filename <- paste(im_folder, "WLSE/omsev/2022/temporal_chi_wlse_multiple_q.pdf", sep = "")
ggsave(filename, width = 20, height = 15, units = "cm")

quantiles <- c(seq(0.9, 0.99, by = 0.01), 0.992, 0.995)
tmax <- 10

for (q in quantiles) {
    chimat_dt_mean <- temporal_chi(rain_new, tmax, quantile = q, mean = TRUE, zeros=F)
    # remove lag 0
    chimat_dt_mean <- chimat_dt_mean[-1]


    # plot mean chi
    df_chi <- data.frame(lag = c(1:tmax) * 5, chi = chimat_dt_mean)
    wlse_temp <- get_estimate_variotemp(df_chi, 
                                        weights = "exp", summary = F, lag_unit = 5)

    c2 <- as.numeric(wlse_temp[[1]])
    beta2 <- as.numeric(wlse_temp[[2]])
    alpha2 <- as.numeric(wlse_temp[[3]])

    df <- df_chi
    dftemp <- data.frame(lag = log(df$lag), chi = zeta(df$chi))

    p <- ggplot(dftemp, aes(x = lag, y = chi)) +
            geom_point(color = btfgreen, size = 4) +
            btf_theme +
            xlab(TeX(r"($\log(\tau)$)")) +
            ylab(TeX(r"($\zeta(\widehat{\chi}(0,\tau))$)")) +
            geom_line(aes(x = lag, y = alpha2 * lag + c2),
                        alpha = 0.6, color = "darkred", linewidth = 1.5)
    # save the plot
    filename <- paste(im_folder, "WLSE/omsev/2022/temporal_chi_",
                        as.character(q), ".pdf", sep = "")
    ggsave(filename, plot = p, width = 20, height = 15, units = "cm")
}



################################################################################
# DATA 2025 -------------------------------------------------------------
################################################################################

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



rain <- rain.all5[c(1, 6:ncol(rain.all5) - 1)]
typeof(rain)
class(rain)
colnames(rain)

rownames(rain) <- rain$dates
rain_new <- rain[-1]


# TEMPORAL CHI -----------------------------------------------------------------

q <- 0.98 # quantile
tmax <- 10
nsites <- ncol(rain_new)
chimat_dtlag <- temporal_chi(rain_new, tmax, quantile = q, mean = FALSE, zeros=F)


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

# Save the plot
filename <- paste(im_folder, "WLSE/omsev/2025/bp_temporal_chi_", q,
                  ".pdf", sep = "")
ggsave(filename, width = 20, height = 15, units = "cm")

# for multiple quantiles
qs <- seq(0.9, 0.995, by = 0.01)
tmax <- 10  # par exemple

# Initialiser une liste pour stocker les résultats
chi_list <- lapply(qs, function(q) {
  chi_vals <- temporal_chi(rain_new, tmax, quantile = q, mean = TRUE, zeros = FALSE)
  chi_vals <- chi_vals[-1]  # enlever lag 0
  data.frame(
    lag = (1:tmax) * 5,
    chi = chi_vals,
    q = q
  )
})

# Combiner tous les data.frames en un seul
chi_df <- bind_rows(chi_list)

# Tracer
chitemp_multi <- ggplot(chi_df, aes(x = lag, y = chi, color = factor(q))) +
  geom_line(alpha = 0.8, linewidth = 1.2) +
  xlab(TeX(r"($\tau$ (minutes))")) +
  ylab(TeX(r"($\widehat{\chi}(0,\tau)$)")) +
  ylim(0, 1) +
  scale_color_brewer(palette = "BuPu", name = "Quantile") +
  theme_minimal() +
  theme(panel.grid.major = element_line(color = "white"),
        panel.grid.minor = element_line(color = "white")) +
  btf_theme

chitemp_multi

# Save the plot
filename <- paste(im_folder, "WLSE/omsev/2025/temporal_chi_multi_q.pdf", sep = "")
ggsave(filename, width = 20, height = 15, units = "cm")



# WLSE transformation ----------------------------------------------------------
quantiles <- c(0.9, 0.95, 0.99, 0.992, 0.995)  # quantiles to consider
tmax <- 10  # adapte si besoin

results <- data.frame()
points_df <- data.frame()
lines_df <- data.frame()


for (q in quantiles) {
  chimat_dt_mean <- temporal_chi(rain_new, tmax, quantile = q, mean = TRUE, zeros = FALSE)
  df_chi <- data.frame(lag = c(0:tmax) * 5, chi = chimat_dt_mean)

  wlse_temp <- get_estimate_variotemp(df_chi, weights = "exp", summary = F, lag_unit = 5)
  
  c2 <- as.numeric(wlse_temp[[1]])
  beta2 <- as.numeric(wlse_temp[[2]])
  alpha2 <- as.numeric(wlse_temp[[3]])
  signif_beta <- wlse_temp[[4]]
  signif_alpha <- wlse_temp[[5]]

  coef_table <- data.frame(
    quantile = q,
    parameter = c("beta", "alpha"),
    value = c(beta2, alpha2),
    significance = c(signif_beta, signif_alpha)
  )
  results <- rbind(results, coef_table)

  df_chi_nz <- df_chi[-1, ]
  dftemp <- data.frame(
    lag = log(df_chi_nz$lag),
    chi = zeta(df_chi_nz$chi),
    q = factor(q)
  )
  points_df <- rbind(points_df, dftemp)

  lag_seq <- seq(0, tmax * 5, by = 1)
  lag_seq <- log(lag_seq[lag_seq > 0])  # log transformation
  fit_line <- data.frame(
    lag = lag_seq,
    chi = alpha2 * lag_seq + c2,
    q = factor(q)
  )
  lines_df <- rbind(lines_df, fit_line)
}

results

# get table with quantiles, beta and alpha as columns
results <- results %>%
  pivot_wider(names_from = parameter, values_from = value) %>%
  mutate(
    quantile = as.character(quantile),
    beta = round(beta, 3),
    alpha = round(alpha, 3)
  ) %>%
  select(quantile, beta, alpha, significance)


library(ggplot2)
library(latex2exp)  # pour TeX() si utilisé

ggplot() +
  geom_point(data = points_df, aes(x = lag, y = chi, color = q), size = 3, alpha = 0.7) +
  geom_line(data = lines_df, aes(x = lag, y = chi, color = q), size = 1.2, alpha = 0.8) +
  btf_theme +  # ton thème perso
  xlab(TeX(r"($\log(\tau)$)")) +
  ylab(TeX(r"($\zeta(\widehat{\chi}(0,\tau))$)")) +
  scale_color_brewer(palette = "Purples") +  # ou autre palette adaptée
  theme(
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10),
    panel.grid.major = element_line(color = "white"),
    panel.grid.minor = element_line(color = "white")
  )

# Save the plot
filename <- paste(im_folder, "WLSE/omsev/2025/temporal_chi_wlse_multiple_q.pdf", sep = "")
ggsave(filename, width = 20, height = 15, units = "cm")

quantiles <- c(seq(0.9, 0.99, by = 0.01), 0.992, 0.995)
tmax <- 10

for (q in quantiles) {
    chimat_dt_mean <- temporal_chi(rain_new, tmax, quantile = q, mean = TRUE, zeros=F)
    # remove lag 0
    chimat_dt_mean <- chimat_dt_mean[-1]


    # plot mean chi
    df_chi <- data.frame(lag = c(1:tmax) , chi = chimat_dt_mean)
    wlse_temp <- get_estimate_variotemp(df_chi, 
                                    weights = "exp", summary = F, lag_unit = 1)

    c2 <- as.numeric(wlse_temp[[1]])
    beta2 <- as.numeric(wlse_temp[[2]])
    alpha2 <- as.numeric(wlse_temp[[3]])

    df <- df_chi
    dftemp <- data.frame(lag = log(df$lag), chi = zeta(df$chi))

    p <- ggplot(dftemp, aes(x = lag, y = chi)) +
            geom_point(color = btfgreen, size = 4) +
            btf_theme +
            xlab(TeX(r"($\log(\tau)$)")) +
            ylab(TeX(r"($\zeta(\widehat{\chi}(0,\tau))$)")) +
            geom_line(aes(x = lag, y = alpha2 * lag + c2),
                        alpha = 0.6, color = "darkred", linewidth = 1.5)
    # save the plot
    filename <- paste(im_folder, "WLSE/omsev/2025/temporal_chi_index_",
                        as.character(q), ".pdf", sep = "")
    ggsave(filename, plot = p, width = 20, height = 15, units = "cm")
}


########################################################################################
# COMEPHORE DATA -------------------------------------------------------------
########################################################################################

# get rain data from comephore
filename_com <- paste0(data_folder, "comephore/zoom_5km.csv")
comephore_raw <- read.csv(filename_com, sep = ",")
filename_loc <- paste0(data_folder, "comephore/loc_px_zoom_5km.csv")
loc_px <- read.csv(filename_loc, sep = ",")

df_comephore <- comephore_raw
head(df_comephore)

# Take only data after 2007
colnames(df_comephore)[1] <- "date"
df_comephore <- df_comephore[df_comephore$date >= "2008-01-01", ]

# put date in index
rownames(df_comephore) <- df_comephore$date
comephore <- df_comephore[-1] # remove dates column
rain_new <- comephore
# Get distances matrix
dist_mat <- get_dist_mat(loc_px)
nsites <- nrow(loc_px)


# TEMPORAL CHI -----------------------------------------------------------------

q <- 0.98 # quantile
tmax <- 10
nsites <- ncol(rain_new)
chimat_dtlag <- temporal_chi(rain_new, tmax, quantile = q, mean = FALSE, zeros=F)


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

# Save the plot
filename <- paste(im_folder, "WLSE/comephore/bp_temporal_chi_", q,
                  ".pdf", sep = "")
ggsave(filename, width = 20, height = 15, units = "cm")

# for multiple quantiles
qs <- seq(0.9, 0.995, by = 0.01)
tmax <- 10  # par exemple

# Initialiser une liste pour stocker les résultats
chi_list <- lapply(qs, function(q) {
  chi_vals <- temporal_chi(rain_new, tmax, quantile = q, mean = TRUE, zeros = FALSE)
  chi_vals <- chi_vals[-1]  # enlever lag 0
  data.frame(
    lag = (1:tmax),
    chi = chi_vals,
    q = q
  )
})

# Combiner tous les data.frames en un seul
chi_df <- bind_rows(chi_list)

# Tracer
chitemp_multi <- ggplot(chi_df, aes(x = lag, y = chi, color = factor(q))) +
  geom_line(alpha = 0.8, linewidth = 1.2) +
  xlab(TeX(r"($\tau$ (minutes))")) +
  ylab(TeX(r"($\widehat{\chi}(0,\tau)$)")) +
  ylim(0, 1) +
  scale_color_brewer(palette = "BuPu", name = "Quantile") +
  theme_minimal() +
  theme(panel.grid.major = element_line(color = "white"),
        panel.grid.minor = element_line(color = "white")) +
  btf_theme

chitemp_multi

# Save the plot
filename <- paste(im_folder, "WLSE/comephore/temporal_chi_multi_q.pdf", sep = "")
ggsave(filename, width = 20, height = 15, units = "cm")



# WLSE transformation ----------------------------------------------------------
quantiles <- c(0.9, 0.95, 0.99, 0.992, 0.995, 0.998)  # quantiles to consider
# quantiles <- c(seq(0.996, 0.9995, by = 0.0005))
tmax <- 10  # adapte si besoin

results <- data.frame()
points_df <- data.frame()
lines_df <- data.frame()


for (q in quantiles) {
  chimat_dt_mean <- temporal_chi(rain_new, tmax, quantile = q, mean = TRUE, zeros = TRUE)
  df_chi <- data.frame(lag = c(0:tmax), chi = chimat_dt_mean)

  wlse_temp <- get_estimate_variotemp(df_chi, weights = "exp", summary = F, lag_unit = 1)
  
  c2 <- as.numeric(wlse_temp[[1]])
  beta2 <- as.numeric(wlse_temp[[2]])
  alpha2 <- as.numeric(wlse_temp[[3]])
  signif_beta <- wlse_temp[[4]]
  signif_alpha <- wlse_temp[[5]]

  coef_table <- data.frame(
    quantile = q,
    parameter = c("beta", "alpha"),
    value = c(beta2, alpha2),
    significance = c(signif_beta, signif_alpha)
  )
  results <- rbind(results, coef_table)

  df_chi_nz <- df_chi[-1, ]
  dftemp <- data.frame(
    lag = log(df_chi_nz$lag),
    chi = zeta(df_chi_nz$chi),
    q = factor(q)
  )
  points_df <- rbind(points_df, dftemp)

  lag_seq <- seq(0, tmax, by = 1)
  lag_seq <- log(lag_seq[lag_seq > 0])  # log transformation
  fit_line <- data.frame(
    lag = lag_seq,
    chi = alpha2 * lag_seq + c2,
    q = factor(q)
  )
  lines_df <- rbind(lines_df, fit_line)
}

results

# get table with quantiles, beta and alpha as columns
results_2 <- results %>%
  pivot_wider(names_from = parameter, values_from = value) %>%
  mutate(
    quantile = as.character(quantile),
    beta = round(beta, 3),
    alpha = round(alpha, 3)
  ) %>%
  select(quantile, beta, alpha, significance)


library(ggplot2)
library(latex2exp)  # pour TeX() si utilisé

ggplot() +
  geom_point(data = points_df, aes(x = lag, y = chi, color = q), size = 3, alpha = 0.7) +
  geom_line(data = lines_df, aes(x = lag, y = chi, color = q), size = 1.2, alpha = 0.8) +
  btf_theme +  # ton thème perso
  xlab(TeX(r"($\log(\tau)$)")) +
  ylab(TeX(r"($\zeta(\widehat{\chi}(0,\tau))$)")) +
  scale_color_brewer(palette = "Purples") +  # ou autre palette adaptée
  theme(
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10),
    panel.grid.major = element_line(color = "white"),
    panel.grid.minor = element_line(color = "white")
  )

# Save the plot
filename <- paste(im_folder, "WLSE/comephore/temporal_chi_wlse_multiple_q.pdf", sep = "")
ggsave(filename, width = 20, height = 15, units = "cm")

quantiles <- c(seq(0.9, 0.99, by = 0.01), 0.992, 0.995)
tmax <- 10

for (q in quantiles) {
    chimat_dt_mean <- temporal_chi(rain_new, tmax, quantile = q, mean = TRUE, zeros=F)
    # remove lag 0
    chimat_dt_mean <- chimat_dt_mean[-1]


    # plot mean chi
    df_chi <- data.frame(lag = c(1:tmax) , chi = chimat_dt_mean)
    wlse_temp <- get_estimate_variotemp(df_chi, 
                                    weights = "exp", summary = F, lag_unit = 1)

    c2 <- as.numeric(wlse_temp[[1]])
    beta2 <- as.numeric(wlse_temp[[2]])
    alpha2 <- as.numeric(wlse_temp[[3]])

    df <- df_chi
    dftemp <- data.frame(lag = log(df$lag), chi = zeta(df$chi))

    p <- ggplot(dftemp, aes(x = lag, y = chi)) +
            geom_point(color = btfgreen, size = 4) +
            btf_theme +
            xlab(TeX(r"($\log(\tau)$)")) +
            ylab(TeX(r"($\zeta(\widehat{\chi}(0,\tau))$)")) +
            geom_line(aes(x = lag, y = alpha2 * lag + c2),
                        alpha = 0.6, color = "darkred", linewidth = 1.5)
    # save the plot
    filename <- paste(im_folder, "WLSE/comephore/temporal_chi_",
                        as.character(q), ".pdf", sep = "")
    ggsave(filename, plot = p, width = 20, height = 15, units = "cm")
}
