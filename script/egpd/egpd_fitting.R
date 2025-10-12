# libraries
# remove objects from the environment
rm(list = ls())

# clear console
cat("\014")

# load libraries
library(generain)
source("./script/load_libraries.R")

################################################################################
# DATA ---------------------------------------------------------------------
################################################################################
folder_rain <- paste0(data_folder, "omsev/omsev_5min/")
load(file = paste0(folder_rain, "rain_mtp_5min_2019_2024.RData"))
rain <- rain.all5[c(1, 6:(ncol(rain.all5)-1))] 
rownames(rain) <- rain$dates
rain_new <- rain[-1] # remove dates column

# get rain data from comephore
filename_com <- paste0(data_folder, "comephore/zoom_5km.csv")
comephore_raw <- read.csv(filename_com, sep = ",")
df_comephore <- comephore_raw
head(df_comephore)

# remova "hydro", "brives" and  "cines" column
index_col <- which(colnames(rain_new) != "hydro" & colnames(rain_new) != "brives" & colnames(rain_new) != "cines")
rain_sub <- rain_new[ ,index_col]
colnames(rain_sub)

# Take only data 3 years of data from COMEPHORE
colnames(df_comephore)[1] <- "date"
df_comephore <- df_comephore[df_comephore$date >= "2019-01-01", ]
length(df_comephore[,1]) / 24 / 365
rain_com <- df_comephore[-1] # remove dates column
ncol(rain_com)

################################################################################
# EGPD ---------------------------------------------------------------------
################################################################################

library(ggplot2)
library(dplyr)

## ---------- Fonction utilitaire pour compter fréquences et tracer ----------
make_freq_barplot <- function(data, low, high, dataset_name, im_folder, rain_unit = "mm/5min") {
  # filtrer et arrondir
  df <- data[data > low & data < high]
  df <- data.frame(Rain = round(df, 4))
  
  # compter occurrences exactes
  freq_df <- df %>%
    group_by(Rain) %>%
    summarise(count = n(), .groups = "drop")
  
  # labels = valeurs uniques
  labels_name <- freq_df$Rain
  
  # barplot exact
  hist_plot <- ggplot(freq_df, aes(x = Rain, y = count)) +
    geom_col(fill = "lightblue", color = "black", width = 0.01) +
    geom_text(aes(label = labels_name), vjust = -0.5, size = 2.5) +
    labs(
      title = paste0(dataset_name, " rain > ", low, " and < ", high, " ", rain_unit),
      x = paste("Rain (", rain_unit, ")", sep = ""), y = "Count"
    ) +
    btf_theme
  
  # sauvegarde du pdf
  filename <- paste0(im_folder, "EGPD/barplot_rain_", dataset_name, "_", low, "_", high, ".pdf")
  ggsave(filename = filename, plot = hist_plot, device = "pdf", bg = "transparent",
         width = 8, height = 6)
  
  return(list(freq_df = freq_df, labels_name = labels_name, plot = hist_plot))
}
##

compute_rmse_nrmse <- function(y, left_censore) {
  inits <- init_values(y, 0)
  sigma_0 <- inits[1]
  xi_0 <- inits[2]
  kappa_0 <- 1

  fit <- tryCatch({
    fit.extgp(y, model = 1, method = "mle",
              init = c(kappa_0, sigma_0, xi_0),
              censoring = c(left_censore, Inf),
              plots = FALSE, confint = FALSE, ncpus = 1)
  }, error = function(e) return(NULL))

  if (is.null(fit)) return(data.frame(RMSE = NA, NRMSE = NA))

  param <- fit$fit$mle
  probs <- ppoints(length(y))
  q_theo <- qextgp(probs, type = 1, kappa = param[1], sigma = param[2], xi = param[3])
  q_emp <- sort(y)

  rmse <- sqrt(mean((q_emp - q_theo)^2))
  nrmse <- (rmse / mean(y))

  return(data.frame(RMSE = rmse, NRMSE = nrmse))
}


# Fonction pour générer le tableau pour un vecteur de censures
rmse_nrmse_table <- function(y, censure_values) {
  results <- lapply(censure_values, function(c) {
    res <- compute_rmse_nrmse(y, c)
    res$Censoring <- c
    return(res)
  })
  do.call(rbind, results)
}

## ---------- COMEPHORE ----------
# Example for one site
folder_censure <- paste0(im_folder, "/EGPD/COMEPHORE/censure_choice/")
if (!dir.exists(folder_censure)) {
  dir.create(folder_censure, recursive = TRUE)
}

folder_fit <- paste0(im_folder, "/EGPD/COMEPHORE/sites/")
if (!dir.exists(folder_fit)) {
  dir.create(folder_fit, recursive = TRUE)
}

# Example loop over all pixels
rain_com <- na.omit(rain_com)
rain_data <- rain_com  # Data frame with all pixels
censure_values <- seq(0.1, 1, by = 0.1)
all_censures <- list()
all_results <- list()
for (i in seq_along(rain_data)) {
  site_name <- colnames(rain_data)[i]
  cat("Processing site:", site_name, "\n")
  y_site <- na.omit(rain_data[[i]])
  y_site <- y_site[y_site > 0]
  if (length(y_site) < 20) next

  # Compute RMSE and NRMSE for all censure values
  rmse_vals <- sapply(censure_values, function(c) compute_rmse(y_site, c))
  nrmse_vals <- sapply(censure_values, function(c) compute_nrmse(y_site, c))
  table_cens <- data.frame(Censoring = censure_values, RMSE = rmse_vals, NRMSE = nrmse_vals)

  # Save table
  write.csv(table_cens, file = paste0(folder_censure, site_name, ".csv"),
                                      row.names = FALSE)

  # Select best censure (minimum RMSE or NRMSE)
  best_cens_rmse <- table_cens$Censoring[which.min(table_cens$RMSE)]
  best_cens_nrmse <- table_cens$Censoring[which.min(table_cens$NRMSE)]

  # Fit EGPD and save QQ-plot for RMSE
  rmse_folder <- paste0(folder_fit, "RMSE/")
  # fit_rmse <- process_site(y_site, site_name, rmse_folder, best_cens_rmse)

  # Fit EGPD and save QQ-plot for NRMSE
  nrmse_folder <- paste0(folder_fit, "NRMSE/")
  # fit_nrmse <- process_site(y_site, site_name, nrmse_folder, best_cens_nrmse)

  # all_results[[i]] <- list(RMSE_fit = fit_rmse, NRMSE_fit = fit_nrmse)
  all_censures[[i]] <- c(BestCens_RMSE = best_cens_rmse, BestCens_NRMSE = best_cens_nrmse)
}

# histogram of best censure values for RMSE
best_censures_rmse_com <- sapply(all_censures, function(res) res[["BestCens_RMSE"]])
hist_best_cens_rmse <- ggplot(data.frame(Best_Cens = best_censures_rmse_com), aes(x = Best_Cens)) +
  geom_histogram(binwidth = 0.1, fill = btfgreen, color = "black", alpha = 0.5) +
  labs(
    title = "",
    x = "Best left-censoring", y = "Count"
  ) +
  btf_theme
hist_best_cens_rmse

filename <- paste0(im_folder, "EGPD/COMEPHORE/hist_best_cens_rmse.pdf")
ggsave(filename = filename, plot = hist_best_cens_rmse, device = "pdf", bg = "transparent",
       width = 8, height = 6)

best_censures_nrmse <- sapply(all_censures, function(res) res[["BestCens_NRMSE"]])
hist_best_cens_nrmse <- ggplot(data.frame(Best_Cens = best_censures_nrmse), aes(x = Best_Cens)) +
  geom_histogram(binwidth = 0.1, fill = btfgreen, color = "black", alpha = 0.5) +
  labs(
    title = "",
    x = "Best left-censoring", y = "Count"
  ) +
  btf_theme
hist_best_cens_nrmse

filename <- paste0(im_folder, "EGPD/COMEPHORE/hist_best_cens_nrmse.pdf")
ggsave(filename = filename, plot = hist_best_cens_nrmse, device = "pdf", bg = "transparent",
       width = 8, height = 6)

## ---------- OMSEV ----------
# Example for one site
folder_censure <- paste0(im_folder, "/EGPD/OMSEV/censure_choice/")
if (!dir.exists(folder_censure)) {
  dir.create(folder_censure, recursive = TRUE)
}

folder_fit <- paste0(im_folder, "/EGPD/OMSEV/sites/")
if (!dir.exists(folder_fit)) {
  dir.create(folder_fit, recursive = TRUE)
}


# Do fit for poly site to compare possible censure values
rain_poly <- rain_sub$poly
folder_poly <- paste0(im_folder, "/EGPD/OMSEV/poly_compa/")
y_poly <- na.omit(rain_poly)
y_poly <- y_poly[y_poly > 0]
# possible step values in poly
y_max08 <- y_poly[y_poly < 0.8]
possible_values <- sort(unique(round(y_max08, 2)))
censures_values <- seq(min(possible_values), max(possible_values), by = 0.01)

table_cens_poly <- rmse_nrmse_table(y_poly, censures_values)
write.csv(table_cens_poly, file = paste0(folder_censure, "poly.csv"),
          row.names = FALSE)
censure_values_poly <- censures_values
for (c in censure_values_poly) {
  cat("Censure:", c, "\n")
  fit_poly <- process_site(y_poly, "poly", folder_poly, best_cens = c, R = 1000)
}

rain_omsev <- na.omit(as.data.frame(rain_sub))
rain_data <- rain_omsev  # Data frame with all pixels
censure_values <- seq(0.22, 0.8, by = 0.01)
# possible_values <- sort(unique(round(rain_omsev[rain_omsev < 0.8 & rain_omsev > 0], 4)))
# censure_values <- possible_values

all_results <- list()
for (i in seq_along(rain_data)) {
  site_name <- colnames(rain_data)[i]
  cat("Processing site:", site_name, "\n")
  y_site <- na.omit(rain_data[[i]])
  y_site <- y_site[y_site > 0]
  if (length(y_site) < 20) next

  # Compute RMSE and NRMSE for all censure values
  rmse_vals <- sapply(censure_values, function(c) compute_rmse(y_site, c))
  nrmse_vals <- sapply(censure_values, function(c) compute_nrmse(y_site, c))
  table_cens <- data.frame(Censoring = censure_values, RMSE = rmse_vals,
                            NRMSE = nrmse_vals)

  # Save table
  write.csv(table_cens, file = paste0(folder_censure, site_name, ".csv"),
            row.names = FALSE)

  # Select best censure (minimum RMSE or NRMSE)
  best_cens_rmse <- table_cens$Censoring[which.min(table_cens$RMSE)]
  best_cens_nrmse <- table_cens$Censoring[which.min(table_cens$NRMSE)]

  # Fit EGPD and save QQ-plot for RMSE
  rmse_folder <- paste0(folder_fit, "RMSE/")
  fit_rmse <- process_site(y_site, site_name, rmse_folder, best_cens_rmse)

  # Fit EGPD and save QQ-plot for NRMSE
  nrmse_folder <- paste0(folder_fit, "NRMSE/")
  fit_nrmse <- process_site(y_site, site_name, nrmse_folder, best_cens_nrmse)

  all_results[[i]] <- list(RMSE_fit = fit_rmse, NRMSE_fit = fit_nrmse)
}

# histogram of best censure values for RMSE
best_censures_rmse_omsev <- sapply(all_results, function(res) res$RMSE_fit$BestCens)
hist_best_cens_rmse <- ggplot(data.frame(Best_Cens = best_censures_rmse_omsev), aes(x = Best_Cens)) +
  geom_histogram(binwidth = 0.001, fill = btfgreen, color = "black", alpha = 0.5) +
  labs(
    title = "",
    x = "Best left-censoring (RMSE)", y = "Count"
  ) +
  btf_theme
hist_best_cens_rmse

filename <- paste0(im_folder, "EGPD/OMSEV/hist_best_cens_rmse.pdf")
ggsave(filename = filename, plot = hist_best_cens_rmse, device = "pdf", bg = "transparent",
       width = 8, height = 6)

best_censures_nrmse_omsev <- sapply(all_results, function(res) res$NRMSE_fit$BestCens)
hist_best_cens_nrmse <- ggplot(data.frame(Best_Cens = best_censures_nrmse_omsev), aes(x = Best_Cens)) +
  geom_histogram(binwidth = 0.001, fill = btfgreen, color = "black", alpha = 0.5) +
  labs(
    title = "",
    x = "Best left-censoring (NRMSE)", y = "Count"
  ) +
  btf_theme
hist_best_cens_nrmse

filename <- paste0(im_folder, "EGPD/OMSEV/hist_best_cens_nrmse.pdf")
ggsave(filename = filename, plot = hist_best_cens_nrmse, device = "pdf", bg = "transparent",
       width = 8, height = 6)


## ---------- Fusion of results COMEPHORE and OMSEV ----------

params_subrain <- get_egpd_estimates(
  rain_com,
  left_censoring = best_censures_rmse_com
)


params_subrain <- get_egpd_estimates(
  rain_omsev,
  left_censoring = best_censures_rmse_omsev
)


df_long_com <- get_df_long_params_egpd(params_com_mtp) %>%
  mutate(Dataset = "COMEPHORE")

df_long_ohsm <- get_df_long_params_egpd(params_subrain) %>%
  mutate(Dataset = "OMSEV")

df_combined <- bind_rows(df_long_com, df_long_ohsm)

params_boxplot <- ggplot(df_combined, aes(x = Variable, y = Value, fill = Dataset)) +
  geom_boxplot(alpha = 0.8, color = "#0000006c") +
  labs(x = "Parameters", y = "Distribution") +
  btf_theme +
  scale_x_discrete(labels = c(
    Xi = TeX(r"($\widehat{\xi}$)"),
    Sigma = TeX(r"($\widehat{\sigma}$)"),
    Kappa = TeX(r"($\widehat{\kappa}$)")
  )) +
  scale_fill_manual(values = c("COMEPHORE" = "#a55c5cbe", "OMSEV" = btfgreen)) +
  theme(
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 14)
  ) +
  ylim(0, 2)
  
# x11()

ggsave(
  filename = paste0(im_folder, "EGPD/params_estim_2DS_CensoreRMSE_not_continuous.pdf"),
  plot = params_boxplot, device = "pdf", bg = "transparent",
  width = 8, height = 6
)

params_boxplot_facet <- ggplot(df_combined, aes(x = Variable, y = Value, fill = Dataset)) +
      geom_boxplot(alpha = 0.8, color = "#0000006c") +
      labs(x = "Parameters", y = "Distribution") +
      btf_theme +
      scale_x_discrete(labels = c(
            Xi = TeX(r"($\widehat{\xi}$)"),
            Sigma = TeX(r"($\widehat{\sigma}$)"),
            Kappa = TeX(r"($\widehat{\kappa}$)")
      )) +
      scale_fill_manual(values = c("COMEPHORE" = "#a55c5cbe", "OMSEV" = btfgreen)) +
      facet_wrap(~Dataset, scales = "fixed") + 
      theme(
            legend.position = "none",
            strip.text = element_text(size = 14, face = "bold")
      ) +
      ylim(0, 1.7)

ggsave(
  filename = paste0(im_folder, "EGPD/params_estim_facet_byDataset_CensoreRMSE_not_continuous.pdf"),
  plot = params_boxplot_facet, device = "pdf", bg = "transparent",
  width = 10, height = 6
)


params_boxplot_com <- ggplot(df_long_com, aes(x = Variable, y = Value, fill = Dataset)) +
  geom_boxplot(alpha = 0.8, color = "#0000006c") +
  labs(x = "Parameters", y = "Distribution") +
  btf_theme +
  scale_x_discrete(labels = c(
    Xi = TeX(r"($\widehat{\xi}$)"),
    Sigma = TeX(r"($\widehat{\sigma}$)"),
    Kappa = TeX(r"($\widehat{\kappa}$)")
  )) +
  scale_fill_manual(values = c("COMEPHORE" = "#a55c5cbe")) +
  theme(
    legend.position = "none"
  )

ggsave(
  filename = paste0(im_folder, "EGPD/params_estim_COMEPHORE_CensoreRMSE_not_continuous.pdf"),
  plot = params_boxplot_com, device = "pdf", bg = "transparent",
  width = 8, height = 6
)

params_boxplot_ohsm <- ggplot(df_long_ohsm, aes(x = Variable, y = Value, fill = Dataset)) +
  geom_boxplot(alpha = 0.8, color = "#0000006c") +
  labs(x = "Parameters", y = "Distribution") +
  btf_theme +
  scale_x_discrete(labels = c(
    Xi = TeX(r"($\widehat{\xi}$)"),
    Sigma = TeX(r"($\widehat{\sigma}$)"),
    Kappa = TeX(r"($\widehat{\kappa}$)")
  )) +
  scale_fill_manual(values = c("OMSEV" = btfgreen)) +
      theme(
            legend.position = "none"
      )

ggsave(
  filename = paste0(im_folder, "EGPD/params_estim_OMSEV_CensoreRMSE_not_continuous.pdf"),
  plot = params_boxplot_ohsm, device = "pdf", bg = "transparent",
  width = 8, height = 6
)



results <- list()
sites_name <- colnames(rain_sub)

save_path <- paste0(im_folder, "EGPD/OMSEV/sites/")
if (!dir.exists(save_path)) {
  dir.create(save_path, recursive = TRUE)
}

for (site_name in sites_name) {
  cat("Processing:", site_name, "\n")
  y_raw <- as.data.frame(na.omit(rain_sub[[site_name]]))
  index_site <- which(colnames(rain_sub) == site_name)
  best_cens <- round(best_censures_omsev[index_site], 3)
  site_result <- tryCatch({
    process_site(y_raw, site_name,
                 save_path = save_path, R = 1000, best_cens = best_cens)
  }, error = function(e) {
    warning(paste("Failed for site:", site_name))
    NULL
  })
  if (!is.null(site_result)) {
    results[[site_name]] <- site_result
  }
}

df_all_results <- do.call(rbind, results)
print(df_all_results)
