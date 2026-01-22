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
library(boot)
library(spdep)

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

################################################################################

# rain.all5 is the data frame name for 5 min data
head(rain)
rownames(rain) <- rain$dates
rain$dates <- as.Date(rain$dates)
library(dplyr)
rain <- rain[-1]  # remove dates column
rain.wo.na <- drop_na(rain[c(1, which(colnames(rain) == "cefe"))])
y <- rain.wo.na$cnrs[rain.wo.na$cnrs > 0]

# Precision 
# min for 1 tipping = 0.215252
rain_chu1_pos <- rain$chu1[rain$chu1 > 0]
min_precision <- min(rain_chu1_pos, na.rm = TRUE)
sort(unique(rain_chu1_pos))


results <- list()
sites_name <- colnames(rain)

# Get possible censoring values according to precision
min_censuring <- round(min_precision, 4)
k_prec_censuring <- seq(min_censuring, min_censuring * 4, by = min_censuring)

save_path <- paste0(im_folder,"/EGPD/OMSEV/2019_2024/censoring_022_to_060/")
if (!dir.exists(save_path)) {
  dir.create(save_path, recursive = TRUE)
}
# Fit EGPD for each site
for (site_name in sites_name) {
  cat("Processing:", site_name, "\n")
  y_raw <- as.data.frame(na.omit(rain[[site_name]]))
  site_result <- tryCatch({
    process_site(y = y_raw, site_name =  site_name, 
                 save_path = save_path, R = 1000, forced_cens = k_prec_censuring)
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
# save results
write.csv(df_all_results, file = paste0(im_folder,"/EGPD/OMSEV/2019_2024/egpd_results.csv"), 
            row.names = FALSE)

# remove brives, hydro, cines
df_all_results <- df_all_results %>%
  filter(!Site %in% c("brives", "hydro", "cines"))

df_plot <- df_all_results %>% arrange(xi) 

ggplot(df_plot, aes(x = reorder(Site, xi), y = xi)) +
  geom_point(size = 3, color = "#0072B2") +
  geom_errorbar(aes(ymin = xi_low, ymax = xi_high), width = 0.2, color = "#0072B2") +
  coord_flip() +
  theme_minimal() +
  xlab("Site") + ylab("Parameter") +
  btf_theme 
foldername <- paste0("./thesis/resources/images/EGPD/OMSEV/")
ggsave(paste0(foldername, "xi_by_station.png"), width = 8, height = 6, dpi = 400)


ggplot(df_plot, aes(x = reorder(Site, sigma), y = sigma)) +
  geom_point(size = 3, color = "#009E73") +
  geom_errorbar(aes(ymin = sigma_low, ymax = sigma_high), width = 0.2, color = "#009E73") +
  coord_flip() +
  theme_minimal() +
  xlab("Site") + ylab("Parameter") +
  btf_theme 
 ggsave(paste0(foldername, "sigma_by_station.png"), width = 8, height = 6, dpi = 400)


ggplot(df_plot, aes(x = reorder(Site, kappa), y = kappa)) +
  geom_point(size = 3, color = "#D55E00") +
  geom_errorbar(aes(ymin = kappa_low, ymax = kappa_high), width = 0.2, color = "#D55E00") +
  coord_flip() +
  theme_minimal() +
    xlab("Site") + ylab("Parameter") +
  btf_theme 
ggsave(paste0(foldername, "kappa_by_station.png"), width = 8, height = 6, dpi = 400)

################################################################################
# Statistical tests -----------------------------------------------------------
################################################################################

library(spdep)
coords <- cbind(location_gauges$Longitude, location_gauges$Latitude)
# around 5 years
sites_name_5y <- c("cefe", "cnrs", "crbm", "iem", "mse",
                   "poly", "um", "archiw")

# around 3 years
sites_name_3y <- c("archie", "chu1", "chu2", "chu3",
                   "chu4", "chu5", "chu6", "chu7",
                   "um35")

# around 2 years
sites_name_2y <- c("cines", "hydro", "brives")

# add duration column
df_all_results$duration <- dplyr::case_when(
  df_all_results$Site %in% sites_name_5y ~ "5 ans",
  df_all_results$Site %in% sites_name_3y ~ "3 ans",
  df_all_results$Site %in% sites_name_2y ~ "2 ans",
  TRUE ~ "Inconnue"
)

# remove brives, hydro, cines
df_all_results <- df_all_results %>%
  filter(!Site %in% c("brives", "hydro", "cines"))

# do kruskal wallis test -------------------------------------------------
kruskal.test(xi ~ duration, data = df_all_results)
kruskal.test(sigma ~ duration, data = df_all_results)
kruskal.test(kappa ~ duration, data = df_all_results)

# Adjust parameters by duration means
df_results <- df_all_results %>%
  group_by(duration) %>%
  mutate(
    xi_adj = xi - mean(xi, na.rm = TRUE),
    sigma_adj = sigma - mean(sigma, na.rm = TRUE),
    kappa_adj = kappa - mean(kappa, na.rm = TRUE)
  ) %>% ungroup()

colnames(location_gauges) <- c("Site", "Longitude", "Latitude")


# Do Moran test -----------------------------------------------------------
df_geo <- left_join(df_results, location_gauges, by = "Site")
coords_mat <- as.matrix(df_geo[, c("Longitude", "Latitude")])
nb <- knn2nb(knearneigh(coords_mat, k = 4))   # k = 4 nearest neighbors
lw <- nb2listw(nb, style = "W")
moran.test(df_geo$xi_adj, lw)
moran.test(df_geo$sigma_adj, lw)
moran.test(df_geo$kappa_adj, lw)

# Kolmogorov-Smirnov test ---------------------------------------------------
rain_mat <- rain[, sapply(rain, is.numeric)]
sites <- colnames(rain_mat)
# remove brives, hydro, cines
rain_mat <- rain_mat %>%
  select(-c(brives, hydro, cines))
sites <- colnames(rain_mat)
rain_long <- stack(rain_mat)
names(rain_long) <- c("Rain", "Site")
rain_long <- subset(rain_long, !is.na(Rain) & Rain > 0)

ks_results <- combn(sites, 2, function(pair) {
  site1 <- rain_long$Rain[rain_long$Site == pair[1]]
  site2 <- rain_long$Rain[rain_long$Site == pair[2]]
  
  test <- suppressWarnings(ks.test(site1, site2))
  
  data.frame(
    Site1 = pair[1],
    Site2 = pair[2],
    p.value = test$p.value,
    statistic = test$statistic
  )
}, simplify = FALSE)

ks_df <- bind_rows(ks_results)
head(ks_df)

mean(ks_df$p.value > 0.05)










library(dplyr)
library(tidyr)

# Recharger les données avec dates si pas déjà fait
rain_df <- rain.all5 %>%
  select(dates, all_of(df_all_results$Site)) %>%
  mutate(dates = as.Date(dates))

# Garder seulement les dates présentes dans toutes les stations
rain_sync <- rain_df %>%
  filter(rowSums(!is.na(select(., -dates))) == ncol(.) - 1)

# Long format
rain_long_sync <- rain_sync %>%
  pivot_longer(-dates, names_to = "Site", values_to = "Rain") %>%
  filter(Rain > 0)

# KS synchronisé
sites_sync <- unique(rain_long_sync$Site)

ks_sync_results <- combn(sites_sync, 2, function(pair) {
  x1 <- rain_long_sync$Rain[rain_long_sync$Site == pair[1]]
  x2 <- rain_long_sync$Rain[rain_long_sync$Site == pair[2]]
  suppressWarnings(ks.test(x1, x2))$p.value
})

mean(ks_sync_results > 0.05)



library(dplyr)
library(tidyr)
library(kSamples)
library(twosamples)

# 1) Séries synchronisées sur la période commune à toutes les stations
rain_df <- rain.all5 %>%
  select(dates, all_of(sites)) %>%
  mutate(dates = as.Date(dates))

rain_sync <- rain_df 
rain_long_sync <- rain_sync %>%
  pivot_longer(-dates, names_to = "Site", values_to = "Rain") %>%
  filter(Rain > 0)

sites_sync <- unique(rain_long_sync$Site)

# 2) CvM en paires
u <- quantile(rain_long_sync$Rain, 0.9)  # exemple

cvm_results <- combn(sites_sync, 2, function(pair) {
  x1 <- rain_long_sync$Rain[rain_long_sync$Site == pair[1]]
  x2 <- rain_long_sync$Rain[rain_long_sync$Site == pair[2]]
  x1_tail <- x1[x1 > u]
  x2_tail <- x2[x2 > u]

  cvmtest <- twosamples::cvm_test(x1_tail, x2_tail, nboots = 10000)
  adtest <- twosamples::ad_test(x1_tail, x2_tail, nboots = 20000)
  data.frame(
    Site1 = pair[1],
    Site2 = pair[2],
    cvm_statistic = cvmtest[1],
    cvm_p.value = cvmtest[2],
    ad_statistic = adtest[1],
    ad_p.value = adtest[2]
  )
}, simplify = FALSE) |> bind_rows()

mean_cvm_ok <- mean(cvm_results$cvm_p.value > 0.05)
mean_cvm_ok

mean_ad_ok <- mean(cvm_results$ad_p.value > 0.05)
mean_ad_ok

x1 <- rain_long_sync$Rain[rain_long_sync$Site == "cefe"]
x2 <- rain_long_sync$Rain[rain_long_sync$Site == "cnrs"]
test <- twosamples::cvm_test(x1, x2, nboots = 10000)
test[1]


u <- 0.65

x1_tail <- x1[x1 > u]
x2_tail <- x2[x2 > u]

twosamples::ad_test(x1_tail, x2_tail, nboots = 20000)
twosamples::cvm_test(x1_tail, x2_tail, nboots = 10000)
