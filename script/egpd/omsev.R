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
# save in csv
# filename_omsev <- paste0(data_folder,
#                          "omsev/omsev_5min/rain_mtp_5min_2019_2025.csv")
# # write.csv(rain, filename_omsev, row.names = FALSE)
# rain <- read.csv(filename_omsev)
# head(rain)
# tail(rain)

# rain_tilljan2025 <- rain[rain$dates < as.POSIXct("2025-02-01"), ] 
# tail(rain_tilljan2025)
filename_omsev <- paste0(data_folder,
                         "omsev/omsev_5min/rain_mtp_5min_2019_2024.RData")

load(filename_omsev)
rain <- rain.all_save
# rain.all5 <- rain.all5[, c(1, 6:(ncol(rain.all5)))]

# compare rain.all5 and rain dataframes on same period
rain$dates <- as.POSIXct(rain$dates)

# write.csv(rain_tilljan2025, filename_omsev, row.names = FALSE)
# rain <- read.csv(filename_omsev)
# # head(rain)
tail(rain)
# keep only dates before january 2025
# rain <- rain[rain$dates < as.POSIXct("2025-11-30"), ]
rownames(rain) <- rain$dates
# remove rain$dates
rain$dates <- NULL

################################################################################

# rain.wo.na <- drop_na(rain[c(1, which(colnames(rain) == "cefe"))])
# y <- rain.wo.na$cnrs[rain.wo.na$cnrs > 0]


rain_chu1_pos <- rain$cnrs[rain$cnrs > 0]
min_precision <- min(rain_chu1_pos, na.rm = TRUE)
sort(unique(rain_chu1_pos))


results <- list()
sites_name <- colnames(rain)

# remove year 2024
# rain <- rain[!grepl("2025", rownames(rain)), ]
# check if there are any NA values in the data

# Get possible censoring values according to precision
min_censuring <- round(min_precision, 2)
k_prec_censuring <- seq(min_censuring, min_censuring * 4, by = min_censuring)
# k_prec_censuring <- seq(0.22, 0.44, by = 0.01)
# p <- min_censuring
# k_prec_censuring <- c(p, 1.5 * p, 2 * p, 2.5 * p, 3 * p)

save_path <- paste0(im_folder,"/EGPD/OMSEV/2019_2024/precision_new/")

if (!dir.exists(save_path)) {
  dir.create(save_path, recursive = TRUE)
}
# Fit EGPD for each site
for (site_name in sites_name[-1]) {
  cat("Processing:", site_name, "\n")
  y <- as.data.frame(na.omit(rain[[site_name]]))
  site_result <- tryCatch({
    process_site(y = y, site_name =  site_name,
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
# write.csv(df_all_results, file = paste0(im_folder,"/EGPD/OMSEV/2019_2025/egpd_results.csv"), 
#             row.names = FALSE)

df_all_results <- read.csv(paste0(im_folder,"/EGPD/OMSEV/2019_2024/egpd_results.csv"))










# remove brives, hydro, cines
df_all_results <- df_all_results %>%
  filter(!Site %in% c("brives", "hydro", "cines"))

df_plot <- df_all_results %>% arrange(xi) 

foldername <- paste0(im_folder,"/EGPD/OMSEV/2019_2025/")
ggplot(df_plot, aes(x = reorder(Site, kappa), y = xi)) +
  geom_point(size = 3, color = "#0072B2") +
  geom_errorbar(aes(ymin = xi_low, ymax = xi_high), width = 0.2, color = "#0072B2") +
  coord_flip() +
  theme_minimal() +
  xlab("Site") + ylab("Parameter") +
  btf_theme 
ggsave(paste0(save_path, "xi_by_station.png"), width = 8, height = 6, dpi = 400)


ggplot(df_plot, aes(x = reorder(Site, kappa), y = sigma)) +
  geom_point(size = 3, color = "#009E73") +
  geom_errorbar(aes(ymin = sigma_low, ymax = sigma_high), width = 0.2, color = "#009E73") +
  coord_flip() +
  theme_minimal() +
  xlab("Site") + ylab("Parameter") +
  btf_theme 
 ggsave(paste0(save_path, "sigma_by_station.png"), width = 8, height = 6, dpi = 400)


ggplot(df_plot, aes(x = reorder(Site, kappa), y = kappa)) +
  geom_point(size = 3, color = "#D55E00") +
  geom_errorbar(aes(ymin = kappa_low, ymax = kappa_high), width = 0.2, color = "#D55E00") +
  coord_flip() +
  theme_minimal() +
    xlab("Site") + ylab("Parameter") +
  btf_theme 
ggsave(paste0(save_path, "kappa_by_station.png"), width = 8, height = 6, dpi = 400)


# Do boxplots -------------------------------------------------------
df_melted <- df_all_results %>%
  select(Site, xi, sigma, kappa) %>%
  pivot_longer(cols = c(xi, sigma, kappa), names_to = "Parameter", values_to = "Value")
latex_names <- c(
    xi = TeX("$\\widehat{\\xi}$"),
    sigma = TeX("$\\widehat{\\sigma}$"),
    kappa = TeX("$\\widehat{\\kappa}$")
  )
ggplot(df_melted, aes(x = Parameter, y = Value)) +
  geom_boxplot(fill=btfgreen, alpha = 0.7) +
  theme_minimal() +
  xlab("Parameter") + ylab("Value") +
  scale_x_discrete(labels = latex_names) +
  theme(legend.position = "none") +
  ylim(0, 1.5) +
  btf_theme 
  
ggsave(paste0(save_path, "boxplot_parameters.pdf"), width = 8, height = 6, dpi = 400)
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

location_gauges <- location_gauges %>%
  select(codestation, Longitude, Latitude)
colnames(location_gauges) <- c("Site", "Longitude", "Latitude")


# Do Moran test -----------------------------------------------------------
df_geo <- left_join(df_results, location_gauges, by = "Site")
coords_mat <- as.matrix(df_geo[, c("Longitude", "Latitude")])
nb <- knn2nb(knearneigh(coords_mat, k = 4))   # k = 4 nearest neighbors
lw <- nb2listw(nb, style = "W")
moran.test(df_geo$xi_adj, lw)
moran.test(df_geo$sigma_adj, lw)
moran.test(df_geo$kappa_adj, lw)
moran.test(df_geo$kappa_adj, lw, alternative = "two.sided")
moran.test(df_geo$xi_adj, lw, alternative = "two.sided")
moran.test(df_geo$sigma_adj, lw, alternative = "two.sided")


moran.plot(df_geo$sigma_adj, lw)
moran.plot(df_geo$xi_adj, lw)

ggplot(df_geo) +
  geom_point(aes(Longitude, Latitude, color = xi_adj), size = 4) +
  scale_color_viridis_c() +
  coord_fixed() 

# Kolmogorov-Smirnov test ---------------------------------------------------
rain_mat <- rain[, sapply(rain, is.numeric)]
sites <- colnames(rain_mat)
# remove brives, hydro, cines
rain_mat <- rain_mat %>%
  select(-c(brives, hydro, cines))
sites <- colnames(rain_mat)
rain_long <- stack(rain_mat)
names(rain_long) <- c("Rain", "Site")
rain_long <- subset(rain_long, !is.na(Rain))

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


library(kSamples)


# liste des échantillons par site
samples <- split(rain_long$Rain, rain_long$Site)

# test AD k-samples
ad_res <- ad.test(samples)

ad_res

# on positive rain only
samples_pos <- lapply(samples, function(x) x[x > 2])
ad_res_pos <- ad.test(samples_pos)

ad_res_pos