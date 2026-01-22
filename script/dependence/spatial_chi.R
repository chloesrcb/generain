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
                         "omsev/omsev_5min/rain_mtp_5min_2019_2024.RData")

load(filename_omsev)
rain <- rain.all5[c(1, 6:ncol(rain.all5))]
typeof(rain)
class(rain)
colnames(rain)
rownames(rain) <- rain$dates
rain_new <- rain[-1] # remove dates column


# SPATIAL CHI ------------------------------------------------------------------


build_spatial_classes <- function(dist_mat, nb_classes = 8, method = c("quantile", "fixed"), bin_width = 200) {
  method <- match.arg(method)
  
  # Étape 1: Extraire la partie triangulaire supérieure sans diagonale
  dist_vals <- dist_mat[upper.tri(dist_mat, diag = FALSE)]
  dist_vals <- dist_vals[dist_vals > 0]  # par sécurité

  # Étape 2: Créer les bornes des classes
  if (method == "quantile") {
    breaks <- unique(quantile(dist_vals, probs = seq(0, 1, length.out = nb_classes + 1)))
  } else {
    breaks <- seq(from = min(dist_vals), to = max(dist_vals) + bin_width, by = bin_width)
  }

  # Étape 3: Créer la matrice rad_mat avec les classes
  rad_mat <- matrix(NA, nrow = nrow(dist_mat), ncol = ncol(dist_mat))

  for (i in 2:length(breaks)) {
    # Assigne la classe uniquement à la partie supérieure
    in_bin <- dist_mat >= breaks[i - 1] & dist_mat < breaks[i] & upper.tri(dist_mat)
    rad_mat[in_bin] <- breaks[i]  # ou (breaks[i-1]+breaks[i])/2
  }

  return(list(
    rad_mat = rad_mat,
    breaks = breaks
  ))
}

res <- build_spatial_classes(dist_mat, nb_classes = 10)
rad_mat <- res$rad_mat
breaks <- res$breaks
barplot(table(rad_mat), main = "Nombre de paires par classe de distance", xlab = "Lag (m)", ylab = "Count")


# Get dataframe for the histogram plot
df_hist <- data.frame(
  Interval = names(table(rad_mat)),
  Count = as.numeric(table(rad_mat))
)

df_hist$lower_bounds <- breaks[-length(breaks)]
df_hist$upper_bounds <- breaks[-1]
df_hist$Interval_str <- paste0(round(df_hist$lower_bounds,0), " - ", round(df_hist$upper_bounds,0))

# # Histogram
histradius <- ggplot(df_hist, aes(x = Interval, y = Count)) +
  btf_theme +
  geom_bar(stat = "identity", fill = btfgreen, alpha = 0.8) +
  xlab("Spatial lag") +
  ylab("Pair count") +
  theme(axis.text.x = element_text(angle = 45)) +
  scale_x_discrete(labels = df_hist$Interval_str) +
  scale_y_continuous(breaks = c(0, 4, 6, 8, 10, 12))

histradius




# remove brives, hydro and cines from rain_new
rain_new <- rain_new[, -c(18, 19, 20)] # remove brives, hydro and cines
chi_by_pair <- compute_all_pairs_chi(rain_new, dist_mat, q = 0.96, zeros = FALSE)
head(chi_by_pair)

# plot chi for each distance
plot <- ggplot(chi_by_pair, aes(x = distance, y = chi)) +
  geom_point() +
  xlab("Distance (m)") +
  ylab("Chi value") +
  ggtitle("Chi values for each pair of stations") +
  btf_theme

plot

# Exemple avec des bins réguliers


agg_chi <- aggregate_chi_by_distance(chi_by_pair, nb_bins = 12, method = "quantile")

# Affichage
print(agg_chi)

# Graphique
plot(agg_chi$distance_mid, agg_chi$chi,
     type = "b", pch = 19, col = "darkgreen",
     xlab = "Distance (bin midpoint)", ylab = "Mean Chi")



q <- 0.96 # quantile
# plot chi for each distance
# rain_new <- rain
chispa_df <- spatial_chi(rad_mat, rain_new,
                         q = q, zeros = F)

chispa_df2 <- spatial_chi(rad_mat, rain_new,
                         q = q, zeros = F, breaks = breaks, mid = T)


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

chispa_df$lagspa <- chispa_df$lagspa / 1000 # convert to km
wlse_spa <- get_estimate_variospa(chispa_df, weights = "exp", summary = TRUE)
c1 <- as.numeric(wlse_spa[[1]])
beta1 <- as.numeric(wlse_spa[[2]])
alpha1 <- as.numeric(wlse_spa[[3]])
signif_beta <- wlse_spa[[4]]
signif_alpha <- wlse_spa[[5]]

etachispa_df <- data.frame(chi = zeta(chispa_df$chi),
                           lagspa = log(chispa_df$lagspa))


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
  ylim(-2, 0)

chispa_eta_estim



quantiles <- seq(0.9, 0.99, by = 0.01)
results_spa <- data.frame()

for (q in quantiles) {
  # Spatial chi for each quantile
  print(q)
  chispa_df <- spatial_chi(rad_mat, rain_new, breaks = breaks, quantile = q, zeros = F, mid = T)
  # Remove lag 0
#   chispa_df <- chispa_df[chispa_df$lagspa != 0, ]
  chispa_df$lagspa <- chispa_df$lagspa / 1000# convert to km
  # Estimate parameters using WLSE
  wlse_spa <- get_estimate_variospa(chispa_df, weights = "exp", summary = T)
  c1 <- as.numeric(wlse_spa[[1]])
  beta1 <- as.numeric(wlse_spa[[2]])
  alpha1 <- as.numeric(wlse_spa[[3]])
  # signif_beta <- wlse_spa[[4]]
  # signif_alpha <- wlse_spa[[5]]
  
  # results_spa <- rbind(results_spa,
  #                      data.frame(quantile = q,
  #                                 beta = beta1,
  #                                 alpha = alpha1,
  #                                 signif_beta = signif_beta,
  #                                 signif_alpha = signif_alpha))
  
  # Add transformation
  etachispa_df <- data.frame(
    chi = zeta(chispa_df$chi),
    lagspa = log(chispa_df$lagspa)
  )
  
  # Plot points + regression
  p <- ggplot(etachispa_df, aes(x = lagspa, y = chi)) +
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
    geom_line(aes(y = alpha1 * lagspa + c1),
              alpha = 0.6, color = "darkred", size = 1.5) +
    ggtitle(paste0("Quantile = ", q))
  
  # save plot
  # filename <- paste(im_folder, "WLSE/omsev/2022/spatial_chi_",
  #                       as.character(q), ".pdf", sep = "")
  # ggsave(filename, plot = p, width = 8, height = 6, dpi = 300)
}


print(results_spa)



# total count of exedences for all stations

chi_slag <- c()
chi_val <- NA
par(mfrow = c(2, 3))
q <- 0.998
step_rad <- 100 
hmax <- 1600
# TODO
for (h in seq(step_rad, hmax, step_rad)){
  # station number inside h lag
  indices <- data.frame(which(rad_mat == h, arr.ind = TRUE))
  # count_Nh <- sum(df_rad$value == h)
  print(paste0("h = ", h))
  nb_pairs <- dim(indices)[1]
  if (nb_pairs == 0) {
    chi_slag <- c(chi_slag, NA)
  } else {
    # get index pairs
    ind_s1 <- indices$row
    ind_s2 <- indices$col
    # mutual_excess <- c()
    # quants <- c()
    # q_vals <- c()
    chi_val <- c()
    for (i in 1:nb_pairs){
      rain_cp <- drop_na(rain_new[, c(ind_s1[i], ind_s2[i])])
      # # colnames(rain_cp)
      # quant <- quant_mat[ind_s1[i], ind_s2[i]]
      # q_val <- get_quantile_level(quant, rain_cp[, 1], rain_cp[, 2])
      # count <- sum(rain_cp[, 1] > q_val & rain_cp[, 2] > q_val)
      # count <- count_excess[ind_s1[i], ind_s2[i]]
      # q_vals <- c(q_vals, q_val)
      # quants <- c(quants, quant)
      # mutual_excess <- c(mutual_excess, count)
      # qut <- quant_mat[ind_s1[i], ind_s2[i]]
      cp_dt <- chiplot(rain_cp[], xlim=c(0.995, 0.9999), qlim = c(0.99, 0.9997),
                       which = 1,
                       main1=paste0("s", ind_s1[i], ", s", ind_s2[i], ", lag ", h))
      abline(v=q)
      quantile <- which(cp_dt$quantile > q)[1]
      count <- sum(rain_cp[, 1] > q & rain_cp[, 2] > q)
      print(count)
      chi_cp_dt <- cp_dt$chi[, 2]
      chi_val <- c(chi_val, chi_cp_dt[quantile])
      
      
    }
    chi_slag <- c(chi_slag, mean(chi_val))
  }
  print(chi_slag)
}












# for multiple quantiles
qs <- seq(0.9, 0.99, by = 0.01)
hmax <- 10  # par exemple

# Initialiser une liste pour stocker les résultats
chi_list <- lapply(qs, function(q) {
  chi_vals <- spatial_chi(rad_mat, rain_new, quantile = q, zeros = FALSE, mid = T)
#   chi_vals <- chi_vals[-1]  # enlever lag 0
  data.frame(
    lag = chi_vals$lagspa / 1000,  # convertir en km
    chi = chi_vals$chi,
    q = q
  )
})

# Combiner tous les data.frames en un seul
chi_df <- bind_rows(chi_list)

# Tracer
chispa_multi <- ggplot(chi_df, aes(x = lag, y = chi, color = factor(q))) +
  geom_line(alpha = 1, linewidth = 1.2) +
xlab(TeX(r"($\|h\|$ (km))")) +
  ylab(TeX(r"($\widehat{\chi}(h,0)$)")) +
  ylim(0, 1) +
  scale_color_brewer(palette = "BuPu", name = "Quantile") +
  theme_minimal() +
  theme(panel.grid.major = element_line(color = "white"),
        panel.grid.minor = element_line(color = "white")) +
  btf_theme

chispa_multi

# Save the plot
filename <- paste(im_folder, "WLSE/omsev/2022/spatial_chi_12intervals_multi_q_zeros.pdf", sep = "")
ggsave(filename, width = 20, height = 15, units = "cm")






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


# WLSE transformation ----------------------------------------------------------



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


# WLSE transformation ----------------------------------------------------------
