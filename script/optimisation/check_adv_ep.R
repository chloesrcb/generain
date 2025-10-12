# Get all files in the folder "R"
functions_folder <- "./R"
files <- list.files(functions_folder, full.names = TRUE)
# # load all functions in files
invisible(lapply(files, function(f) source(f, echo = FALSE)))
library(latex2exp)
library(sf)
library(dplyr)
library(ggplot2)
library(gridExtra)

# Load libraries and set theme
source("./script/load_libraries.R")
source("./script/optimisation/config_com.R")

q <- 0.97
delta <- 30
min_spatial_dist <- 10
adv_filename <- paste0(data_folder, "comephore/adv_estim/advection_results_q",
                       q * 100, "_delta", delta, "_dmin", min_spatial_dist,
                       ".csv")
adv_df <- read.csv(adv_filename, sep = ",")


filename_com_res <- paste(data_folder, "/comephore/optim_results/lalpha/free_eta/combined_optim_results.csv", sep = "")
com_results <- read.csv(filename_com_res)
# remove na
com_results <- na.omit(com_results)
# get estimates from comephore

params_com <- com_results[com_results$q == q*100 &
                           com_results$delta == delta &
                           com_results$dmin == min_spatial_dist, ]
init_params_com <- c(params_com$beta1, params_com$beta2,
                     params_com$alpha1, params_com$alpha2,
                     params_com$eta1, params_com$eta2)

V_episodes <- adv_df %>%
  filter(!is.na(mean_dx_kmh) & !is.na(mean_dy_kmh)) %>%
  select(mean_dx_kmh, mean_dy_kmh)
colnames(V_episodes) <- c("vx", "vy")

adv_x_1 <- (abs(V_episodes$vx)^1) * sign(V_episodes$vx) * 1
adv_y_1 <- (abs(V_episodes$vy)^1) * sign(V_episodes$vy) * 1
head(adv_x_1)
max(adv_x_1)
which.max(adv_x_1)
V_episodes[which.max(adv_x_1), ]

eta1_com <- params_com$eta1
eta2_com <- params_com$eta2
adv_x_com <- (abs(V_episodes$vx)^eta2_com) * sign(V_episodes$vx) * eta1_com
adv_y_com <- (abs(V_episodes$vy)^eta2_com) * sign(V_episodes$vy) * eta1_com
head(adv_x_com)
which.max(adv_x_com)
max(adv_x_com)
V_episodes[which.max(adv_x_com), ]

adv_df_1 <- data.frame(adv_x = adv_x_1, adv_y = adv_y_1)
adv_df_com <- data.frame(adv_x = adv_x_com, adv_y = adv_y_com)
adv_df_1$set <- "Barycentric advection"
adv_df_com$set <- "Transformed advection"
adv_df_all <- rbind(adv_df_1, adv_df_com)

wind_df_plot <- ggplot(adv_df_all, aes(x = 0, y = 0)) +
  geom_segment(aes(xend = adv_x, yend = adv_y), 
               arrow = arrow(length = unit(0.1, "cm")),
               size = 0.5) +
  geom_point(aes(x = adv_x, y = adv_y), size = 1) +
  btf_theme +
  coord_equal() +
  xlab("Advection in x (km/h)") +
  ylab("Advection in y (km/h)") +
  ggtitle("") +
  theme(plot.title = element_text(hjust = 0.5)) +
  facet_wrap(~set)

wind_df_plot

# save plot
foldername <- paste0(data_folder, "comephore/adv_estim/")
if (!dir.exists(foldername)) {
  dir.create(foldername, recursive = TRUE)
}
ggsave(filename = paste0(foldername, "adv_comparison_etas_adv_q", q*100,
                         "_delta", delta,
                         "_dmin", min_spatial_dist, ".png"),
       plot = wind_df_plot,
       width = 8, height = 4)

# remove adv_x and adv_y > 300
# adv_x_1[adv_x_1 > 300] <- NA
# adv_y_1[adv_y_1 > 300] <- NA
adv_x_com[abs(adv_x_com) > 300] <- NA
adv_y_com[abs(adv_y_com) > 300] <- NA

# number of removed rows
nb_outliers <- nrow(adv_df_com) - nrow(adv_df_com[!is.na(adv_x_com) & !is.na(adv_y_com), ])

adv_df_1 <- data.frame(adv_x = adv_x_1, adv_y = adv_y_1)
adv_df_com <- data.frame(adv_x = adv_x_com, adv_y = adv_y_com)
adv_df_1$set <- "Barycentric advection"
adv_df_com$set <- "Transformed advection"
adv_df_all <- rbind(adv_df_1, adv_df_com)

# remove outliers from adv_df_all
adv_df_all <- adv_df_all[!is.na(adv_df_all$adv_x) & !is.na(adv_df_all$adv_y), ]
wind_df_plot <- ggplot(adv_df_all, aes(x = 0, y = 0)) +
  geom_segment(aes(xend = adv_x, yend = adv_y), 
               arrow = arrow(length = unit(0.1, "cm")),
               size = 0.5) +
  geom_point(aes(x = adv_x, y = adv_y), size = 1) +
  btf_theme +
  coord_equal() +
  xlab("Advection in x (km/h)") +
  ylab("Advection in y (km/h)") +
  ggtitle(paste0("Removed outliers: ", nb_outliers)) +
  theme(plot.title = element_text(hjust = 0.5)) +
  facet_wrap(~set) 
wind_df_plot

# save plot
foldername <- paste0(data_folder, "comephore/adv_estim/")
if (!dir.exists(foldername)) {
  dir.create(foldername, recursive = TRUE)
}
ggsave(filename = paste0(foldername, "adv_comparison_etas_adv_q", q*100,
                         "_delta", delta,
                         "_dmin", min_spatial_dist, "_no_outliers.png"),
       plot = wind_df_plot,
       width = 8, height = 4)



