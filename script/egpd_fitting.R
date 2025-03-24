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

load("./data/PluvioMontpellier_1min/rain_mtp_5min_2019_2022.RData")
rain <- rain.all5[c(1, 6:ncol(rain.all5))]
rownames(rain) <- rain$dates
rain_new <- rain[-1] # remove dates column

df_comephore <- read.csv("./data/comephore/inside_HF.csv", sep = ",")
# loc_px <- read.csv("./data/comephore/loc_pixels_HF.csv", sep = ",")

# count number of non na data by column
count_list <- c()
for (i in 1:ncol(rain_new)) {
   count_list[i] <- sum(!is.na(rain_new[,i]))
}

# convert count 5 min list into count year
count_hour <- count_list / 12
count_day <- count_hour / 24
count_year <- count_day / 365

# keep only columns with at least 2.5 year
index_col <- which(count_year > 2.5)
rain_sub <- rain_new[ ,index_col]

# Take only data 3 years of data from COMEPHORE
colnames(df_comephore)[1] <- "date"
df_comephore <- df_comephore[df_comephore$date >= "2018-01-01", ]
length(df_comephore[,1]) / 24 / 365
rain_com <- df_comephore[-1] # remove dates column
ncol(rain_com)


# rain <- rain[rain$dates >= "2018-01-01"]
rownames(rain) <- rain$dates
rain_new <- rain[-1] # remove dates column
colnames(rain_sub)
################################################################################
# EGPD ---------------------------------------------------------------------
################################################################################

# get censoring
censores <- seq(0.2, 5, 0.1)
df_score_com <- choose_censore(rain_com, censores, n_samples = 100)
censores <- seq(0.2, 0.5, 0.01)
df_score_ohsm <- choose_censore(rain_sub, censores, n_samples = 100)

params_com_mtp <- get_egpd_estimates(rain_com,
                            left_censoring = df_score_com$censoreRMSE)

params_subrain <- get_egpd_estimates(rain_sub,
                            left_censoring = 0.3)
df_long_com <- get_df_long_params_egpd(params_com_mtp)
df_long_com$Dataset <- "COMEPHORE"

df_long_ohsm <- get_df_long_params_egpd(params_subrain)
df_long_ohsm$Dataset <- "OHSM"

df_combined <- dplyr::bind_rows(df_long_com, df_long_ohsm)

# Create the boxplot for all three parameters by data source using ggplot
params_boxplot <- ggplot(df_combined, aes(x = Variable, y = Value,
                         fill = Dataset)) +
  geom_boxplot(alpha = 0.8, color = "#0000006c") +
  labs(x = "Parameters", y = "Distribution") +
  btf_theme +
  scale_x_discrete(labels = c(Xi = TeX(r"($\widehat{\xi}$)"),
                              Sigma = TeX(r"($\widehat{\sigma}$)"),
                              Kappa = TeX(r"($\widehat{\kappa}$)"))) +
  scale_fill_manual(values = c("COMEPHORE" = "#a55c5cbe", "OHSM" = btfgreen)) +
  theme(legend.text = element_text(size = 12),
        legend.title = element_text(size = 14)) +
  ylim(0, 1.5)
params_boxplot

ggsave(filename = paste0(im_folder, "EGPD/params_estim_2DS_censoreRMSE.png"),
       plot = params_boxplot, device = "png", bg = "transparent",
       width = 8, height = 6)


# ggsave(filename = paste0(im_folder, "EGPD/OHSM/params_estim_boxplot_maxperiod.png"),
#          plot = com_boxplot, device = "png", bg = "transparent",
#          width = 8, height = 6)

kappa_com <- mean(params_com_mtp$kappa)
sigma_com <- mean(params_com_mtp$sigma)
xi_com <- mean(params_com_mtp$xi)
kappa_ohsm <- mean(params_subrain$kappa)
sigma_ohsm <- mean(params_subrain$sigma)
xi_ohsm <- mean(params_subrain$xi)
