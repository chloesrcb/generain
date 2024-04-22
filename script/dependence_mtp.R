library(generain)
setwd("./script")
# load libraries
source("load_libraries.R")

################################################################################
# LOCATION ---------------------------------------------------------------------
################################################################################
# get location of each rain gauge
location_gauges <- read.csv("../data/PluvioMontpellier_1min/pluvio_mtp_loc.csv")
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
load("../data/PluvioMontpellier_1min/rain_mtp_5min_2019_2022.RData")
rain <- rain.all5[c(1, 6:ncol(rain.all5))]

################################################################################
# QUANTILE ---------------------------------------------------------------------
################################################################################

# get a matrix of high quantiles for all pair
q <- 0.97 # quantile
rownames(rain) <- rain$dates
rain_new <- rain[-1] # remove dates column
list_count_quant <- quantile_matrix(q, rain_new, qlim = FALSE, zeros = FALSE,
                                    count_min = 50) # with removing zeros
quant_mat <- list_count_quant[1][[1]]
count_mat <- list_count_quant[2]

# 0.998 equiv to 0.96 without 0

################################################################################
# EXTREMOGRAM ------------------------------------------------------------------
################################################################################
# Spatial and temporal CHIPLOT and CHI value for all pairs

# TEMPORAL CHI -----------------------------------------------------------------

tmax <- 10
# Temporal chi with spatial lag fixed at 0
# compute chiplot of every site with itself but lagged in time
start_time <- Sys.time()
chimat_dtlag <- temporal_chi(rain_new, tmax, quantile = 0.998, mean = TRUE)
end_time <- Sys.time()
elapsed_time <- end_time - start_time
print(elapsed_time)

# every chi lagged mean
par(mfrow = c(1, 1))
colnames(chi_df_dt) <- c(1:tmax) # temporal lags from 1 to tmax
rownames(chi_df_dt) <- c(1:17) # stations

# get h axis in minutes ie x5 minutes
df <- data.frame(lag = c(1:tmax) * 5, chi = chi_df_dt)
ggplot(df, aes(x = lag, y = chi)) +
  geom_point(color = btfgreen) +
  btf_theme +
  xlab("Temporal lag") +
  ylab(TeX(r"($\hat{\chi}$)"))

# boxplot all stations values for chi temp
# Reshape data using gather function
df_gathered <- chi_df_dt %>% gather(key = "variable", value = "value")
df_gathered$group <- factor(as.integer(df_gathered$variable),
                            levels = seq(1, 20))

