
# load libraries
source("load_libraries.R")
library(units)

# load global functions
library(generain)

# load data
df_comephore <- read.csv("../data/comephore/inside_mtp.csv", sep = ",")
loc_px <- read.csv("../data/comephore/loc_pixels_mtp.csv", sep = ",")

comephore <- df_comephore[-1] # remove dates column
# Get distances matrix
dist_mat <- get_dist_mat(loc_px)
df_dist <- reshape_distances(dist_mat)


################################################################################
# QUANTILE ---------------------------------------------------------------------
################################################################################

# get a matrix of high quantiles for all pair
q <- 0.99 # quantile
list_count_quant <- quantile_matrix(q, comephore, qlim = TRUE, zeros = TRUE,
                                    count_min = 50) # with removing zeros
quant_mat <- list_count_quant[1][[1]]
count_mat <- list_count_quant[2]

################################################################################
# EXTREMOGRAM
################################################################################

# Temporal chi with spatial lag fixed at 0 -------------------------------------
tmax <- 10
nsites <- length(loc_px$pixel_name)
# Temporal chi with spatial lag fixed at 0
# compute chiplot of every site with itself but lagged in time
start_time <- Sys.time()
chimat_dtlag <- temporal_chi(comephore, quantile = 0.99, tmax = tmax,
                             mean = FALSE)
end_time <- Sys.time()
elapsed_time <- end_time - start_time
print(elapsed_time)

# every chi lagged mean
par(mfrow = c(1, 1))
chi_df_dt <- data.frame(chimat_dtlag)
colnames(chi_df_dt) <- c(1:tmax) # temporal lags from 1 to tmax
rownames(chi_df_dt) <- c(1:nsites) # stations

chimat_dt_mean <- temporal_chi(comephore, tmax, quantile = 0.99, mean = TRUE)
# get h axis in minutes ie x5 minutes
df <- data.frame(lag = c(1:tmax), chi = chimat_dt_mean)
ggplot(df, aes(x = lag, y = chi)) +
  geom_point(color = btfgreen) +
  btf_theme +
  xlab("Temporal lag") +
  ylab(TeX(r"($\hat{\chi}$)"))

# boxplot all stations values for chi temp
# Reshape data using gather function
df_gathered <- chi_df_dt %>% gather(key = "variable", value = "value")
df_gathered$group <- factor(as.integer(df_gathered$variable),
                            levels = seq(1, tmax))

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
                   limits = c(0, 1)) +
  theme(panel.grid.major = element_line(color = "white"),
        panel.grid.minor = element_line(color = "white"))

chitemp

wlse_temp <- get_estimate_variotemp(chimat_dt_mean, tmax, nsites,
                                    weights = "exp", summary = TRUE)

alpha2 <- wlse_temp[[2]]
beta2 <- wlse_temp[[1]]
c2 <- log(beta2)

dftemp <- data.frame(lag = log(df$lag), chi = eta(df$chi))

chitemp_eta_estim <- ggplot(dftemp, aes(x = lag, y = chi)) +
  geom_point(color = btfgreen, size = 4) +
  btf_theme +
  xlab(TeX(r"($\log(\tau)$)")) +
  ylab(TeX(r"($\eta(\widehat{\chi}(0,\tau))$)")) +
  geom_line(aes(x = lag, y = alpha2 * lag + c2),
            alpha = 0.6, color = "darkred", linewidth = 1.5)

chitemp_eta_estim

# Spatial chi ------------------------------------------------------------------

q <- 0.99
chispa_df <- spatial_chi_alldist(df_dist, data_rain = comephore, quantile = q,
                                hmax = 6000, comephore = TRUE)

etachispa_df <- data.frame(chi = eta(chispa_df$chi),
                           lagspa = log(chispa_df$lagspa))

chispa_plot <- ggplot(chispa_df, aes(lagspa, chi)) +
  btf_theme +
  geom_point(col = btfgreen, size = 4) +
  xlab(TeX(r"($h$)")) +
  ylab(TeX(r"($\widehat{\chi}(h, 0)$)")) +
  theme(axis.ticks = element_blank(),
        axis.line = element_blank(),
        panel.border = element_blank(),
        legend.text = element_text(size = 15),
        legend.title = element_text(size = 15),
        legend.position = "right",
        legend.margin = margin(0.5, 0.5, 0.5, 0, "cm"),
        panel.grid = element_line(color =  "#5c595943")) +
  ylim(0, 1)

chispa_plot


chispa_eta <- ggplot(etachispa_df, aes(lagspa, chi)) +
  btf_theme +
  geom_point(col = "darkred", size = 4) +
  xlab(TeX(r"($\log(h)$)")) +
  ylab(TeX(r"($\eta(\widehat{\chi}(h, 0))$)")) +
  theme(axis.ticks = element_blank(),
        axis.line = element_blank(),
        panel.border = element_blank(),
        legend.text = element_text(size = 15),
        legend.title = element_text(size = 15),
        legend.position = "right",
        legend.margin = margin(0.5, 0.5, 0.5, 0, "cm"),
        panel.grid = element_line(color = "#5c595943"))

chispa_eta


wlse_spa <- get_estimate_variospa(chispa_df, weights = "exp", summary = TRUE)

alpha1 <- wlse_spa[[2]]
beta1 <- wlse_spa[[1]]
c1 <- log(beta1)

chispa_eta_estim <- ggplot(etachispa_df, aes(lagspa, chi)) +
  btf_theme +
  geom_point(col = btfgreen, size = 4) +
  xlab(TeX(r"($\log(h)$)")) +
  ylab(TeX(r"($\eta(\widehat{\chi}(h, 0))$)")) +
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

################################################################################
# VARIOGRAM
################################################################################

# estimates
param <- c(beta1, beta2, alpha1, alpha2, 0.001, 0.001)

h_vect <- get_h_vect(df_dist, hmax = 6000)
tau_vect <- 1:10

excesses_com <- empirical_excesses(comephore, quantile = 0.99, tau = tau_vect,
                df_dist = df_dist, h_vect = h_vect)

parscale <- c(1,1,1,1,1,1)  # scale parameters
# get the variogram
result <- optim(par = param, fn = neg_ll,
                  excesses = excesses_com, control = list(parscale = parscale),
                  h_vect = h_vect, tau = tau_vect,
                  df_dist = df_dist,
                  quantile = 0.99)
library(optimx)
lower.bound <- c(1e-6, 1e-6, 1e-6, 1e-6, -Inf, -Inf)
upper.bound <- c(Inf, Inf, 1.999, 1.999, Inf, Inf)
result <- optimr(par = param, method = "Rcgmin",
                  gr = "grfwd", fn = function(par) {
                  neg_ll(par, excesses = excesses_com, quantile = 0.99,
                        h_vect = h_vect, tau = tau_vect, df_dist = df_dist,
                        simu = comephore)
                  }, lower = lower.bound, upper = upper.bound,
                  control = list(parscale = parscale, maxit = 1000))
