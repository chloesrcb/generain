# remove objects from the environment
rm(list = ls())

# clear console
cat("\014")

# Load libraries and set theme
source("./script/load_libraries.R")

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

load("workspace.RData")
# library(generain)

# LOAD DATA ####################################################################
filename_com <- paste0(data_folder, "comephore/comephore_2008_2024.csv")
# filename_com <- paste0(data_folder, "comephore/zoom_5km.csv")
comephore_raw <- read.csv(filename_com, sep = ",")
filename_loc <- paste0(data_folder, "comephore/loc_px_zoom_5km.csv")
loc_px <- read.csv(filename_loc, sep = ",")

# remove p43 from loc_px
loc_px <- loc_px[loc_px$pixel_name != "p43", ]

df_comephore <- as.data.frame(comephore_raw)
df_comephore$datetime <- as.POSIXct(df_comephore$datetime,
                                format = "%Y-%m-%d %H:%M:%S", tz = "UTC")
head(df_comephore)

# Take only data after 2007
colnames(df_comephore)[1] <- "date"
df_comephore <- df_comephore[df_comephore$date >= "2008-01-01", ]
tail(df_comephore)
# put date in index
rownames(df_comephore) <- df_comephore$date
comephore <- df_comephore[-1] # remove dates column


# # get wind data
filename_wind <- paste0(data_folder, "wind/data_gouv/wind_mtp.csv")
wind_mtp <- read.csv(filename_wind)

# if datetime is just a date put 00:00:00 as time
if (nchar(wind_mtp$datetime[1]) == 10) {
  wind_mtp$datetime <- paste(wind_mtp$datetime, "00:00:00")
}


# Convert datetime to POSIXct
wind_mtp$datetime <- as.POSIXct(wind_mtp$datetime,
                                format = "%Y-%m-%d %H:%M:%S", tz = "UTC")

# Get only data after 2008
wind_mtp <- wind_mtp[wind_mtp$datetime >= "2008-01-01", ]

# Apply function to the DD column
wind_mtp$cardDir <- sapply(wind_mtp$DD, convert_to_cardinal, nb_cardinal = 8)
wind_mtp$cardDir <- as.character(wind_mtp$cardDir)  # Ensure it's character
wind_mtp$cardDir[is.na(wind_mtp$DD)] <- NA
# Check if NA values are properly handled
summary(wind_mtp)
plot(wind_mtp$FF)
# 1 nd = 0,514 m/s
head(wind_mtp$cardDir)

# get year 2024
# wind_mtp <- wind_mtp[wind_mtp$datetime >= "2024-01-01", ]

# FF are in m/s 1/10 so we convert to m/s
# wind_mtp$FF <- wind_mtp$FF 
wind_mtp$FF_kmh <- wind_mtp$FF * 3.6  # convert to km/h

plot(wind_mtp$FF_kmh, type = "l", col = "blue",
     main = "Wind speed in km/h", xlab = "Time", ylab = "Wind speed (km/h)")

# DISTANCE AND COORDS ##########################################################

# Get distances matrix
# dist_mat <- get_dist_mat(loc_px)

# Get number of sites
nsites <- nrow(loc_px)

# Get coords
sites_coords <- loc_px[, c("Longitude", "Latitude")]
rownames(sites_coords) <- loc_px$pixel_name
sites_coords_sf <- st_as_sf(sites_coords, coords = c("Longitude", "Latitude"),
                            crs = 4326)

sites_coords_sf <- st_transform(sites_coords_sf, crs = 2154)

coords_m <- st_coordinates(sites_coords_sf)

grid_coords_km <- sites_coords
grid_coords_m <- sites_coords
grid_coords_m$x_m <- (coords_m[, "X"] - min(coords_m[, "X"]))
grid_coords_m$y_m <- (coords_m[, "Y"] - min(coords_m[, "Y"]))
grid_coords_km$x_km <- (coords_m[, "X"] - min(coords_m[, "X"])) / 1000
grid_coords_km$y_km <- (coords_m[, "Y"] - min(coords_m[, "Y"]))  / 1000




p1 <- ggplot(grid_coords_km, aes(x = Longitude, y = Latitude)) +
  geom_point(color = "blue", size = 2) +
  geom_text(aes(label = rownames(grid_coords_km)), hjust = -0.2, size = 1.5) +
  coord_fixed() +
  ggtitle("GPS coordinates WGS84") +
  theme_minimal()

p2 <- ggplot(grid_coords_km, aes(x = x_km, y = y_km)) +
  geom_point(color = "red") +
  geom_text(aes(label = rownames(grid_coords_km)), size = 1.5, hjust = -0.2) +
  coord_fixed() +
  theme_minimal() +
  ggtitle("Transformed coordinates in km")

p_coords <- grid.arrange(p1, p2, ncol = 2)

# save plot
# filename <- paste(im_folder, "optim/comephore/coords_transformation.png",
#                   sep = "")
# ggsave(plot = p_coords, filename = filename, width = 20, height = 15,
#        units = "cm")

# get distance matrix
grid_coords_m <- grid_coords_m[, c("x_m", "y_m")]
grid_coords_km <- grid_coords_km[, c("x_km", "y_km")]
colnames(grid_coords_m) <- c("Longitude", "Latitude")
colnames(grid_coords_km) <- c("Longitude", "Latitude")
dist_mat <- get_dist_mat(grid_coords_m, latlon = FALSE)

# Spatial chi
df_dist <- reshape_distances(dist_mat)
# df_dist$value <- round(df_dist$value / 1000, 1) * 1000  # / 1000 in km
df_dist_km <- df_dist
df_dist_km$value <- round(df_dist$value / 1000, 1)


# Spatial chi WLSE #############################################################
h_vect <- sort(unique(df_dist_km$value))
h_vect <- h_vect[h_vect > 0]  # remove 0
hmax <- h_vect[10] # 10th value

hmax <- 7

q_no0_spa <- 0.99
chispa_df <- spatial_chi_alldist(df_dist_km, data_rain = comephore,
                            quantile = q_no0_spa, hmax = hmax, zeros = FALSE)

# chispa_df <- spatial_chi_extremogram(df_dist_km, comephore, q_no0_spa,
#                         hmax = hmax, zeros = FALSE)

etachispa_df <- data.frame(chi = eta(chispa_df$chi),
                           lagspa = log(chispa_df$lagspa))


chispa_plot <- ggplot(chispa_df, aes(lagspa, chi)) +
  btf_theme +
  geom_point(col = btfgreen) +
  xlab(TeX(r"($h$)")) +
  ylab(TeX(r"($\widehat{\chi}(h, 0)$)")) +
  ylim(0, 1)

chispa_plot

# save plot
filename <- paste(im_folder, "WLSE/comephore/full_spatial_chi_", q_no0_spa,
                  ".png", sep = "")
ggsave(filename, width = 20, height = 15, units = "cm")

# WLSE
wlse_spa <- get_estimate_variospa(chispa_df, weights = "exp", summary = TRUE)
print(wlse_spa)

c1 <- wlse_spa[[1]]
beta1 <- wlse_spa[[2]]
alpha1 <- wlse_spa[[3]]

chispa_eta_estim <- ggplot(etachispa_df, aes(lagspa, chi)) +
  btf_theme +
  geom_point(col = btfgreen) +
  xlab(TeX(r"($\log(h)$)")) +
  ylab(TeX(r"($\zeta(\widehat{\chi}(h, 0))$)")) +
  geom_line(aes(x = lagspa, y = alpha1 * lagspa + c1), alpha = 0.6,
            color = "darkred", linewidth = 0.5)

chispa_eta_estim

# save plot
filename <- paste(im_folder, "WLSE/comephore/full_spatial_chi_eta_estim_", q_no0_spa,
                  ".png", sep = "")
ggsave(filename, width = 20, height = 15, units = "cm")


# Temporal chi WLSE ############################################################

tmax <- 10
q_no0_temp <- 0.98 # quantile for temporal chi

# Temporal chi with spatial lag fixed at 0
# compute chiplot of every site with itself but lagged in time
chimat_dtlag <- temporal_chi(comephore, quantile = q_no0_temp, tmax = tmax,
                             mean = FALSE, zeros = T)

chi_df_dt <- data.frame(chimat_dtlag)
colnames(chi_df_dt) <- c(0:tmax) # temporal lags from 1 to tmax
rownames(chi_df_dt) <- c(1:nsites) # stations
# remove lag 0
chi_df_dt <- chi_df_dt[, -1] # remove lag 0

# boxplot all pixels values for chi temp
# Reshape data using gather function
df_gathered <- chi_df_dt %>% gather(key = "variable", value = "value")
df_gathered$group <- factor(as.integer(df_gathered$variable),
                            levels = seq(0, tmax))

# Plot boxplots
chitemp <- ggplot(df_gathered, aes(x = group, y = value)) +
  geom_boxplot(fill = btfgreen, alpha = 0.4, outlier.color = btfgreen) +
  btf_boxplot_theme +
  xlab(TeX(r"($\tau$ (hours))")) +
  ylab(TeX(r"($\widehat{\chi}(0,\tau)$)"))

chitemp

# save plot
filename <- paste(im_folder, "WLSE/comephore/full_temporal_chi_boxplot_",
                q_no0_temp,  ".png", sep = "")
ggsave(filename, width = 20, height = 15, units = "cm")

q_no0_temp <- 0.95 # quantile for temporal chi
tmax = 20
# Mean of chi
chimat_dt_mean <- temporal_chi(comephore, tmax, quantile = q_no0_temp,
                               mean = TRUE, zeros = T)
df_chi <- data.frame(lag = c(0:tmax), chi = chimat_dt_mean)

df_chi_not0 <- df_chi[df_chi$lag > 0, ]
wlse_temp <- get_estimate_variotemp(df_chi_not0, weights = "exp", summary = TRUE)
print(wlse_temp)
c2 <- as.numeric(wlse_temp[[1]])
beta2 <- as.numeric(wlse_temp[[2]])
alpha2 <- as.numeric(wlse_temp[[3]])
print(beta2)

dftemp <- data.frame(lag = log(df_chi_not0$lag), chi = eta(df_chi_not0$chi))

# remove first row
# dftemp <- dftemp[-1, ]

chitemp_eta_estim <- ggplot(dftemp, aes(x = lag, y = chi)) +
  geom_point(color = btfgreen) +
  btf_theme +
  xlab(TeX(r"($\log(\tau)$)")) +
  ylab(TeX(r"($\zeta(\widehat{\chi}(0,\tau))$)")) +
  geom_line(aes(x = lag, y = alpha2 * lag + c2),
            alpha = 0.5, color = "darkred", linewidth = 0.5)

chitemp_eta_estim

# save plot
filename <- paste(im_folder, "WLSE/comephore/full_temporal_chi_eta_estim_",
                q_no0_temp, ".png", sep = "")
ggsave(filename, width = 20, height = 15, units = "cm")

# Result WLSE
df_result <- data.frame(beta1 =  beta1,
                        beta2 = beta2,
                        alpha1 = alpha1,
                        alpha2 = alpha2)

colnames(df_result) <- c("beta1", "beta2", "alpha1", "alpha2")

kable(df_result, format = "latex") %>%
  kable_styling(bootstrap_options = c("striped", "hover", "condensed",
  "responsive"), latex_options = "H")


# CHOOSE EXTREME EPISODE FOR R-PARETO ##########################################

q <- 0.998 # quantile
# get central site from sites_coords
# comephore_subset <- comephore[rownames(comephore) >= "2020-01-01", ]
comephore_subset <- comephore
set_st_excess <- get_spatiotemp_excess(comephore_subset, q)

list_s <- set_st_excess$list_s
unique(unlist(list_s))
list_t <- set_st_excess$list_t

comephore[list_t[[1]], list_s[[1]]]

list_u <- set_st_excess$list_u
list_u[[1]]

# Spatio-temporal neighborhood parameters
min_spatial_dist <- 5 # in km
delta <- 30 # in hours
episode_size <- delta # size of the episode
s0t0_set <- get_s0t0_pairs(grid_coords_km, comephore_subset,
                            min_spatial_dist = min_spatial_dist,
                            episode_size = episode_size,
                            set_st_excess = set_st_excess,
                            n_max_episodes = 10000,
                            latlon = FALSE)

selected_points <- s0t0_set

# check that for all s0, t0 we have an excess above corresponding threshold
for (i in 1:length(selected_points$s0)) {
  s0 <- selected_points$s0[i]
  t0 <- selected_points$t0[i]
  u_s0 <- selected_points$u_s0[i]
  # check that the excess is above the threshold
  if (comephore_subset[t0, s0] <= u_s0) {
    stop(paste("Excess is not above threshold for s0 =", s0, "and t0 =", t0))
  }
}


# Threshold histogram
df_threshold <- data.frame(u_s0 = selected_points$u_s0)
breaks <- seq(floor(min(df_threshold$u_s0)), ceiling(max(df_threshold$u_s0)), by = 0.1)

ggplot(df_threshold, aes(x = u_s0)) +
  geom_histogram(breaks = breaks, fill = btfgreen, alpha = 0.5, color = "#5f5d5d") +
  btf_theme +
  xlab(TeX(paste0("Threshold for quantile $q = ", q, "$"))) +
  ylab("Count")
filename <- paste(im_folder, "optim/comephore/3km_threshold_histogram_q",
                  q * 1000, "_delta", delta, "_dmin", min_spatial_dist,
                  ".png", sep = "")
ggsave(filename, width = 20, height = 15, units = "cm")



n_episodes <- length(selected_points$s0)
print(n_episodes)
length(unique(selected_points$s0)) # can be same s0
length(unique(selected_points$t0)) # never same t0?
print(min(selected_points$u_s0)) # min threshold
t0_list <- selected_points$t0
s0_list <- selected_points$s0
u_list <- selected_points$u_s0
list_episodes_points <- get_extreme_episodes(selected_points, comephore_subset,
                                     episode_size = episode_size, unif = FALSE,
                                     beta = 0)

list_episodes <- list_episodes_points$episodes
episode <- list_episodes[[1]]
u <- u_list[[1]] 
# check the episode
head(episode)
library(ggplot2)
library(reshape2)  # for melting wide data to long format

# Convert matrix to data frame
index <- 10
sort(t0_list)
which(t0_list == t0_list[index])
episode_test <- list_episodes[[index]]
df_episode <- as.data.frame(episode_test)
df_episode$Time <- 0:(nrow(df_episode)-1)  # Add a time column
s0_list[index]
u_episode <- selected_points$u_s0[index]
# t0_episode <- t0_list[index]
# Convert from wide to long format
df_long <- melt(df_episode, id.vars = "Time")
head(df_long)
colnames(df_long) <- c("Time", "Pixel", "Value")
ggplot(df_long, aes(x = Time, y = Value, group = Pixel)) +
  geom_line(color = btfgreen) +
  geom_hline(yintercept = u_episode, color = "red", linetype = "dashed") +
  theme_minimal()

# filename <- paste(im_folder, "optim/comephore/extreme_episode", index, "_min", min_spatial_dist,
#                   "km_max", tmax, "h_delta_", delta, ".png", sep = "")
# # filename <- "test.png"
# ggsave(filename, width = 20, height = 15, units = "cm")


list_episodes_unif_points <- get_extreme_episodes(selected_points, comephore_subset,
                                      episode_size = episode_size, unif = TRUE)

list_episodes_unif <- list_episodes_unif_points$episodes

s0_list <- selected_points$s0
t0_list <- selected_points$t0
u_list <- selected_points$u_s0

library(parallel)

tau_vect <- 0:10
tmax <- max(tau_vect)
df_coords <- as.data.frame(grid_coords_m)
# Compute the lags and excesses for each conditional point
list_results <- mclapply(1:length(s0_list), function(i) {
  s0 <- s0_list[i]
  s0_coords <- df_coords[s0, ]
  # t0 <- t0_list[i]
  u <- u_list[i]
  episode <- list_episodes[[i]] / u
  ind_t0_ep <- 0 # index of t0 in the episode
  lags <- get_conditional_lag_vectors(df_coords, s0_coords, ind_t0_ep,
                                tau_vect, latlon = FALSE)
  # lags$hnorm <- lags$hnorm / 1000 # convert to km
  excesses <- empirical_excesses_rpar(episode, 1, lags, t0 = ind_t0_ep)
  lags$tau <- lags$tau * 3600 # convert tau from hours to seconds
  list(lags = lags, excesses = excesses)
}, mc.cores = detectCores() - 1)

list_lags <- lapply(list_results, `[[`, "lags")
list_excesses <- lapply(list_results, `[[`, "excesses")

# save workspace
save.image("workspace.RData")

# load("workspace.RData")
s0 <- s0_list[index]
s0_coords <- df_coords[s0, ]
excesses <- list_excesses[[10]]
excesses$kij
excesses[1:20, ]
df_lags <- list_lags[[1]]
tail(df_lags)

sum_excess_vect <- c()
for (i in 1:length(list_excesses)) {
  excesses <- list_excesses[[i]]
  sum_excess_vect <- c(sum_excess_vect, sum(excesses$kij))
}


# show distribution of sum of excesses
library(ggplot2)
df_excesses <- data.frame(sum_excess = sum_excess_vect)
ggplot(df_excesses, aes(x = sum_excess)) +
  geom_histogram(bins = 30, fill = btfgreen, alpha = 0.5, color = "#5f5d5d") +
  btf_theme +
  xlab("Sum of excesses") +
  ylab("Count") +
  ggtitle("Distribution of sum of excesses")

ggsave(filename = paste0(im_folder, "optim/comephore/sum_excesses_histogram_q_",
                         q * 1000, "_min_spatial_dist_", min_spatial_dist,
                         "km_tmax_", tmax, "h_delta_", delta, ".png"),
       width = 20, height = 15, units = "cm")


# ADD WIND DATA ################################################################

wind_per_episode <- Map(compute_wind_episode_vector_mean, list_episodes,
               MoreArgs = list(wind_df = wind_mtp))

wind_ep_df <- do.call(rbind, wind_per_episode)
head(wind_ep_df)
tail(wind_ep_df)


# OPTIMIZATION #################################################################

# Compute the wind vector for each episode (-FF because it's the wind direction)
# sin and cos are shifted because 0 degree means North
wind_ep_df$vx <- -wind_ep_df$FF_vector_mean * sin(wind_ep_df$DD_vector_mean * pi / 180)
wind_ep_df$vy <- -wind_ep_df$FF_vector_mean * cos(wind_ep_df$DD_vector_mean * pi / 180)

# if values of vx or vy are really close to 0, set them to 0
wind_ep_df$vx[abs(wind_ep_df$vx) < 1e-08] <- 0
wind_ep_df$vy[abs(wind_ep_df$vy) < 1e-08] <- 0

wind_df <- wind_ep_df[, c("vx", "vy")]
colnames(wind_df) <- c("vx", "vy")
# which episode has NA wind values
ind_NA <- which(is.na(wind_df$vx))

if (any(ind_NA > 0)) {
  # remove these episodes
  wind_opt <- wind_df[-ind_NA, ]
  episodes_opt <- list_episodes[-ind_NA]
  lags_opt <- list_lags[-ind_NA]
  excesses_opt <- list_excesses[-ind_NA]
} else {
  wind_opt <- wind_df
  episodes_opt <- list_episodes
  lags_opt <- list_lags
  excesses_opt <- list_excesses
}

init_param <- c(beta1, beta2, alpha1, alpha2, 1, 1)

result <- optim(par = init_param, fn = neg_ll_composite,
        list_lags = lags_opt, list_episodes = episodes_opt,
        list_excesses = excesses_opt, hmax = 7000,
        wind_df = wind_opt,
        latlon = FALSE,
        directional = TRUE,
        fixed_eta1 = FALSE,
        fixed_eta2 = TRUE,
        convert_in_hours = TRUE,
        convert_in_km = TRUE,
        method = "L-BFGS-B",
        lower = c(1e-08, 1e-08, 1e-08, 1e-08, 1e-08, 1e-08),
        upper = c(10, 10, 1.999, 1.999, 10, 10),
        control = list(maxit = 10000,
                      trace = 1,
                      parscale = c(1, 1, 1, 1, 1, 1)),
        hessian = F)

result

# [1] 0.1919899 1.2861977 0.2326258 0.6075999 0.6722028 0.9990000


# EMPIRICAL VARIOGRAM ##########################################################


library(gstat)
library(sp)

# Choisir un instant donné, ex. 2018-01-01 à 12h
t0 <- as.POSIXct("2018-01-01 12:00:00")

# Extraire les valeurs de pluie pour cet instant
row_t0 <- which(df_comephore$date == t0)
rain_values <- df_comephore[row_t0, -1]  # On enlève la colonne date

# Joindre avec coordonnées
rain_df <- data.frame(rain = as.numeric(rain_values), coords_pixels)

# Créer objet spatial
coordinates(rain_df) <- ~x+y
proj4string(rain_df) <- CRS("+init=epsg:2154")

# Calcul du variogramme empirique
vg_emp <- variogram(rain ~ 1, data = rain_df)
plot(vg_emp, main = paste("Variogramme spatial à", t0))


# VARIOGRAM PLOTS ##############################################################

# compute variogram with parameters
df_result <- data.frame(beta1 = result$par[1],
                        beta2 = result$par[2],
                        alpha1 = result$par[3],
                        alpha2 = result$par[4],
                        eta1 = result$par[5],
                        eta2 = result$par[6])

beta1 <- df_result$beta1
beta2 <- df_result$beta2
alpha1 <- df_result$alpha1
alpha2 <- df_result$alpha2
eta1 <- df_result$eta1
eta2 <- df_result$eta2


# Compute advection for each row
adv_df_ep <- wind_df %>%
  mutate(
    adv_x = (abs(vx)^eta1) * sign(vx) * eta2,
    adv_y = (abs(vy)^eta1) * sign(vy) * eta2
  )

# remove NA values
adv_df_ep <- adv_df_ep[!is.na(adv_df_ep$adv_x), ]
head(adv_df_ep)

episode <- list_episodes[[1]]

lags <- list_lags[[1]]

wind_ep_1 <- wind_df[1, ]

adv_ep_1 <- adv_df_ep[1, ]



# Vecteur d'advection pour l'épisode 1
vx <- wind_ep_1$vx
vy <- wind_ep_1$vy

adv_x <- (abs(vx)^eta1) * sign(vx) * eta2
adv_y <- (abs(vy)^eta1) * sign(vy) * eta2

# Recalage des lags
lags_recal <- lags %>%
  mutate(
    hx_recal = s1x - s2x + adv_x * tau,
    hy_recal = s1y - s2y + adv_y * tau,
    hnorm_recal = sqrt(hx_recal^2 + hy_recal^2),
    gamma_model = beta1 * (hnorm_recal^alpha1) + beta2 * (abs(tau)^alpha2)
  )
