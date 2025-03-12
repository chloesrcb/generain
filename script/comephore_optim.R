# remove objects from the environment
rm(list = ls())

# clear console
cat("\014")

# Load libraries and set theme
source("./script/load_libraries.R")

# LOAD DATA ####################################################################
comephore_raw <- read.csv("./data/comephore/zoom_3km.csv", sep = ",")
loc_px <- read.csv("./data/comephore/loc_px_zoom_3km.csv", sep = ",")

df_comephore <- comephore_raw

# Take only data after 2007
df_comephore <- df_comephore[df_comephore$date >= "2008-01-01", ]

# put date in index
rownames(df_comephore) <- df_comephore$date
comephore <- df_comephore[-1] # remove dates column
# Get distances matrix
dist_mat <- get_dist_mat(loc_px)
df_dist <- reshape_distances(dist_mat)
nsites <- nrow(loc_px)

# get wind data
wind_mtp <- read.csv("./data/wind/data_gouv/wind_mtp.csv")

# Convert datetime to POSIXct
wind_mtp$datetime <- as.POSIXct(wind_mtp$datetime,
                                format = "%Y-%m-%d %H:%M:%S", tz = "UTC")

# Apply function to the DD column
wind_mtp$cardDir <- sapply(wind_mtp$DD, convert_to_cardinal)
wind_mtp$cardDir <- as.character(wind_mtp$cardDir)  # Ensure it's character
wind_mtp$cardDir[is.na(wind_mtp$DD)] <- NA
# Check if NA values are properly handled
summary(wind_mtp)

head(wind_mtp$cardDir)

# DATA WITHOUT 0 ###############################################################

# Remove rows with all columns equal to 0
# comephore_no0 <- comephore[rowSums(comephore) != 0, ]

# nrow(comephore)
# nrow(comephore_no0)

# CHOOSE QUANTILE ##############################################################

# get a matrix of high quantiles for all pair
# q <- 0.99 # quantile
# list_count_quant <- quantile_matrix(q, comephore, qlim = TRUE, zeros = TRUE,
#                                     count_min = 150) # with removing zeros
# quant_mat <- list_count_quant[1][[1]]
# count_mat <- list_count_quant[2]
library(ggplot2)
par(mfrow = c(1, 2))
comephore_pair <- comephore[, c(12, 2)]
comephore_pair_no0 <- comephore_pair[rowSums(comephore_pair) > 0, ]
chiplot(comephore_pair_no0, xlim = c(0.85, 1), ylim1 = c(0.5, 1), which = 1,
        qlim = c(0.85, 0.999))
abline(v = 0.96, col = "red", lty = 2)

# Same with all zeros
comephore_pair <- comephore[, c(12, 2)]
chiplot(comephore_pair, xlim = c(0.99, 1), ylim1 = c(0.5, 1), which = 1,
        qlim = c(0.99, 0.9999))
abline(v = 0.998, col = "red", lty = 2)

# count conjoint excesses
q <- 0.96
# uniformize the data
n <- nrow(comephore_pair_no0)
data_unif <- cbind(rank(comephore_pair_no0[, 1]) / (n + 1),
                          rank(comephore_pair_no0[, 2]) / (n + 1))

count_excesses <- sum(data_unif[, 1] > q & data_unif[, 2] > q)
print(count_excesses)

# count conjoint excesses
q <- 0.9975
# uniformize the data
n <- nrow(comephore_pair)
data_unif <- cbind(rank(comephore_pair[, 1]) / (n + 1),
                          rank(comephore_pair[, 2]) / (n + 1))

count_excesses <- sum(data_unif[, 1] > q & data_unif[, 2] > q)
print(count_excesses)

q_no0_spa <- 0.96
# get the corresponding quantile in the data without 0 and with 0
q_value_no0 <- quantile(c(comephore_pair_no0[,1], comephore_pair_no0[,2]), 
                        probs = q_no0_spa, na.rm = TRUE)

q_equiv_prob <- ecdf(c(comephore_pair[,1], comephore_pair[,2]))(q_value_no0)
print(q_equiv_prob)

q_equiv <- quantile(c(comephore[,1], comephore[,15]), probs = q_equiv_prob, 
              na.rm = TRUE)
print(q_equiv)  # quantile and threshold

q_spa <- round(q_equiv_prob, 4)

# Temporal chi
par(mfrow = c(1,2))
rain_nolag <- comephore$p102[1:(length(comephore$p102) - 5)]
rain_lag <- comephore$p102[6:length(comephore$p102)]
comephore_pair <- cbind(rain_nolag, rain_lag)
comephore_pair_no0 <- comephore_pair[rowSums(comephore_pair) > 0, ]
chiplot(comephore_pair_no0, xlim = c(0.8, 0.99), ylim1 = c(0, 1), which = 1,
        qlim = c(0.8, 0.99))
abline(v = 0.96, col = "red", lty = 2)

comephore_pair <- cbind(rain_nolag, rain_lag)
# comephore_pair <- comephore_pair[rowSums(comephore_pair) > 0, ]
chiplot(comephore_pair, xlim = c(0.98, 1), ylim1 = c(0, 1), which = 1,
        qlim = c(0.98, 0.999))
abline(v = 0.995, col = "red", lty = 2)

# count conjoint excesses
q_no0_spa <- 0.96

q_no0_temp <- 0.94
# get the corresponding quantile in the data without 0 and with 0
q_value_no0 <- quantile(c(comephore_pair_no0[,1], comephore_pair_no0[,2]),
                        probs = q_no0_temp, na.rm = TRUE)

q_equiv_prob <- ecdf(c(comephore_pair[,1], comephore_pair[,2]))(q_value_no0)
print(q_equiv_prob)

q_temp <- round(q_equiv_prob, 4)
q_equiv <- quantile(c(comephore[,1], comephore[,15]), probs = q_equiv_prob, 
              na.rm = TRUE)
print(q_equiv)  # quantile and threshold

# BUHL WLSE ####################################################################

# Temporal chi
tmax <- 10
# Temporal chi with spatial lag fixed at 0
# compute chiplot of every site with itself but lagged in time
chimat_dtlag <- temporal_chi(comephore, quantile = q_no0_temp, tmax = tmax,
                             mean = FALSE, zeros = FALSE)

chi_df_dt <- data.frame(chimat_dtlag)
colnames(chi_df_dt) <- c(1:tmax) # temporal lags from 1 to tmax
rownames(chi_df_dt) <- c(1:nsites) # stations

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
filename <- paste(im_folder, "WLSE/comephore/from2008_temporal_chi_boxplot_", q_no0_temp,
                    ".png", sep = "")
ggsave(filename, width = 20, height = 15, units = "cm")

# Mean of chi
chimat_dt_mean <- temporal_chi(comephore, tmax, quantile = 0.93,
                               mean = TRUE, zeros = FALSE)
# get h axis in minutes ie x5 minutes
df_chi <- data.frame(lag = c(1:tmax), chi = chimat_dt_mean)
chitemp_plot <- ggplot(df_chi, aes(x = lag, y = chi)) +
  geom_point(color = btfgreen) +
  btf_theme +
  xlab(TeX(r"($\tau$ (hours))")) +
  ylab(TeX(r"($\widehat{\chi}(0,\tau)$)"))

chitemp_plot

# remove the first value
# df_chi <- df_chi[-1, ]
wlse_temp <- get_estimate_variotemp(df_chi, weights = "exp", summary = TRUE)
print(wlse_temp)
c2 <- wlse_temp[[1]]
beta2 <- wlse_temp[[2]]
alpha2 <- wlse_temp[[3]]

dftemp <- data.frame(lag = log(df_chi$lag), chi = eta(df_chi$chi))

# remove first row
# dftemp <- dftemp[-1, ]

chitemp_eta_estim <- ggplot(dftemp, aes(x = lag, y = chi)) +
  geom_point(color = btfgreen) +
  btf_theme +
  xlab(TeX(r"($\log(\tau)$)")) +
  ylab(TeX(r"($\eta(\widehat{\chi}(0,\tau))$)")) +
  geom_line(aes(x = lag, y = alpha2 * lag + c2),
            alpha = 0.5, color = "darkred", linewidth = 0.5)

chitemp_eta_estim

# save plot
filename <- paste(im_folder, "WLSE/comephore/from2008_temporal_chi_eta_estim_",
                q_no0_temp, ".png", sep = "")
ggsave(filename, width = 20, height = 15, units = "cm")

# Spatial chi
df_dist <- reshape_distances(dist_mat)
df_dist$value <- ceiling(df_dist$value / 100) * 100  # / 1000 in km
df_dist_km <- df_dist
df_dist_km$value <- df_dist$value / 1000

h_vect <- sort(unique(df_dist_km$value))
h_vect <- h_vect[h_vect > 0]  # remove 0
hmax <- h_vect[10] # 10th value

hmax <- 7

q <- q_no0_spa
chispa_df <- spatial_chi_alldist(df_dist_km, data_rain = comephore,
                                 quantile = q, hmax = hmax, zeros = FALSE)

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
filename <- paste(im_folder, "WLSE/comephore/from2008_spatial_chi_", q,
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
  ylab(TeX(r"($\eta(\widehat{\chi}(h, 0))$)")) +
  geom_line(aes(x = lagspa, y = alpha1 * lagspa + c1), alpha = 0.6,
            color = "darkred", linewidth = 0.5)

chispa_eta_estim

# save plot
filename <- paste(im_folder, "WLSE/comephore/from2008_spatial_chi_eta_estim_", q,
                    ".png", sep = "")
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

# Find zero gaps
find_zero_gaps <- function(data) {
  data <- as.matrix(data)  # Ensure matrix format
  
  zero_gaps_list <- lapply(seq_len(ncol(data)), function(col_idx) {
    col_data <- data[, col_idx] == 0  # Identify zero positions
    zero_indices <- which(col_data)  # Get row indices where values are zero
    
    if (length(zero_indices) < 4) return(NULL)  # Ignore short series
    
    # Find consecutive zero sequences
    gap_starts <- zero_indices[c(TRUE, diff(zero_indices) > 1)]  # Start of each gap
    gap_lengths <- rle(diff(c(zero_indices, Inf)))$lengths  # Lengths of gaps
    
    # Keep only gaps of at least 4 zeros
    valid_gaps <- gap_lengths >= 4
    
    if (!any(valid_gaps)) return(NULL)  # No valid gaps

    # Return a dataframe for this column
    data.frame(
      column = colnames(data)[col_idx],
      start = gap_starts[valid_gaps],
      gap_length = gap_lengths[valid_gaps]
    )
  })
  
  zero_gaps_df <- do.call(rbind, zero_gaps_list)
  
  return(zero_gaps_df)
}

zero_gaps_df <- find_zero_gaps(comephore)
print(zero_gaps_df)
max(zero_gaps_df$gap_length)


find_extreme_episodes <- function(data, quantile_value, min_gap = 4) {
  data <- as.matrix(data)  # Conversion en matrice pour traitement plus rapide
  
  # Convertir en valeurs uniformes, colonne par colonne
  data_unif <- apply(data, 2, function(col) rank(col, na.last = "keep") / (length(col) + 1))
  
  # Identifier les indices des événements extrêmes (valeurs supérieures au quantile)
  extreme_indices <- which(data_unif > quantile_value, arr.ind = TRUE)
  
  if (nrow(extreme_indices) == 0) {
    return(NULL)  # Aucun événement extrême trouvé
  }
  
  # Initialisation de la liste pour les épisodes extrêmes
  extreme_episodes <- list()

  for (i in 1:nrow(extreme_indices)) {
    row_i <- extreme_indices[i, 1]
    col_i <- extreme_indices[i, 2]

    # Trouver la séquence de zéros avant et après l'extrême
    # before_extreme <- which(data[1:(row_i - 1), col_i] == 0)  # Zéros avant
    # after_extreme <- which(data[(row_i + 1):nrow(data), col_i] == 0)  # Zéros après

    # Get values almost extreme
    before_extreme <- which(data_unif[1:(row_i - 1), col_i] <= quantile_value)
    after_extreme <- which(data[(row_i + 1):nrow(data), col_i] <= quantile_value)

    # Enregistrer l'épisode uniquement s'il y a un gap de zéro d'au moins 4
    if (length(before_extreme) > 0 & length(after_extreme) > 0) {
      first_after <- min(after_extreme)
      last_before <- max(before_extreme)
      start_extreme <- last_before
      end_extreme <- row_i + first_after
      extreme_episodes[[length(extreme_episodes) + 1]] <- data.frame(
        column = colnames(data)[col_i],
        start = start_extreme,
        end = end_extreme,
        length = end_extreme - start_extreme + 1
      )
    }
  }

  # Si des épisodes extrêmes ont été trouvés, les combiner en un seul dataframe
  if (length(extreme_episodes) > 0) {
    extreme_episodes_df <- do.call(rbind, extreme_episodes)
    return(extreme_episodes_df)
  } else {
    return(NULL)
  }
}


extreme_episodes_df <- find_extreme_episodes(comephore, 0.998)
# p254  18376  18394     19
# comephore[18376:18410, "p254"]
max(extreme_episodes_df$length)
sort(unique(extreme_episodes_df$length))
par(mfrow = c(1, 1))
hist(extreme_episodes_df$length, breaks = 20, col = btfgreen, border = "black")
which.max(extreme_episodes_df$length)
extreme_episodes_df[which.max(extreme_episodes_df$length), ]
comephore[87840:87870, "p250"]



# Spatio-temporal neighborhood parameters
min_spatial_dist <- 5  # in km
delta <- 12 # step for the episode before and after the max value

# Get coords
sites_coords <- loc_px[, c("Longitude", "Latitude")]
rownames(sites_coords) <- loc_px$pixel_name

q <- 0.997 # quantile

# Get the selected episodes
selected_points <- select_extreme_episodes(sites_coords, comephore, q,
                                        min_spatial_dist, delta = delta,
                                        n_max_episodes = 10000,
                                        time_ext = 0)

n_episodes <- length(selected_points$s0)
print(n_episodes)
length(unique(selected_points$s0)) # can be same s0
length(unique(selected_points$t0)) # never same t0?
print(min(selected_points$u_s0)) # min threshold

list_episodes_points <- get_extreme_episodes(selected_points, comephore,
                                      delta = delta, unif = FALSE)

list_episodes <- list_episodes_points$episodes

length(list_episodes)
# Verif first episode
episode <- list_episodes[[1]]
# find col and row of the max value
max_value <- max(episode)
max_value_coords <- which(episode == max_value, arr.ind = TRUE)
s0_col <- colnames(episode)[max_value_coords[2]]
t0_row <- rownames(episode)[max_value_coords[1]]

# selected_points$s0[1]
# selected_points$t0[1]
# rownames(comephore)[selected_points$t0[1]]

# histogram of the number of episodes per s0
df_hist <- data.frame(table(selected_points$s0))
colnames(df_hist) <- c("s0", "count")
ggplot(df_hist, aes(x = count)) +
  geom_histogram(binwidth = 1, fill = btfgreen, color = "black") +
  btf_theme +
  xlab("Number of episodes") +
  ylab("Count")

s0_max_ep <- df_hist$s0[which.max(df_hist$count)]
beta <- 0
list_overlaps <- c()
for (s in unique(selected_points$s0)) {
  overlaps <- check_intervals_overlap(s, selected_points, delta, beta)
  list_overlaps <- c(list_overlaps, overlaps)
}
any(overlaps)



# get month and year for p236
# list_date <- lapply(1:length(list_p236), function(ep) {
#   list_date_t0 <- rownames(list_p236[[ep]])[1]
#   date_t0 <- as.Date(list_date_t0)
#   day <- format(date_t0, "%d")
#   month <- format(date_t0, "%m")
#   year <- format(date_t0, "%Y")
#   list(day = day, month = month, year = year)
# })

# list_date <- do.call(rbind, list_date)
# rownames(list_date) <- c(1:47)
# colnames(list_date) <- c("day", "month", "year")
# year = "2010"
# list_date[list_date[,3] == year, ]
# # 34, 47
list_episodes_unif_points <- get_extreme_episodes(selected_points, comephore,
                                      delta = delta, unif = TRUE)

list_episodes_unif <- list_episodes_unif_points$episodes
# list_episodes[[100]] # check the episode
# nrow(list_episodes[[1]]) # size of episode

# boxplot(list_episodes[[1]])
s0_list <- selected_points$s0
t0_list <- selected_points$t0
u_list <- selected_points$u_s0

tmax <- 10
# Compute the lags and excesses for each conditional point
list_results <- mclapply(1:length(s0_list), function(i) {
  s0 <- s0_list[i]
  s0_coords <- sites_coords[s0, ]
  # t0 <- t0_list[i]
  episode <- list_episodes_unif[[i]]
  u <- u_list[i]
  ind_t0_ep <- delta + 1  # index of t0 in the episode
  lags <- get_conditional_lag_vectors(sites_coords, s0_coords, ind_t0_ep,
                                  tau_max = tmax, latlon = TRUE)
  lags$hx <- lags$hx / 1000  # in km
  lags$hy <- lags$hy / 1000  # in km
  lags$hnorm <- lags$hnorm / 1000  # in km
  excesses <- empirical_excesses(episode, q, lags, type = "rpareto",
                                  t0 = ind_t0_ep)
  list(lags = lags, excesses = excesses)
}, mc.cores = detectCores() - 1)

list_lags <- lapply(list_results, `[[`, "lags")
list_excesses <- lapply(list_results, `[[`, "excesses")

# list_lags[[1]]
# list_excesses[[1]]$kij

s0 <- s0_list[1]
s0_coords <- sites_coords[s0, ]
excesses <- list_excesses[[1]]
# sort excesses by hnorm from smallest to largest
excesses <- excesses[order(excesses$hnorm), ]
excesses$kij
excesses[1:20, ]
# # Initialisation des listes pour les résultats
# all_k_sums <- c()
# all_binomials <- c()  # Liste pour stocker les résultats de la binomiale théorique
# init_param <- c(beta1, beta2, alpha1, alpha2, 1, 1)
# for (i in 1:length(list_excesses)) {
#   excesses <- list_excesses[[i]]
#   df_lags <- list_lags[[i]]
#   wind_vect <- c(wind_ep_df$vx[i], wind_ep_df$vy[i])
#   chi_theorical <- theoretical_chi(init_param, df_lags, latlon = TRUE,
#                                    wind_vect = wind_vect)
#   tau_vect <- unique(excesses$tau)

#   # Pour chaque valeur de hnorm
#   for (hnorm in unique(df_lags$hnorm)) {

#     all_kij_h_tau <- c()  # Liste pour stocker kij pour chaque hnorm et tau
#     for (tau in tau_vect) {
#       # Filtrer les données pertinentes
#       kij <- excesses$kij[excesses$tau == tau & excesses$hnorm == hnorm]
#       nij <- excesses$Tobs[excesses$tau == tau & excesses$hnorm == hnorm]

#       # Ajouter les kij dans la liste
#       all_kij_h_tau <- c(all_kij_h_tau, kij)
#     }

#     sum_k_h <- sum(all_kij_h_tau)  # Somme des kij pour hnorm et tau
#     all_k_sums <- c(all_k_sums, sum_k_h)

#     # Calcul de la probabilité théorique de la binomiale
#     chi_h_tau <- chi_theorical$chi[chi_theorical$hnorm == hnorm]

#     # Calcul de la binomiale théorique pour les essais et probabilité
#     trials <- length(all_kij_h_tau)  # Nombre d'essais
#     binomial_theoretical <- dbinom(sum_k_h, trials, chi_h_tau) # binomial
#     all_binomials <- c(all_binomials, binomial_theoretical)
#   }
# }


# ADD WIND DATA ################################################################


list_episodes_points <- get_extreme_episodes(selected_points, comephore,
                                      delta = delta, unif = FALSE)

list_episodes <- list_episodes_points$episodes

wind_per_episode <- Map(compute_wind_episode, list_episodes, s0_list, u_list,
               MoreArgs = list(wind_df = wind_mtp, delta = delta))

wind_ep_df <- do.call(rbind, wind_per_episode)
head(wind_ep_df)

sum(is.na(wind_ep_df$DD_deg))
sum(is.na(wind_ep_df$DD_t0))
# sum(is.na(wind_ep_df$DD_excess))
sum(is.na(wind_ep_df$FF))

# wind_ep_df$cardDirt0 <- sapply(wind_ep_df$DD_t0, convert_to_cardinal)

wind_ep_df$cardDir <- factor(wind_ep_df$cardDir,
                       levels = c("N", "NE", "E", "SE", "S", "SW", "W", "NW"))


# Define wind speed categories and reorder levels to change legend order
wind_ep_df <- wind_ep_df %>%
  mutate(FF_interval = factor(cut(FF,
                            breaks = c(0, 2, 5, 10, Inf),
                            labels = c("<2", "2-5", "5-10", ">10"),
                            right = FALSE),
                    levels = c(">10", "5-10", "2-5", "<2")))  # Reverse order

# Define wind direction order
wind_ep_df$cardDir <- factor(wind_ep_df$cardDir,
                    levels = c("N", "NE", "E", "SE", "S", "SW", "W", "NW"))

# Custom color palette
custom_colors <- c("#0335258e", "#1a755a8e", "#5b99868e", "#98d6c48e")

# Plot wind rose
ggplot(wind_ep_df, aes(x = cardDir, fill = FF_interval)) +
  geom_bar(position = "stack", width = 1, color = "black") +
  coord_polar(start = 15 * pi / 8) +
  scale_fill_manual(values = custom_colors, name = "Force (m/s)") +
  labs(x = "Direction", y = "Count", title = "") +
  theme_minimal() +
  theme(axis.text.x = element_text(size = 12, face = "bold"),
        legend.title = element_text(size = 12, face = "bold"),
        legend.text = element_text(size = 10))

# Save plot
end_filename <- paste(n_episodes, "_ep_", min_spatial_dist, "_km_",
                      tmax, "_h_delta_", delta, ".png", sep = "")
filename <- paste(im_folder, "wind/datagouv/wind_card_dir_", end_filename,
                    sep = "")
ggsave(filename, width = 20, height = 15, units = "cm")


# Ensure cardDir is a factor with correct order
wind_ep_df$cardDir_t0 <- factor(wind_ep_df$cardDir_t0,
                        levels = c("N", "NE", "E", "SE", "S", "SW", "W", "NW"))

# Plot wind rose
ggplot(wind_ep_df, aes(x = cardDir_t0, fill = FF_interval)) +
  geom_bar(position = "stack", width = 1, color = "black") +
  coord_polar(start = 15 * pi / 8) +
  scale_fill_manual(values = custom_colors, name = "Force (m/s)") +
  labs(x = "Direction in t0", y = "Count", title = "") +
  theme_minimal() +
  theme(axis.text.x = element_text(size = 12, face = "bold"),
        legend.title = element_text(size = 12, face = "bold"),
        legend.text = element_text(size = 10))

# save the wind data
filename <- paste(im_folder, "wind/datagouv/wind_card_dir_t0_", end_filename,
                    sep = "")
ggsave(filename, width = 20, height = 15, units = "cm")

# Plot wind rose for DD_t0
ggplot(wind_ep_df, aes(x = DD_t0)) +
  geom_bar(position = "stack", width = 3, color = "#5048489d", fill = btfgreen,
        alpha = 0.8) +
  coord_polar(start = 15.85 * pi / 8) +
  labs(x = "Direction in t0 (degree)", y = "Count", title = "") +
  theme_minimal() +
  theme(axis.text.x = element_text(size = 12, face = "bold"),
        legend.title = element_text(size = 12, face = "bold"),
        legend.text = element_text(size = 10))

# save the wind data
filename <- paste(im_folder, "wind/datagouv/wind_dir_t0_", end_filename,
                    sep = "")
ggsave(filename, width = 20, height = 15, units = "cm")

# OPTIMIZATION #################################################################

# TODO: TO PUT IN TESTS
# df_lags <- list_lags[[2]]
# excesses <- list_excesses[[2]]
# wind_vect <- c(wind_km_h$vx[2], wind_km_h$vy[2])
# eta1 <- 1
# eta2 <- 1
# adv <- (abs(wind_vect)^eta1) * sign(wind_vect) * eta2
# print(adv)
# init_param <- c(beta1, beta2, alpha1, alpha2, eta1, eta2)
# neg_ll(init_param, df_lags, excesses, wind_vect = wind_vect,
#                       hmax = 7)
# init_adv <- c(beta1, beta2, alpha1, alpha2, 0, 0)
# neg_ll(init_adv, df_lags, excesses, wind_vect = NA,
#                       hmax = 7)
# chi_adv <- theorical_chi(init_adv, df_lags, NA)
# chi_adv[6, ]

# init_param <- c(beta1, beta2, alpha1, alpha2, 1, 1.5)
# neg_ll_composite_rpar(init_param, list_lags, list_excesses, wind_km_h, hmax = 7)
# 51130.74 for eta1 = 1, eta2 = 1
# 51130.67 for eta1 = 0.9, eta2 = 1.5

# convert FF into km/h
# wind_km_h <- wind_ep_df
# wind_km_h$FF <- wind_ep_df$FF * 3.6
# wind_km_h$vx <- wind_km_h$FF * cos(wind_km_h$DD_t0)
# wind_km_h$vy <- wind_km_h$FF * sin(wind_km_h$DD_t0)
# head(wind_km_h)

# Compute the wind vector for each episode (-FF because it's the wind direction)
wind_ep_df$vx <- -wind_ep_df$FF * sin(wind_ep_df$DD_t0 * pi / 180)
wind_ep_df$vy <- -wind_ep_df$FF * cos(wind_ep_df$DD_t0 * pi / 180)
# vx = -FF * sin(DD * pi / 180)  # Conversion de DD en radians
# vy = -FF * cos(DD * pi / 180)  # Conversion de DD en radians

head(wind_ep_df)

wind_df <- wind_ep_df[, c("DD_t0", "vx", "vy")]
colnames(wind_df) <- c("dir", "vx", "vy")
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

# init_param <- c(beta1, 1.15, alpha1, 0.76, 1, 1)
eta1 <- 1
a_ratio <- 1
phi <- 0.5
beta1 = 0.01
beta2 = 0.8
alpha1 = 1.5
alpha2 = 0.8
init_param <- c(beta1, beta2, alpha1, alpha2, 1, 1)

# q <- 1
result <- optim(par = init_param, fn = neg_ll_composite,
        list_lags = lags_opt,
        list_excesses = excesses_opt, hmax = 7,
        wind_df = wind_opt,
        latlon = TRUE,
        directional = TRUE,
        method = "L-BFGS-B",
        lower = c(1e-08, 1e-08, 1e-08, 1e-08, 1e-08, 1e-08),
        upper = c(Inf, Inf, 1.999, 1.999, Inf, Inf),
        control = list(maxit = 10000,
                      trace = 1),
        hessian = F)
result



# Check the convergence
if (result$convergence != 0) {
  warning("The optimization did not converge.")
}

# Extract the results
df_result <- data.frame(quantile = q,
                        delta = delta,
                        minDist = min_spatial_dist,
                        beta1 =  round(result$par[1], 3),
                        beta2 = round(result$par[2], 3),
                        alpha1 = round(result$par[3], 3),
                        alpha2 = round(result$par[4], 3),
                        eta1 = round(result$par[5], 3),
                        eta2 = round(result$par[6], 3))

kable(df_result, format = "latex") %>%
  kable_styling(bootstrap_options = c("striped", "hover", "condensed",
  "responsive"), latex_options = "H")


# fix parameters
neg_ll_composite_fixpar <- function(beta1, beta2, alpha1, alpha2, eta1, a_ratio, phi, list_lags,
                    list_excesses, wind_df, hmax = NA, latlon = TRUE) {
  params <- c(beta1, beta2, alpha1, alpha2, eta1, a_ratio, phi)
  nll_composite <- neg_ll_composite(params, list_lags,
                                      list_excesses, wind_df,
                                      hmax = hmax, latlon = latlon,
                                      directional = TRUE)
  return(nll_composite)
}

library(bbmle)
result <- mle2(neg_ll_composite_fixpar,
              start = list(beta1 = init_param[1],
                           beta2 = init_param[2],
                           alpha1 = init_param[3],
                           alpha2 = init_param[4],
                           eta1 = 1,
                           a_ratio = 1,
                           phi=0.5),
              data = list(list_lags = lags_opt,
                  list_excesses = excesses_opt, hmax = 7,
                  wind_df = wind_opt,
                  method = "L-BFGS-B",
                  lower = c(1e-08, 1e-08, 1e-08, 1e-08, 1e-08, 1e-08),
                  upper = c(Inf, Inf, 1.999, 1.999, Inf, Inf)),
              control = list(maxit = 10000),
              fixed = list(beta2 = init_param[2],
                           alpha2 = init_param[4]))
result


# Vector of test values for eta1
eta1_values <- seq(1e-08, 2, length.out = 50)

# Store for negative log-likelihood values
neg_log_likelihoods <- numeric(length(eta1_values))

# Loop over different values of eta1
for (i in seq_along(eta1_values)) {
  params <- c(beta1, beta2, alpha1, 1, eta1_values[i], 1)
  neg_ll <- neg_ll_composite(params, list_lags = lags_opt,
                              list_excesses = excesses_opt, hmax = 7,
                              wind_df = wind_opt)
  neg_log_likelihoods[i] <- neg_ll
}

par(mfrow = c(1, 1))
# Plotting the likelihood function
plot(eta1_values, neg_log_likelihoods, type = "l", col = "blue", lwd = 2,
     xlab = expression(eta[1]), ylab = "- Log-likelihood",
     main = "Likelihood profile as a function of eta1")


eta1_values <- seq(0.001, 1.5, length.out = 30)
beta2_values <- seq(0.1, 3, length.out = 30)

likelihood_matrix <- matrix(NA, nrow = length(eta1_values), ncol = length(beta2_values))

for (i in seq_along(eta1_values)) {
  for (j in seq_along(beta2_values)) {
    params <- c(beta1, beta2_values[j], alpha1, alpha2, eta1_values[i], 1)
    likelihood_matrix[i, j] <- neg_ll_composite(params, list_lags = lags_opt,
                                                list_excesses = excesses_opt, hmax = 7,
                                                wind_df = wind_opt)
  }
}

contour(eta1_values, beta2_values, likelihood_matrix, xlab = expression(eta[1]),
        ylab = expression(beta[2]), main = "Log-likelihood contours for eta1 and beta2")


# GENERATE TABLE FROM TEST #####################################################

library(knitr)
library(kableExtra)

# Configurations
q_values <- c(0.997, 0.998, 0.999)
delta_values <- c(10, 12, 15)
min_spatial_dist_values <- c(3, 5)
tau_ranges <- list(-10:10)
wind_direction_type <- c("t_0", "freq")

# Init results table
results_table <- list()
n_ep_list <- c()

for (q in q_values) {
  for (delta in delta_values) {
    for (min_spatial_dist in min_spatial_dist_values) {
      for (tau_range in tau_ranges) {
        for (wind_dir in wind_direction_type) {
          # print configuration
          print(paste("q:", q, "delta:", delta, "minDist:", min_spatial_dist,
                      "taumax:", max(tau_range), "direction:", wind_dir))

          print("Select (s0, t0) pairs for extreme episodes")
          #  Select (s0, t0) pairs for extreme episodes
          selected_points <- select_extreme_episodes(sites_coords, comephore, q,
                                              min_spatial_dist, delta = delta,
                                              n_max_episodes = 10000,
                                              time_ext = 0)
          n_episodes <- length(selected_points$s0)
          n_ep_list <- c(n_ep_list, n_episodes)
          # Get the extreme episodes
          list_episodes_unif_points <- get_extreme_episodes(selected_points, comephore,
                                          delta = delta, unif = TRUE)
          list_episodes_unif <- list_episodes_unif_points$episodes
          s0_list <- selected_points$s0 # List of s0
          t0_list <- selected_points$t0 # List of t0
          u_list <- selected_points$u_s0 # List of threshold
          tmax <- max(tau_range)
          print("Compute the lags and excesses for each conditional point")
          # Get the lags and excesses for each conditional point
          list_results <- mclapply(1:length(s0_list), function(i) {
            s0 <- s0_list[i]
            s0_coords <- sites_coords[s0, ]
            episode <- list_episodes_unif[[i]]
            # u <- u_list[i]
            ind_t0_ep <- delta  # index of t0 in the episode
            lags <- get_conditional_lag_vectors(sites_coords, s0_coords,
                                      ind_t0_ep, tau_max = tmax, latlon = TRUE)
            lags$hx <- lags$hx / 1000  # in km
            lags$hy <- lags$hy / 1000  # in km
            lags$hnorm <- lags$hnorm / 1000  # in km
            excesses <- empirical_excesses(episode, q, lags, type = "rpareto",
                                           t0 = ind_t0_ep)
            list(lags = lags, excesses = excesses)
          }, mc.cores = detectCores() - 1)

          list_lags <- lapply(list_results, `[[`, "lags")
          list_excesses <- lapply(list_results, `[[`, "excesses")

          list_episodes_points <- get_extreme_episodes(selected_points, comephore,
                                          delta = delta, unif = FALSE)
          list_episodes <- list_episodes_points$episodes
          # Get wind components for each episode
          wind_per_episode <- Map(compute_wind_episode, list_episodes, s0_list,
                            u_list,
                            MoreArgs = list(wind_df = wind_mtp, delta = delta))

          wind_ep_df <- do.call(rbind, wind_per_episode)

          if (wind_dir == "t_0") { # Use wind direction at t0
            wind_ep_df$vx <- -wind_ep_df$FF * cos(wind_ep_df$DD_t0 * pi / 180)
            wind_ep_df$vy <- -wind_ep_df$FF * sin(wind_ep_df$DD_t0 * pi / 180)
          } else { # Use most frequent wind direction
            wind_ep_df$vx <- -wind_ep_df$FF * cos(wind_ep_df$DD * pi / 180)
            wind_ep_df$vy <- -wind_ep_df$FF * sin(wind_ep_df$DD * pi / 180)
          }

          # Get wind vectors in a dataframe
          wind_df <- wind_ep_df[, c("vx", "vy")]

          # Verify if there are NA values
          ind_NA <- which(is.na(wind_df$vx))

          if (any(ind_NA > 0)) { # Remove episodes with NA values
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

          # Initial parameters
          init_param <- c(beta1, beta2, alpha1, alpha2, 1, 1)

          print("Optimization")
          # Optimisation
          result <- optim(par = init_param, fn = neg_ll_composite,
                          list_lags = lags_opt,
                          list_excesses = excesses_opt, hmax = 7,
                          wind_df = wind_opt,
                          latlon = TRUE,
                          directional = TRUE,
                          method = "L-BFGS-B",
                          lower = c(1e-08, 1e-08, 1e-08, 1e-08, 1e-08, 1e-08),
                          upper = c(Inf, Inf, 1.999, 1.999, Inf, Inf),
                          control = list(maxit = 10000,
                                         trace = 1),
                          hessian = F)

          # Check the convergence
          if (result$convergence != 0) { # No convergence
            df_result <- data.frame(
              quantile = q,
              delta = delta,
              minDist = min_spatial_dist,
              tau = tmax,
              direction = wind_dir,
              beta1 = "NC",
              beta2 = "NC",
              alpha1 = "NC",
              alpha2 = "NC",
              eta1 = "NC",
              eta2 = "NC"
            )
          } else { # Convergence
            # Extract the results
            df_result <- data.frame(
              quantile = q,
              delta = delta,
              minDist = min_spatial_dist,
              tau = tmax,
              direction = wind_dir,
              beta1 = round(result$par[1], 3),
              beta2 = round(result$par[2], 3),
              alpha1 = round(result$par[3], 3),
              alpha2 = round(result$par[4], 3),
              eta1 = round(result$par[5], 3),
              eta2 = round(result$par[6], 3)
            )
          }

          print(df_result)
          # Add the results to the table
          results_table[[length(results_table) + 1]] <- df_result
        }
      }
    }
  }
}

# Combine the results into a single dataframe
final_results <- do.call(rbind, results_table)

# Latex table
kable(final_results, format = "latex",
      caption = "Optimization results with wind components as covariates") %>%
  kable_styling(bootstrap_options = c("striped", "hover", "condensed",
                                            "responsive"), latex_options = "H")


# ADVECTION ESTIMATION #########################################################

# get first result
df_result <- final_results[31, ]
df_result

# Estimated advection vector with wind components
adv_hat <- abs(wind_df) ^ as.numeric(df_result$eta1) * sign(wind_df) *
                            as.numeric(df_result$eta2)
print(adv_hat)

# remove NA values
adv_hat <- adv_hat[!is.na(adv_hat$vx), ]

# get wind FF and DD from adv_hat
adv_hat$FF_hat_kmh <- sqrt(adv_hat$vx^2 + adv_hat$vy^2) # Wind magnitude in km/h
adv_hat$FF_hat_ms <- adv_hat$FF_hat_kmh / 3.6  # Convert to m/s
adv_hat$FF_hat_mms <- adv_hat$FF_hat_ms * 1000  # Convert to mm/s
adv_hat$DD_hat <- atan2(adv_hat$vy, adv_hat$vx) * 180 / pi
head(adv_hat)

adv_hat$DD_hat <- ifelse(adv_hat$DD_hat < 0, adv_hat$DD_hat + 360,
                        adv_hat$DD_hat)
head(adv_hat)

adv_hat$cardDir <- sapply(adv_hat$DD_hat, convert_to_cardinal)
head(adv_hat)

# Define wind speed categories and reorder levels to change legend order
wind_plot_df <- adv_hat %>%
  mutate(FF_interval = factor(cut(FF_hat_mms,
                            breaks = c(0, 100, 120, 150, 200, Inf),
                            labels = c("<100", "100-120", "120-150", "150-200", ">200"),
                            right = FALSE),
            levels = c("<100", "100-120", "120-150", "150-200", ">200")))  # Reverse order

# Define wind direction order
wind_plot_df$cardDir <- factor(wind_plot_df$cardDir,
                    levels = c("N", "NE", "E", "SE", "S", "SW", "W", "NW"))

# remove NA values
wind_plot_df <- wind_plot_df[!is.na(wind_plot_df$cardDir), ]
# Custom color palette
custom_colors <- c("#0335258e", "#1a755a8e", "#5b99868e", "#98d6c48e")

# Plot wind rose
advplot <- ggplot(wind_plot_df, aes(x = cardDir, fill = FF_interval)) +
  geom_bar(position = "stack", width = 1, color = "black") +
  coord_polar(start = 0) +
  scale_fill_manual(values = custom_colors, name = "Force (mm/s)") +
  labs(x = "Direction", y = "Count", title = "") +
  theme_minimal() +
  theme(axis.text.x = element_text(size = 12, face = "bold"),
        legend.title = element_text(size = 12, face = "bold"),
        legend.text = element_text(size = 10))

file_path <- paste(im_folder, "optim/comephore/estimated_adv_",
                    df_result$quantile*1000, "_", df_result$delta, "_",
                    df_result$minDist, "_", df_result$tau,
                    "_", df_result$direction, ".png", sep = "")
ggsave(file_path, plot = advplot, width = 20, height = 15, units = "cm")



# VARIOGRAM PLOTS ##############################################################

# compute variogram with parameters
result <- df_result
# convert last column into numeric
result[, c("beta1", "beta2", "alpha1", "alpha2", "eta1", "eta2")] <-
  lapply(result[, c("beta1", "beta2", "alpha1", "alpha2", "eta1", "eta2")],
         as.numeric)

beta1 <- result$beta1
beta2 <- result$beta2
alpha1 <- result$alpha1
alpha2 <- result$alpha2
eta1 <- result$eta1
eta2 <- result$eta2

# Compute advection for each row
adv_df <- wind_df %>%
  mutate(
    adv_x = (abs(vx)^eta1) * sign(vx) * eta2,
    adv_y = (abs(vy)^eta1) * sign(vy) * eta2
  )
  
# remove NA values
adv_df <- adv_df[!is.na(adv_df$adv_x), ]

list_lags <- lapply(seq_along(list_lags), function(i) {
  list_lags[[i]] %>%
    mutate(vx = wind_df$vx[i], vy = wind_df$vy[i])
})

params <- c(beta1, beta2, alpha1, alpha2, adv_df$adv_x, adv_df$adv_y)

list_chi <- lapply(seq_along(list_lags), function(i) {
  theoretical_chi(params, list_lags[[i]], latlon = TRUE)
})

head(list_chi[[1]])

df_chi_all <- bind_rows(list_chi, .id = "episode") 
head(df_chi_all)


library(ggplot2)

df_chi_all$hlagabs <- abs(df_chi_all$hlag)

ggplot(df_chi_all, aes(x = hlagabs, y = vario)) +
  geom_line() +
  facet_wrap(~ tau, scales = "free_x",
      labeller = labeller(tau = function(x) paste0("tau = ", x))) +
  labs(
    title = "Directional variogram",
    x = "Spatial lag",
    y = "Variogram"
  ) +
  theme_minimal()

n_episodes <- length(unique(df_chi_all$episode))
# save the plot
end_filename <- paste(df_result$quantile*1000, "_", min_spatial_dist, "km_",
                     wind_dir, "_delta", delta, ".png", sep = "")
filename <- paste(im_folder, "optim/comephore/variogram_rpareto_all",
                    "_", end_filename, sep = "")
ggsave(filename, width = 20, height = 15, units = "cm")



# Plot only for tau = 10
df_chi_tau10 <- df_chi_all[df_chi_all$tau == 10, ]

ggplot(df_chi_tau10, aes(x = hlagabs, y = vario)) +
  geom_line() +
  labs(
    title = "",
    x = "Spatial lag",
    y = "Variogram"
  ) +
  theme_minimal()

# save the plot
filename <- paste(im_folder, "optim/comephore/variogram_rpareto_tau10",
                    "_", end_filename, sep = "")
ggsave(filename, width = 20, height = 15, units = "cm")


# Plot only for tau = 10
df_chi_tau10 <- df_chi_all[df_chi_all$tau == 0, ]

ggplot(df_chi_tau10, aes(x = hlagabs, y = vario)) +
  geom_line() +
  labs(
    title = "",
    x = "Spatial lag",
    y = "Variogram"
  ) +
  theme_minimal()

# save the plot
filename <- paste(im_folder, "optim/comephore/variogram_rpareto_tau0",
                    "_", end_filename, sep = "")
ggsave(filename, width = 20, height = 15, units = "cm")
