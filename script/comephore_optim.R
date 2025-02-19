# remove objects from the environment
rm(list = ls())

# clear console
cat("\014")

# Load libraries and set theme
source("./script/load_libraries.R")


# LOAD DATA ####################################################################
df_comephore <- read.csv("./data/comephore/inside_mtp.csv", sep = ",")
loc_px <- read.csv("./data/comephore/loc_pixels_mtp.csv", sep = ",")

# Take only data after 2007
df_comephore <- df_comephore %>% filter(date > as.Date("2007-12-31"))
# put date in index
rownames(df_comephore) <- df_comephore$date
comephore <- df_comephore[-1] # remove dates column
# Get distances matrix
dist_mat <- get_dist_mat(loc_px)
df_dist <- reshape_distances(dist_mat)
nsites <- nrow(loc_px)




# CHOOSE QUANTILE ##############################################################

q <- 0.998

# BUHL WLSE ####################################################################

# Temporal chi
tmax <- 10
# Temporal chi with spatial lag fixed at 0
# compute chiplot of every site with itself but lagged in time
chimat_dtlag <- temporal_chi(comephore, quantile = q, tmax = tmax,
                             mean = FALSE)

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
filename <- paste(im_folder, "WLSE/comephore/temporal_chi_boxplot_", q,
                    ".png", sep = "")
ggsave(filename, width = 20, height = 15, units = "cm")

# Mean of chi
chimat_dt_mean <- temporal_chi(comephore, tmax, quantile = q,
                               mean = TRUE)
# get h axis in minutes ie x5 minutes
df <- data.frame(lag = c(1:tmax), chi = chimat_dt_mean)
chitemp_plot <- ggplot(df, aes(x = lag, y = chi)) +
  geom_point(color = btfgreen) +
  btf_theme +
  xlab(TeX(r"($\tau$ (hours))")) +
  ylab(TeX(r"($\widehat{\chi}(0,\tau)$)"))

chitemp_plot

wlse_temp <- get_estimate_variotemp(chimat_dt_mean, tmax, nsites,
                                    weights = "exp", summary = T)
print(wlse_temp)
alpha2 <- wlse_temp[[2]]
beta2 <- wlse_temp[[1]]
c2 <- log(beta2)

dftemp <- data.frame(lag = log(df$lag), chi = eta(df$chi))

chitemp_eta_estim <- ggplot(dftemp, aes(x = lag, y = chi)) +
  geom_point(color = btfgreen) +
  btf_theme +
  xlab(TeX(r"($\log(\tau)$)")) +
  ylab(TeX(r"($\eta(\widehat{\chi}(0,\tau))$)")) +
  geom_line(aes(x = lag, y = alpha2 * lag + c2),
            alpha = 0.5, color = "darkred", linewidth = 0.5)

chitemp_eta_estim

# save plot
filename <- paste(im_folder, "WLSE/comephore/temporal_chi_eta_estim_", q,
                    ".png", sep = "")
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

q <- 0.998
chispa_df <- spatial_chi_alldist(df_dist_km, data_rain = comephore,
                                 quantile = q, hmax = hmax)

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
filename <- paste(im_folder, "WLSE/comephore/spatial_chi_", q,
                    ".png", sep = "")
ggsave(filename, width = 20, height = 15, units = "cm")

# WLSE
wlse_spa <- get_estimate_variospa(chispa_df, weights = "exp", summary = T)
print(wlse_spa)

alpha1 <- wlse_spa[[2]]
beta1 <- wlse_spa[[1]]
c1 <- log(beta1)

chispa_eta_estim <- ggplot(etachispa_df, aes(lagspa, chi)) +
  btf_theme +
  geom_point(col = btfgreen) +
  xlab(TeX(r"($\log(h)$)")) +
  ylab(TeX(r"($\eta(\widehat{\chi}(h, 0))$)")) +
  geom_line(aes(x = lagspa, y = alpha1 * lagspa + c1), alpha = 0.6,
            color = "darkred", linewidth = 0.5)

chispa_eta_estim

# save plot
filename <- paste(im_folder, "WLSE/comephore/spatial_chi_eta_estim_", q,
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

# select_extreme_episodes <- function(sites_coords, data, quantile,
#                                     min_spatial_dist, min_time_dist,
#                                     delta, n_max_episodes, beta = 0) {
#   # get maximal value
#   max_value <- max(na.omit(data))
#   max_indices <- which(data == max_value, arr.ind = TRUE)
#   data_temp <- data # copy of the data

#   # initialize the selected points list
#   selected_points <- data.frame(s0 = character(), t0 = integer())

#   # transform the data to a uniform
#   # data_unif <- data
#   # for (i in 1:ncol(data)) {
#   #   data_unif[, i] <- rank(data[, i]) / (nrow(data) + 1)
#   # }
#   # data_temp <- data_unif
#   i <- 0
#   nb_episode <- 0 # number of episodes (selected points)
#   while (nb_episode < n_max_episodes && max_value > quantile) {
#     best_candidate <- NULL

#     for (idx in 1:nrow(max_indices)) {
#       s0_candidate <- colnames(data)[max_indices[idx, 2]]
#       t0_candidate <- max_indices[idx, 1]

#       if (nrow(selected_points) == 0) { # first episode
#         best_candidate <- data.frame(s0 = s0_candidate, t0 = t0_candidate)
#         nb_episode <- nb_episode + 1
#         break
#       } else {
#         # get distances and time differences within the selected points
#         distances <- distHaversine(sites_coords[selected_points$s0, ],
#                                    sites_coords[s0_candidate, ]) / 1000
#         time_differences <- abs(selected_points$t0 - t0_candidate)

#         # check if the candidate is valid given the minimum spatial and
#         # time distances
#         if (all(distances >= min_spatial_dist) ||
#                                       all(time_differences >= min_time_dist)) {
#           best_candidate <- data.frame(s0 = s0_candidate, t0 = t0_candidate)
#           nb_episode <- nb_episode + 1
#           break # stop the loop
#         } else { # not valid
#           # "remove" the candidate from the data to avoid selecting it again
#           data_temp[t0_candidate, s0_candidate] <- -Inf
#         }
#       }
#     }

#     # add the best candidate to the selected points
#     selected_points <- rbind(selected_points, best_candidate)

#     # "remove" the episode from the data
#     t_inf <- best_candidate$t0 - (delta - 1) - beta
#     t_sup <- best_candidate$t0 + (delta - 1) + beta
#     data_temp[t_inf:t_sup, best_candidate$s0] <- -Inf

#     # get the new maximal value
#     max_value <- max(na.omit(data_temp))
#     max_indices <- which(data_temp == max_value, arr.ind = TRUE)
#     i <- i + 1
#   }
#   return(selected_points)
# }

select_extreme_episodes <- function(sites_coords, data, quantile,
                                    min_spatial_dist, min_time_dist,
                                    delta, n_max_episodes, beta = 0) {
  # Convert data to a matrix for fast access
  data <- as.matrix(data)
  # n_sites <- ncol(data)

  # Extract site names
  site_names <- colnames(data)

  # Compute distance matrix (ensure named rows/columns)
  dist_matrix <- as.matrix(distm(sites_coords[site_names, ],
                                        fun = distHaversine)) / 1000

  rownames(dist_matrix) <- site_names
  colnames(dist_matrix) <- site_names

  # Get initial max values
  max_value <- max(na.omit(data))
  max_indices <- which(data == max_value, arr.ind = TRUE)

  # Store selected points
  selected_points <- data.table(s0 = character(), t0 = integer())

  # Logical mask for invalid times
  invalid_time_mask <- matrix(FALSE, nrow(data), ncol(data))

  nb_episode <- 0

  while (nb_episode < n_max_episodes && max_value > quantile) {
    best_candidate <- NULL

    for (idx in seq_len(nrow(max_indices))) {
      s0_candidate <- site_names[max_indices[idx, 2]]
      t0_candidate <- max_indices[idx, 1]

      if (nrow(selected_points) == 0) { # First selection
        best_candidate <- data.table(s0 = s0_candidate, t0 = t0_candidate)
        nb_episode <- nb_episode + 1
        break
      } else {
        # Convert site names to indices for distance lookup
        selected_sites <- selected_points$s0
        valid_indices <- which(site_names %in% selected_sites)

        # Compute distances only if there are previous selections
        if (length(valid_indices) > 0) {
          distances <- dist_matrix[selected_sites, s0_candidate, drop = FALSE]
          time_differences <- abs(selected_points$t0 - t0_candidate)

          # Validate the candidate
          if (all(distances >= min_spatial_dist) || 
                            all(time_differences >= min_time_dist)) {
            best_candidate <- data.table(s0 = s0_candidate, t0 = t0_candidate)
            nb_episode <- nb_episode + 1
            # print(nb_episode)
            break
          } else {
            # Mark the candidate as invalid
            invalid_time_mask[t0_candidate, max_indices[idx, 2]] <- TRUE
          }
        }
      }
    }

    if (!is.null(best_candidate)) {
      # Store selected episode
      selected_points <- rbindlist(list(selected_points, best_candidate))

      # Remove nearby temporal data
      t_inf <- max(1, best_candidate$t0 - (delta - 1) - beta)
      t_sup <- min(nrow(data), best_candidate$t0 + (delta - 1) + beta)
      invalid_time_mask[t_inf:t_sup, 
                        which(site_names == best_candidate$s0)] <- TRUE
    }

    # Update max selection by ignoring invalid positions
    masked_data <- data
    masked_data[invalid_time_mask] <- -Inf
    max_value <- max(masked_data)
    max_indices <- which(masked_data == max_value, arr.ind = TRUE)
  }

  return(selected_points)
}



get_extreme_episodes <- function(selected_points, data, delta, beta = 0) {
  # initialize the list of episodes
  episodes <- list()
  for (i in 1:nrow(selected_points)) {
    # s0 <- selected_points$s0[i]
    t0 <- selected_points$t0[i]
    t_inf <- t0 - (delta - 1) - beta
    t_sup <- t0 + (delta) + beta
    episode <- data[t_inf:t_sup, ] # get the episode
    episodes <- append(episodes, list(episode))
  }
  return(episodes)
}

min_spatial_dist <- 5  # in km
min_time_dist <- 5  # in hours
delta <- 12 # step for the episode before and after the max value

# Get coords
sites_coords <- loc_px[, c("Longitude", "Latitude")]
rownames(sites_coords) <- loc_px$pixel_name

q <- 0.998 # quantile
threshold <- 1000 # initialize the threshold to a high value
# get minimal threshold for q quantile for all columns in comephore
for (i in 1:ncol(comephore)) {
  threshold <- min(threshold, quantile(comephore[, i], probs = q, na.rm = TRUE))
}

# Get the selected episodes
selected_points <- select_extreme_episodes(sites_coords, comephore, threshold,
                                        min_spatial_dist, min_time_dist,
                                        delta = delta, n_max_episodes = 10000)

n_episodes <- length(selected_points$s0)
length(unique(selected_points$s0)) # can be same s0
length(unique(selected_points$t0)) # never same t0

list_episodes <- get_extreme_episodes(selected_points, comephore,
                                      delta = delta)
list_episodes[[100]]
nrow(list_episodes[[1]]) # size of episode

# boxplot(list_episodes[[1]])
s0_list <- selected_points$s0
t0_list <- selected_points$t0

# Compute the lags and excesses for each conditional point
list_results <- mclapply(1:length(s0_list), function(i) {
  s0 <- s0_list[i]
  s0_coords <- sites_coords[s0, ]
  t0 <- t0_time_list[i]
  episode <- list_episodes[[i]]
  ind_t0 <- delta
  lags <- get_conditional_lag_vectors(sites_coords, s0_coords, ind_t0,
                                  tau_max = tmax, latlon = TRUE)
  lags$hx <- lags$hx / 1000  # in km
  lags$hy <- lags$hy / 1000  # in km
  lags$hnorm <- lags$hnorm / 1000  # in km
  excesses <- empirical_excesses(episode, threshold, lags, type = "rpareto",
                                  threshold = TRUE, t0 = ind_t0)
  list(lags = lags, excesses = excesses)
}, mc.cores = detectCores() - 1)

list_lags <- lapply(list_results, `[[`, "lags")
list_excesses <- lapply(list_results, `[[`, "excesses")

# ADD WIND DATA ################################################################

# get wind data
wind_mtp <- read.csv("./data/wind/data_gouv/wind_mtp.csv")

# Convert datetime to POSIXct
wind_mtp$datetime <- as.POSIXct(wind_mtp$datetime,
                                format = "%Y-%m-%d %H:%M:%S", tz = "UTC")

convert_to_cardinal <- function(degrees) {
  if (is.na(degrees)) {
    return(NA)  # Return NA if the input is NA
  }

  directions <- c("N", "NE", "E", "SE", "S", "SW", "W", "NW", "N")
  breaks <- c(0, 22.5, 67.5, 112.5, 157.5, 202.5, 247.5, 292.5, 337.5, 360)

  return(directions[findInterval(degrees, breaks, rightmost.closed = TRUE)])
}

# Apply function to the DD column
wind_mtp$cardDir <- sapply(wind_mtp$DD, convert_to_cardinal)
wind_mtp$cardDir <- as.character(wind_mtp$cardDir)  # Ensure it's character
wind_mtp$cardDir[is.na(wind_mtp$DD)] <- NA
# Check if NA values are properly handled
summary(wind_mtp)

head(wind_mtp$cardDir)

get_mode_dir <- function(x) {
  x <- na.omit(x)  # Remove NA values
  uniq_values <- unique(x)  # Get unique values
  freq_table <- tabulate(match(x, uniq_values))  # Count occurrences
  mode_value <- uniq_values[which.max(freq_table)]  # Most frequent value
  return(mode_value)
}

# Function to compute the mean wind speed and direction for each episode
compute_wind_episode <- function(episode, delta, quantile) {
  timestamps <- as.POSIXct(rownames(episode),
                            format = "%Y-%m-%d %H:%M:%S", tz = "UTC")

  # get column name of the max value
  max_value <- max(episode, na.rm = TRUE)
  s0 <- names(episode)[which(sapply(episode,
                      function(col) max(col, na.rm = TRUE)) == max_value)]

  # Subset the wind data for the episode
  wind_subset <- wind_mtp %>%
    filter(datetime %in% timestamps)

  # Time when there is excess above quantile in s0
  s0_excess_time <- which(episode[, s0] > quantile)
  wind_subset_excess <- wind_mtp %>%
    filter(datetime %in% timestamps[s0_excess_time])

  # Compute the mean wind speed and direction for the episode if there is data
  if (nrow(wind_subset) > 0) {
    FF <- wind_subset$FF[1 + delta]
    cardDir <- get_mode_dir(wind_subset$cardDir)
    DD <- mean(wind_subset$DD[wind_subset$cardDir == cardDir])
    cardDir_excess <- get_mode_dir(wind_subset_excess$cardDir)
    DD_excess <- mean(wind_subset_excess$DD[wind_subset_excess$cardDir ==
                                                            cardDir_excess])
    DD_excess <- mean(wind_subset_excess$DD)
    DD_t0 <- wind_subset$DD[1 + delta]
  } else {
    FF <- NA
    DD <- NA
    cardDir <- NA
    cardDir_excess <- NA
    DD_excess <- NA
    DD_t0 <- NA
  }

  return(data.frame(FF = FF, DD = DD, cardDir = cardDir,
                    cardDir_excess = cardDir_excess,
                    DD_excess = DD_excess, DD_t0 = DD_t0))
}

wind_per_episode <- lapply(list_episodes, compute_wind_episode, delta = delta,
                          quantile = threshold)

wind_ep_df <- do.call(rbind, wind_per_episode)
head(wind_ep_df)

wind_ep_df$cardDirt0 <- sapply(wind_ep_df$DD_t0, convert_to_cardinal)

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
                      min_time_dist, "_h_delta_", delta, ".png", sep = "")
filename <- paste(im_folder, "wind/datagouv/wind_card_dir_", end_filename,
                    sep = "")
ggsave(filename, width = 20, height = 15, units = "cm")


# Ensure cardDir is a factor with correct order
wind_ep_df$cardDirt0 <- factor(wind_ep_df$cardDirt0,
                        levels = c("N", "NE", "E", "SE", "S", "SW", "W", "NW"))

# Plot wind rose
ggplot(wind_ep_df, aes(x = cardDirt0, fill = FF_interval)) +
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

# Mean by excess in episode
wind_ep_df$cardDir_excess <- factor(wind_ep_df$cardDir_excess,
                        levels = c("N", "NE", "E", "SE", "S", "SW", "W", "NW"))

# Plot wind rose
ggplot(wind_ep_df, aes(x = cardDir_excess, fill = FF_interval)) +
  geom_bar(position = "stack", width = 1, color = "black") +
  coord_polar(start = 15 * pi / 8) +
  scale_fill_manual(values = custom_colors, name = "Force (m/s)") +
  labs(x = "Direction in excess (for s0)", y = "Count", title = "") +
  theme_minimal() +
  theme(axis.text.x = element_text(size = 12, face = "bold"),
        legend.title = element_text(size = 12, face = "bold"),
        legend.text = element_text(size = 10))

# save the wind data
filename <- paste(im_folder,
            "wind/datagouv/wind_card_dir_excess_", end_filename, sep = "")
ggsave(filename, width = 20, height = 15, units = "cm")

wind_ep_df$DD_excess <- as.numeric(wind_ep_df$DD_excess)

# Plot wind rose
ggplot(wind_ep_df, aes(x = DD_excess)) +
  geom_bar(color = "#797474", fill = btfgreen, alpha = 0.5, width = 5) +
  coord_polar(start = 15.85 * pi / 8) +
  labs(x = "Mean wind direction in excess for s0 (in degree)", y = "Count",
        title = "") +
  theme_minimal() +
  theme(axis.text.x = element_text(size = 12, face = "bold"))

# save the wind data
filename <- paste(im_folder, "wind/datagouv/wind_dir_excess_", end_filename,
                    sep = "")
ggsave(filename, width = 20, height = 15, units = "cm")

# Plot count of cardinal directions for cardDir and DD_t0
# wind_ep_df$Match <- wind_ep_df$cardDir == wind_ep_df$cardDir_excess

# # Bar plot of matches vs mismatches
# ggplot(wind_ep_df, aes(x = Match)) +
#   geom_bar(aes(fill = Match), color = "#797474", alpha = 0.7, width = 0.5) +
#   scale_fill_manual(values = c("TRUE" = btfgreen, "FALSE" = "tomato")) +
#   labs(x = "Cardinal direction t0 == Most frequent cardinal direction",
#         y = "Count", title = "") +
#   theme_minimal() +
#   theme(axis.text.x = element_text(size = 12, face = "bold"))

# filename <- paste(im_folder, "wind/datagouv/wind_dir_match.png", sep = "")
# ggsave(filename, width = 30, height = 15, units = "cm")

# OPTIMIZATION #################################################################

# Haversine distance between (x1, y1) and (x2, y2)
haversine <- function(x1, y1, x2, y2) {
  R <- 6371  # Earth radius in km
  dx <- (x2 - x1) * pi / 180 # Convert to radians
  dy <- (y2 - y1) * pi / 180 # Convert to radians

  a <- sin(dy / 2)^2 + cos(y1 * pi / 180) * cos(y2 * pi / 180) * sin(dx / 2)^2
  c <- 2 * atan2(sqrt(a), sqrt(1 - a))

  distance <- R * c  # Distance in km
  return(distance)
}

theorical_chi <- function(params, df_lags, wind_vect = NA) {
  beta1 <- params[1]
  beta2 <- params[2]
  alpha1 <- params[3]
  alpha2 <- params[4]
  if (all(is.na(wind_vect))) {
    adv <- params[5:6]
  } else {
    eta1 <- params[5]
    eta2 <- params[6]
    adv <- (abs(wind_vect)^eta1) * sign(wind_vect) * eta2
  }

  chi_df <- df_lags[c("s1", "s2", "tau")]
  # Get vario and chi for each lagtemp in meters
  lat_m_per_deg <- 111132.92 - 559.82 * cos(2 * df_lags$s1y * pi / 180) +
                  1.175 * cos(4 * df_lags$s1y * pi / 180)

  lon_m_per_deg <- 111412.84 * cos(df_lags$s1y * pi / 180) -
                    93.5 * cos(3 * df_lags$s1y * pi / 180)

  chi_df$s1xv <- df_lags$s1x * lon_m_per_deg
  chi_df$s1yv <- df_lags$s1y * lat_m_per_deg
  chi_df$s2xv <- df_lags$s2x * lon_m_per_deg + adv[1] * df_lags$tau
  chi_df$s2yv <- df_lags$s2y * lat_m_per_deg + adv[2] * df_lags$tau

  chi_df$hnormV <- (sqrt((chi_df$s2xv - chi_df$s1xv)^2 +
                        (chi_df$s2yv - chi_df$s1yv)^2)) / 1000

  chi_df$vario <- (2 * beta1) * chi_df$hnormV^alpha1 +
                  (2 * beta2) * abs(chi_df$tau)^alpha2

  chi_df$chi <- 2 * (1 - pnorm(sqrt(0.5 * chi_df$vario)))
  return(chi_df)
}

neg_ll <- function(params, data, df_lags, quantile, excesses, wind_vect = NA,
                      hmax = NA,  s0 = NA, t0 = NA, threshold = FALSE) {
  if (length(params) == 4) {
    params <- c(params, 0, 0)
  } else if (length(params) != 6) {
    stop("The number of initial parameters must be 4 or 6.")
  }

  Tmax <- nrow(data) # number of total observations

  if (all(!is.na(s0))) { # if we have a conditioning location
    p <- 1 # sure excess for r-Pareto process in (s0,t0)
  } else {
    # number of marginal excesses
    nmarg <- get_marginal_excess(data, quantile, threshold)
    p <- nmarg / Tmax # probability of marginal excesses
  }

  # Bounds for the parameters
  lower.bound <- c(1e-08, 1e-08, 1e-08, 1e-08)
  upper.bound <- c(Inf, Inf, 1.999, 1.999)
  if (length(params) == 6 && any(is.na(wind_vect))) {
    lower.bound <- c(lower.bound, -Inf, -Inf)
    upper.bound <- c(upper.bound, Inf, Inf)
  } else if (length(params) == 6 && all(!is.na(wind_vect))) {
    lower.bound <- c(lower.bound, 1e-08, 1e-08)
    upper.bound <- c(upper.bound, Inf, Inf)
  }

  # Check if the parameters are in the bounds
  if (any(params < lower.bound) || any(params > upper.bound)) {
    return(1e50)
  }

  chi <- theorical_chi(params, df_lags, wind_vect) # get chi matrix
  ll_df <- df_lags # copy the dataframe
  ll_df$kij <- excesses$kij # number of excesses
  ll_df$Tobs <- excesses$Tobs
  ll_df$hnormV <- chi$hnormV
  ll_df$chi <- chi$chi
  ll_df$chi <- ifelse(ll_df$chi <= 0, 1e-10, ll_df$chi)
  ll_df$pchi <- 1 - p * ll_df$chi
  ll_df$pchi <- ifelse(ll_df$pchi <= 0, 1e-10, ll_df$pchi)

  # number of non-excesses
  ll_df$non_excesses <- ll_df$Tobs - ll_df$kij
  ll_df$ll <- ll_df$kij * log(ll_df$chi) +
              ll_df$non_excesses * log(ll_df$pchi)
  if (!is.na(hmax)) {
    ll_df <- ll_df[ll_df$hnorm <= hmax, ]
  }

  nll <- -sum(ll_df$ll, na.rm = TRUE)
  return(nll)
}

neg_ll_composite_rpar <- function(params, list_episodes, list_lags, quantile,
                    list_excesses, wind_df, hmax = NA, s0_list = NA,
                    t0_list = NA, threshold = FALSE) {
  if (all(is.na(wind_df))) {
    wind_vect <- NA
  }
  m <- length(list_excesses) # number of r-pareto processes
  nll_composite <- 0 # composite negative log-likelihood
  # print(params)
  for (i in 1:m) {
    # extract lags and excesses from i-th r-pareto process from data
    df_lags <- list_lags[[i]]
    excesses <- list_excesses[[i]]
    episode <- list_episodes[[i]]
    if (!all(is.na(wind_df))) {
      wind_vx <- wind_df$vx[i]
      wind_vy <- wind_df$vy[i]
      wind_vect <- c(wind_vx, wind_vy)
    }
    s0 <- s0_list[i]
    t0 <- t0_list[i]
    nll_i <- neg_ll(params, episode, df_lags, quantile, wind_vect = wind_vect,
                    hmax = hmax, excesses = excesses, s0 = s0, t0 = t0,
                    threshold = threshold)
    nll_composite <- nll_composite + nll_i
  }
  print(nll_composite)
  return(nll_composite)
}


# Get wind vector coordinates  (vx, vy)
wind_ep_df$vx <- wind_ep_df$FF * cos(wind_ep_df$DD_t0)
wind_ep_df$vy <- wind_ep_df$FF * sin(wind_ep_df$DD_t0)
head(wind_ep_df)

wind_df <- wind_ep_df[, c("vx", "vy")]

init_param <- c(beta1, beta2, alpha1, alpha2, 1, 1)
# q <- 1
result <- optim(par = init_param, fn = neg_ll_composite_rpar,
        list_episodes = list_episodes, quantile = threshold,
        list_lags = list_lags,
        list_excesses = list_excesses, hmax = 7, s0_list = s0_list,
        wind_df = wind_df,
        t0_list = t0_list, threshold = TRUE,
        method = "L-BFGS-B",
        lower = c(1e-08, 1e-08, 1e-08, 1e-08, 1e-08, 1e-08),
        upper = c(Inf, Inf, 1.999, 1.999, Inf, Inf),
        control = list(maxit = 10000),
        hessian = F)

estimated_params_noadv <- result$par

# Check the convergence
if (result$convergence != 0) {
  warning("The optimization did not converge.")
}

# Extract the results
df_result <- data.frame(beta1 =  result$par[1],
                        beta2 = result$par[2],
                        alpha1 = result$par[3],
                        alpha2 = result$par[4],
                        eta1 = result$par[5],
                        eta2 = result$par[6])


df_result <- data.frame(beta1 = 6.12,
                        beta2 = 0.14,
                        alpha1 = 0.01,
                        alpha2 = 0.24,
                        eta1 = 1.37,
                        eta2 = 1.25)

colnames(df_result) <- c("beta1", "beta2", "alpha1", "alpha2", "eta1",
                                "eta2")

kable(df_result, format = "latex") %>%
  kable_styling(bootstrap_options = c("striped", "hover", "condensed",
  "responsive"), latex_options = "H")

# VARIOGRAM PLOTS ##############################################################

# compute variogram with parameters
result <- df_result
num_ep <- 200
df_lags_s0t0 <- list_lags[[num_ep]]
wind_s0t0 <- wind_df[num_ep, ]
generate_variogram_plots_rpareto(result, df_lags_s0t0, wind_s0t0)

# save the plot
filename <- paste(im_folder, "optim/comephore/variogram_rpareto_DD_t0_ep", num_ep, 
                "_", end_filename, sep = "")
ggsave(filename, width = 20, height = 15, units = "cm")
