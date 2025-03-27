# Load libraries and set theme
source("./script/load_libraries.R")
devtools::load_all() # load the last version of the package
# load("workspace.RData")

# SIMULATION #################################################################

# Simulayion r-Pareto data
data_folder <- "./data/simulations_rpar/rpar_001_04_15_1_3_1/sim_9s_30t_s0_1_1"

# get all files in the folder
files <- list.files(data_folder, full.names = TRUE)
# keep only 1000 files
files <- files[1:1000]
# length(files)
# print(read.table(files[[69001]], header = TRUE, sep = ","))

# concat all files
df_rpareto <- do.call(rbind, lapply(files, read.table, header = TRUE, sep = ","))

# create a grid 5x5
grid <- expand.grid(1:3, 1:3)
colnames(grid) <- c("Longitude", "Latitude")
rownames(grid) <- colnames(df_rpareto)

filename <- paste0(im_folder, "optim/rpar/rpareto_9s_4_first_episodes.png")
png(filename)
plot(df_rpareto$S1, type = "l", col = btfgreen, xlim=c(1, 120), ylim =c(0, 20),
     xlab = "Time", ylab = "Rainfall (mm)")
abline(v = 30, col = "red", lty = 2)
abline(v=60, col = "blue", lty = 2)
abline(v=90, col = "purple", lty = 2)
abline(v=120, col = "black", lty = 2)
dev.off()


filename <- paste0(im_folder, "optim/rpar/rpareto_25s_4_random_episodes.png")
png(filename)
plot(df_rpareto$S1, type = "l", col = btfgreen, xlim = c(15000, 15000 + 120),
      ylim = c(0, 85),
     xlab = "Time", ylab = "Rainfall (mm)")
abline(v = 15000, col = "#be8ba5", lty = 2)
abline(v = 15000 + 30, col = "red", lty = 2)
abline(v= 15000 + 60, col = "blue", lty = 2)
abline(v= 15000 + 90, col = "purple", lty = 2)
abline(v= 15000 + 120, col = "black", lty = 2)
dev.off()


df_rpareto$S1[1:10]
q <- 0.998 # quantile
min_spatial_dist <- 2 # pixels
delta <- 12 # step for the episode before and after the max value
threshold <- 1
colnames(df_rpareto)
# get quantile for the threshold
quantile(df_rpareto$S1, probs = 0.82, na.rm = TRUE)

# Get the selected episodes
selected_points <- select_extreme_episodes(grid, df_rpareto, 0.82,
                                        min_spatial_dist, delta = delta,
                                        n_max_episodes = 100,
                                        time_ext = 0)
list_episodes_points <- get_extreme_episodes(selected_points, comephore,
                                      delta = delta, unif = FALSE)

list_episodes <- list_episodes_points$episodes


library(ggplot2)
library(reshape2)  # for melting wide data to long format

# Convert matrix to data frame
index <- 21
sort(t0_list)
which(t0_list == t0_list[index]-1)
episode_test <- list_episodes[[index]]
df_episode <- as.data.frame(episode_test)
df_episode$Time <- 1:nrow(df_episode)  # Add a time column
s0_list[index]
# Convert from wide to long format
df_long <- melt(df_episode, id.vars = "Time", variable.name = "Series", value.name = "Value")

ggplot(df_long, aes(x = Time, y = Value, group = Series)) +
  geom_line(color = btfgreen) +
  geom_vline(xintercept = delta + 1, color = "red", linetype = "dashed") + 
  labs(title = "Extreme Episode", x = "Time", y = "Value") +
  annotate("text", x = delta + 1.5, y = 0, label = expression(t[0]),
           color = "red", vjust = 0, size = 7) +
  theme_minimal()

filename <- paste(im_folder, "optim/comephore/extreme_episode", index, "_min", min_spatial_dist,
                  "km_max", tmax, "h_delta_", delta, ".png", sep = "")
filename <- "test.png"
ggsave(filename, width = 20, height = 15, units = "cm")

#plot df_rpareto de S17
ggplot(df_rpareto, aes(y = S17, x = 1:nrow(df_rpareto))) +
  geom_line(color = btfgreen) +
  ylab("Count")

# save the plot
filename <- 'test2.png'
ggsave(filename, width = 20, height = 15, units = "cm")

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
  xlab("Site") +
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


tau_vect <- -10:10
tmax <- max(tau_vect)
# Compute the lags and excesses for each conditional point
list_results <- mclapply(1:length(s0_list), function(i) {
  s0 <- s0_list[i]
  s0_coords <- sites_coords[s0, ]
  # t0 <- t0_list[i]
  episode <- list_episodes_unif[[i]]
  u <- u_list[i]
  ind_t0_ep <- delta + 1  # index of t0 in the episode
  lags <- get_conditional_lag_vectors(sites_coords, s0_coords, ind_t0_ep,
                                  tau_vect, latlon = TRUE)
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

