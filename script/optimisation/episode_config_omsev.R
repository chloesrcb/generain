# remove objects from the environment
rm(list = ls())

# clear console
cat("\014")

source("./script/load_libraries.R")
source("./script/optimisation/config_omsev.R")
ncores <- detectCores() - 1

# Get all files in the folder "R"
functions_folder <- "./R"
files <- list.files(functions_folder, full.names = TRUE)
invisible(lapply(files, function(f) source(f, echo = FALSE)))

library(latex2exp)
library(sf)
library(dplyr)
library(ggplot2)
library(gridExtra)
library(parallel)

# LOAD DATA ####################################################################
# get rain data from omsev
filename_omsev <- paste0(data_folder,
                         "omsev/omsev_5min/rain_mtp_5min_2019_2024.RData")

load(filename_omsev)
rain_omsev <- rain.all5[, c(1, 6:(ncol(rain.all5)-1))]
head(rain_omsev)
filename_loc <- paste0(data_folder,
                           "omsev/loc_rain_gauges.csv")
# get location of each rain gauge
location_gauges <- read.csv(filename_loc)
location_gauges$Station <- c("iem", "mse", "poly", "um", "cefe", "cnrs",
                             "crbm", "archiw", "archie", "um35", "chu1",
                             "chu2", "chu3", "chu4", "chu5", "chu6", "chu7",
                             "cines", "brives", "hydro")


# DISTANCE AND COORDS ##########################################################
nsites <- nrow(location_gauges)
sites_coords <- location_gauges[, c("Longitude", "Latitude")]
rownames(sites_coords) <- location_gauges$Station
sites_coords_sf <- st_as_sf(sites_coords, coords = c("Longitude", "Latitude"),
                            crs = 4326)
sites_coords_sf <- st_transform(sites_coords_sf, crs = 2154)
coords_m <- st_coordinates(sites_coords_sf)
grid_coords_km <- as.data.frame(coords_m / 1000)
colnames(grid_coords_km) <- c("Longitude", "Latitude")
rownames(grid_coords_km) <- rownames(sites_coords)
grid_coords_m <- as.data.frame(coords_m)
colnames(grid_coords_m) <- c("Longitude", "Latitude")
rownames(grid_coords_m) <- rownames(sites_coords)
head(rain_omsev)

# use the same preprocessing as omsev_optim_2024.R
rownames(rain_omsev) <- rain_omsev$dates
rain_omsev <- rain_omsev[, which(colnames(rain_omsev) != "dates")]

# rain_omsev <- rain_omsev[, !(colnames(rain_omsev) %in% c("cines", "hydro", "brives"))]
# location_gauges <- location_gauges[
#   !(location_gauges$Station %in% c("cines", "hydro", "brives")),
# ]
if (!setequal(colnames(rain_omsev), rownames(grid_coords_m))) {
  stop("Station names differ between rain_omsev and grid_coords_m")
}

grid_coords_m <- grid_coords_m[colnames(rain_omsev), , drop = FALSE]
grid_coords_km <- grid_coords_km[colnames(rain_omsev), , drop = FALSE]

# CHOOSE EXTREME EPISODE FOR R-PARETO ##########################################
set_st_excess <- get_spatiotemp_excess(data = rain_omsev, quantile = q,
                                      remove_zeros = TRUE)
# first_ts <- as.POSIXct(rownames(rain_omsev), tz = "UTC")
# starting_year <- year(first_ts)[1]
# starting_year_suffix <- if (starting_year == 2008) "" else paste0("_from", starting_year)

list_s <- set_st_excess$list_s
list_t <- set_st_excess$list_t
list_u <- set_st_excess$list_u
length(list_s)
for (i in seq_along(list_s)) {
  s0 <- list_s[[i]]
  t0 <- list_t[[i]][1]
  u_s0 <- list_u[[i]][1]
  if (rain_omsev[t0, s0] <= u_s0) {
    stop("Excess is not above threshold for s =", s0, " and t =", t0)
  }
}

# rain_omsev <- rain_omsev[, !(colnames(rain_omsev) %in% c("cines", "hydro", "brives"))]
# grid_coords_m <- grid_coords_m[!(rownames(grid_coords_m) %in% c("cines", "hydro", "brives")), ]
# grid_coords_km <- grid_coords_km[!(rownames(grid_coords_km) %in% c("cines", "hydro", "brives")), ]
dmins <- seq(100, 1600, by = 100)  # in m
episode_size <- delta
set_st_excess <- get_spatiotemp_excess(rain_omsev, quantile = q,
                                       remove_zeros = TRUE)
res <- lapply(dmins, function(dm) {
  sel <- get_s0t0_pairs(
    sites_coords = grid_coords_m,
    data = rain_omsev,
    min_spatial_dist = dm,
    episode_size = episode_size,
    set_st_excess = set_st_excess,
    n_max_episodes = 10000,
    latlon = FALSE,
    beta = 0
  )

  data.frame(
    dmin_m = dm,
    n_episodes = nrow(sel)
  )
})

df_tradeoff <- bind_rows(res)

pA <- ggplot(df_tradeoff, aes(x = dmin_m, y = n_episodes)) +
  geom_line(size = 1.1, color = btfgreen) +
  geom_point(size = 2, color = btfgreen) +
  geom_vline(xintercept = min_spatial_dist, linetype = "dashed", color = "red") +
  theme_minimal() +
  labs(
  x = expression(d[min]~"(m)"),
  y = "Number of selected episodes"
  ) + 
  btf_theme

foldername <- paste0(im_folder, "/optim/omsev/choice_config/")
if (!dir.exists(foldername)) {
  dir.create(foldername, recursive = TRUE)
}
filename <- paste0(foldername, "tradeoff_dmin_episodes_delta", delta, ".png")
ggsave(filename, plot = pA, width = 7, height = 5, units = "in", dpi = 300)








library(data.table)
library(dplyr)
library(ggplot2)
library(lubridate)
library(latex2exp)

qs <- c(0.95)
min_spatial_dists <- c(500, 750, 1000, 1200)
deltas <- c(7, 12, 15)

btfgreen <- "#1b9e77"
btf_theme <- theme_minimal(base_size = 14)

# --- Data cleaning ---
rain <- rain[rowSums(!is.na(rain)) > 0, ]

# --- Initialize results table ---
results_summary <- data.frame(
    q = numeric(),
    dmin = numeric(),
    delta = numeric(),
    n_episodes = integer(),
    u_min = numeric(),
    u_max = numeric(),
    stringsAsFactors = FALSE
)

# --- Main loop ---
for (q in qs) {
    for (min_spatial_dist in min_spatial_dists) {
        for (delta in deltas) {

            cat("\n===== Processing for q =", q,
                    ", dmin =", min_spatial_dist,
                    ", delta =", delta, "=====\n")

            # --- Step 1: Excess ---
            set_st_excess <- get_spatiotemp_excess(
                rain, quantile = q, remove_zeros = TRUE
            )

            # --- Step 2: Select episodes ---
            episode_size <- delta
            s0t0_set <- get_s0t0_pairs(
                grid_coords_m, rain,
                min_spatial_dist = min_spatial_dist,
                episode_size = episode_size,
                set_st_excess = set_st_excess,
                n_max_episodes = 10000,
                latlon = FALSE
            )

            if (nrow(s0t0_set) == 0) {
                cat("⚠️ No episode found for these parameters\n")
                results_summary <- rbind(
                    results_summary,
                    data.frame(q = q, dmin = min_spatial_dist, delta = delta,
                                         n_episodes = 0, u_min = NA, u_max = NA)
                )
                next
            }

            selected_points <- s0t0_set %>%
                mutate(t0_date = as.POSIXct(t0_date, tz = "UTC"))
            
            selected_points$t0_date_rounded <- ceiling_date(selected_points$t0_date, "hour")

            datetimes <- unique(selected_points$t0_date)
            datetimes_hour <- unique(selected_points$t0_date_rounded)
            # save datetime list to csv
            datetime_filename <- paste(data_folder, "/omsev/t0_episodes_q", q * 100,
                                    "_delta", delta, "_dmin", min_spatial_dist,
                                    ".csv", sep = "")
            write.csv(data.frame(t0_date = datetimes_hour), datetime_filename, row.names = FALSE)
            # save datetime list to csv
            datetime_filename <- paste(data_folder, "/omsev/t0_5min_episodes_q", q * 100,
                                    "_delta", delta, "_dmin", min_spatial_dist,
                                    ".csv", sep = "")
            write.csv(data.frame(t0_date = datetimes), datetime_filename, row.names = FALSE)

            # --- Statistics ---
            n_episodes <- nrow(selected_points)
            u_min <- min(selected_points$u_s0, na.rm = TRUE)
            u_max <- max(selected_points$u_s0, na.rm = TRUE)

            # --- Add to results table ---
            results_summary <- rbind(
                results_summary,
                data.frame(q = q,
                                     dmin = min_spatial_dist,
                                     delta = delta,
                                     n_episodes = n_episodes,
                                     u_min = u_min,
                                     u_max = u_max)
            )

            # --- (optional) Plots and saving ---
            # ... (your previous ggplot + write.csv code)

            cat("✅ Finished for q =", q,
                    ", dmin =", min_spatial_dist,
                    ", delta =", delta,
                    "→", n_episodes, "episodes\n")
        }
    }
}

# --- Save summary table ---
summary_file <- paste0(data_folder, "/omsev/summary_episodes.csv")
write.csv(results_summary, summary_file, row.names = FALSE)

cat("\n✅ Summary table saved:", summary_file, "\n")
print(results_summary)


# transform in latex table
library(xtable)
latex_table <- xtable(results_summary, digits = 2)
# Combine u_min and u_max into one column formatted as an interval
results_summary$interval <- paste0("[", 
                                   sprintf("%.2f", results_summary$u_min), 
                                   ", ", 
                                   sprintf("%.2f", results_summary$u_max), 
                                   "]")

# Remove old columns
results_summary$u_min <- NULL
results_summary$u_max <- NULL

# Generate LaTeX table
library(xtable)
latex_table <- xtable(results_summary, digits = 2)
print(latex_table, include.rownames = FALSE)

# save latex table to file
table_folder <- "../phd_extremes/thesis/resources/tables/"
latex_file <- paste0(table_folder, "omsev_summary_episodes.tex")

print(latex_table, include.rownames = FALSE, file = latex_file)