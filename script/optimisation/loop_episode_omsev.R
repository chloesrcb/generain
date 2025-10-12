library(data.table)
library(dplyr)
library(ggplot2)
library(lubridate)
library(latex2exp)

qs <- c(0.9, 0.95, 0.97, 0.98, 0.99)
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
