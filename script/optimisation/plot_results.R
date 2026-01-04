# Load libraries and set theme
source("./script/load_libraries.R")
source("./script/optimisation/config.R")

# Get all files in the folder "R"
functions_folder <- "./R"
files <- list.files(functions_folder, full.names = TRUE)
# load all functions in files
invisible(lapply(files, function(f) source(f, echo = FALSE)))
library(latex2exp)

params <- c(0.3, 0.8, 0.2, 0.7) # beta1, beta2, alpha1, alpha2
adv <- c(0.5, 0.2) # advection
ngrid <- 5
temp <- 0:29
s0 <- c(1, 1) # initial condition
t0 <- 0 # time of the first observation
# Configuration
true_param <- c(params, adv)
beta1 <- params[1]
beta2 <- params[2]
alpha1 <- params[3]
alpha2 <- params[4]
spa <- 1:ngrid
nsites <- ngrid^2 # if the grid is squared

# Number of realizations
M <- 100 # number of simulations
m <- 500 # number of extreme episodes
nres <- M * m


# Function to format values correctly, considering decimal places
format_value <- function(x) {
  # Check if x is an integer
  if (x == as.integer(x)) {
    return(sprintf("%d", x))
  } else {
    # Count number of decimal places
    num_decimals <- nchar(sub("^[^.]*\\.", "", as.character(x)))

    if (x >= 1) {
      # If the number is greater than or equal to 1
      if (num_decimals == 1) {
        return(sprintf("%02d", round(x * 10)))
      } else {
        return(sprintf("%03d", round(x * 100)))
      }
    } else {
      # If the number is less than 1
      if (num_decimals == 1) {
        return(sprintf("%02d", round(x * 10)))
      } else {
        return(sprintf("%03d", round(x * 100)))
      }
    }
  }
}
# Apply formatting
param_str <- paste(sapply(params, format_value), collapse = "_")
adv_str <- paste(sapply(adv, format_value), collapse = "_")
s0_str <- paste(sapply(s0, format_value), collapse = "_")
t0_str <- format_value(t0)

# save the result
result_folder <- paste0(data_folder, "optim_results/rpar/")
foldername <- file.path(
  result_folder,
  paste0("rpar_", param_str),
  paste0("sim_", ngrid^2, "s_", length(temp), "t_s0_", s0_str, "_t0_", t0_str),
  s0_type,
  wind_type,
  distance_type,
  eta_type
)
# "../phd_extremes/data/optim_results/rpar//rpar_03_08_02_07_05_02/sim_25s_30t_s0_1_1_t0_0/random_s0/wind/euclidean/free_eta"

name_file <- paste0("/optim_rpar_", M, "simu_", m, "rep_", ngrid^2,
                "s_", length(temp), "t_", param_str, ".csv")


df_bplot <- read.csv(paste0(foldername, name_file), sep = ",")
indices_wout_adv <- c(1, 2, 3, 4, 5, 6)
df_bplot_wout_adv <- df_bplot[, indices_wout_adv]

df_bplot_wout_adv <- stack(df_bplot_wout_adv)

bp_plot <- ggplot(df_bplot_wout_adv, aes(x = ind, y = values)) +
  geom_boxplot() +
  labs(title = "", x = "Parameters", y = "Estimated values") +
  theme_minimal() +
  geom_point(aes(y = true_param[as.numeric(ind)]), color = "red", pch=4) +
  scale_x_discrete(labels = c(TeX("$\\widehat{\\beta}_1$"),
                                TeX("$\\widehat{\\beta}_2$"),
                                TeX("$\\widehat{\\alpha}_1$"),
                                TeX("$\\widehat{\\alpha}_2$"),
                                TeX("$\\widehat{v}_x$"),
                                TeX("$\\widehat{v}_y$"))) 

                                
# df_bplot <- read.csv(paste0(foldername, name_file), sep = ",")
# indices_w_adv <- c(5, 6)
# df_bplot_adv <- df_bplot[, indices_w_adv]

# df_bplot_adv <- stack(df_bplot_adv)

# bp_plot <- ggplot(df_bplot_adv, aes(x = ind, y = values)) +
#   geom_boxplot() +
#   labs(title = "", x = "Parameters", y = "Estimated values") +
#   theme_minimal() +
#   geom_point(aes(y = adv[as.numeric(ind)]), color = "red", pch=4) +
#   scale_x_discrete(labels = c(TeX("$\\widehat{v}_x$"),
#                               TeX("$\\widehat{v}_y$"))) 

# print(bp_plot)

# save plot
foldername <- paste0(im_folder, "optim/rpar/", param_str, "_", adv_str, "/")

if (!dir.exists(foldername)) {
  dir.create(foldername, recursive = TRUE)
}
# setwd("./script")
name_file <- paste0("bp_optim_", M, "simu_", m, "rep_", ngrid^2,
                "s_", length(temp), "t_s0_", s0_str, "_t0_", t0_str, ".png")

ggsave(paste0(foldername, name_file), plot = bp_plot, width = 7, height = 7)

