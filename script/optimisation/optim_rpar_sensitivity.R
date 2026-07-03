
# remove objects from the environment
rm(list = ls())

# clear console
cat("\014")

# LOAD LIBRARIES ###############################################################
library(parallel)
library(abind)
muse <- FALSE
if (muse) {
  # Get the muse folder
  folder_muse <- "/home/serrec/work_rainstsimu/rpareto/"
  setwd(folder_muse)
  # Load libraries and set theme
  source("load_libraries.R")
  im_folder <- "./images"
  source("config.R")
  data_folder <- "./data/"
  num_cores <- 27
} else {
  # Load libraries and set theme
  num_cores <- detectCores() - 1
  source("./script/load_libraries.R")
  source("./script/optimisation/config_rpar.R")
  M <- 50
}

# Get all files in the folder "R"
functions_folder <- "./R"
files <- list.files(functions_folder, full.names = TRUE)
# load all functions in files
invisible(lapply(files, function(f) source(f, echo = FALSE)))
library(latex2exp)

# SIMULATION ###################################################################
result_folder <- paste0(data_folder, "optim_results/rpar/")
if (!dir.exists(result_folder)) {
  dir.create(result_folder, recursive = TRUE)
}


# get random advection for each simulation
set.seed(123)
speed <- exp(rnorm(m, mean = log(1), sd = 0.7))   # large range of speeds
angle <- runif(m, 0, 2*pi)
adv_df <- data.frame(
  adv1 = speed * cos(angle),
  adv2 = speed * sin(angle)
)

set.seed(123)

wind_list <- lapply(1:M, function(i) {
  speed <- exp(rnorm(m, mean = log(1), sd = 0.7))   # large range of speeds
  angle <- runif(m, 0, 2*pi)
  data.frame(
    vx = speed * cos(angle),
    vy = speed * sin(angle)
  )
})

true_param <- c(params, eta1, eta2)
beta1 <- params[1]
beta2 <- params[2]
alpha1 <- params[3]
alpha2 <- params[4]
spa <- 1:ngrid

nsites <- ngrid^2 # if the grid is squared

# Number of realizations
nres <- M * m

# Apply formatting
format_value <- function(x) {
  formatted_values <- sapply(x, function(val) {
    # Check if val is negative
    is_negative <- val < 0
    val <- abs(val)  # Work with the absolute value for formatting

    # Check if val is an integer
    if (val == as.integer(val)) {
      formatted_value <- sprintf("%d", val)
    } else {
      # Count the number of decimals
      num_decimals <- nchar(sub("^[^.]*\\.", "", as.character(val)))

      if (val >= 1) {
        # If the number is greater than or equal to 1
        if (num_decimals == 1) {
          formatted_value <- sprintf("%02d", round(val * 10))
        } else {
          formatted_value <- sprintf("%03d", round(val * 100))
        }
      } else {
        # If the number is less than 1
        if (num_decimals == 1) {
          formatted_value <- sprintf("%02d", round(val * 10))
        } else {
          formatted_value <- sprintf("%03d", round(val * 100))
        }
      }
    }

    if (is_negative) {
      formatted_value <- paste0("neg", formatted_value)
    }

    return(formatted_value)
  })

  # Concatenate values
  return(paste(formatted_values, collapse = "_"))
}
param_str <- format_value(true_param)
# adv_str <- format_value(adv)result_folder
s0_str <- format_value(s0)
t0_str <- format_value(t0)

# Subfolder names
s0_type <- if (random_s0) "random_s0" else "fixed_s0"

eta_type <- if (!is.na(fixed_eta1) && !is.na(fixed_eta2)) {
  "fixed_eta"
} else if (!is.na(fixed_eta1)) {
  "fixed_eta1"
} else if (!is.na(fixed_eta2)) {
  "fixed_eta2"
} else {
  "free_eta"
}

wind_type <- if (use_wind_data) {
  "wind"
} else {
  "no_wind"
}

if (use_wind_data) {
  wind_df <- adv_df
} else {
  wind_df <- NA
  fixed_eta1 <- NA
  fixed_eta2 <- NA
  eta_type <- ""
}

result_folder <- paste0(data_folder, "optim_results/rpar/")
# Folder name to save the data
foldername <- file.path(
  result_folder,
  paste0("rpar_", param_str),
  paste0("sim_", ngrid^2, "s_", length(temp), "t_s0_", s0_str, "_t0_", t0_str),
  s0_type,
  wind_type,
  distance_type,
  eta_type
)

if (!dir.exists(foldername)) {
  message("Folder created: ", foldername)
  dir.create(foldername, recursive = TRUE)
}

tau_vect <- 0:10
u <- 1
sites_coords <- generate_grid_coords(ngrid)
sites_coords <- as.data.frame(sites_coords)

radial_adv_from_wind <- function(wind_df, eta1, eta2, eps = 1e-12) {
  vx <- wind_df$vx
  vy <- wind_df$vy
  r <- sqrt(vx^2 + vy^2)
  scale <- eta1 * (pmax(r, eps)^(eta2 - 1))
  cbind(scale * vx, scale * vy)
}


# Parallel configuration
print("Starting simulations...")

cl <- makeCluster(num_cores)

# Export variables and functions to the cluster
clusterExport(cl, varlist = c(
  "sim_rpareto", "beta1", "beta2", "alpha1", "alpha2",
  "spa", "temp", "adv_df", "t0", "m", "random_s0", "s0",
  "s0_radius", "distance_type",
  "foldername", "generate_grid_coords", "save_simulations",
  "ngrid", "compute_st_gaussian_process",
  "compute_st_variogram", "convert_simulations_to_list",
  "empirical_excesses_rpar", "get_conditional_lag_vectors",
  "t0", "u", "tau_vect", "sites_coords", "eta1", "eta2",
  "radial_adv_from_wind", "wind_list"
))

# Load libraries in the cluster
clusterEvalQ(cl, {
  library(RandomFields)
  library(dplyr)
  library(tidyr)
  library(data.table)
  RandomFields::RFoptions(
    spConform = FALSE,
    allow_duplicated_locations = TRUE,
    storing = FALSE,
    printlevel = 0
  )
})

# Simulate the process in parallel
result_list <- parLapply(cl, 1:M, function(i) {

  wind_df_i <- wind_list[[i]]
  adv_mat_i <- radial_adv_from_wind(wind_df_i, eta1, eta2)  # m x 2

  simu <- sim_rpareto(
    beta1 = beta1, beta2 = beta2,
    alpha1 = alpha1, alpha2 = alpha2,
    x = spa, y = spa, t = temp,
    adv = adv_mat_i,
    t0 = t0,
    nres = m,
    random_s0 = random_s0,
    s0 = s0,
    s0_radius = s0_radius,
    distance = distance_type
  )


  list_rpar <- convert_simulations_to_list(simu$Z, sites_coords)

  list_lags_excesses <- lapply(1:m, function(j) {
    s0_x <- simu$s0_used[[j]][[1]]$x
    s0_y <- simu$s0_used[[j]][[1]]$y
    s0_coords <- sites_coords[sites_coords$Longitude == s0_x &
                              sites_coords$Latitude == s0_y, ]
    lags <- get_conditional_lag_vectors(sites_coords, s0_coords, t0,
                                        tau_vect, latlon = FALSE)
    excesses <- empirical_excesses_rpar(list_rpar[[j]],
                                        df_lags = lags,
                                        t0 = t0, threshold = u)
    lags$tau <- lags$tau
    list(lags = lags, excesses = excesses)
  })

  # Extract lags and excesses from the list of results
  list_lags <- lapply(list_lags_excesses, `[[`, "lags")
  list_excesses <- lapply(list_lags_excesses, `[[`, "excesses")

  return(list(
    list_rpar = list_rpar,
    s0_used = simu$s0_used,
    list_excesses = list_excesses,
    list_lags = list_lags,
    wind_df = wind_df_i
  ))
})
stopCluster(cl) # Stop the cluster

print("Simulations done")
# Get the list of data.frames
list_rpar <- do.call(c, lapply(result_list, function(res) res$list_rpar))
list_excesses <- do.call(c, lapply(result_list,
                                   function(res) res$list_excesses))
list_lags <- do.call(c, lapply(result_list, function(res) res$list_lags))

list_wind <- lapply(result_list, function(res) res$wind_df)


# check excesses
sums_kij <- sapply(list_excesses, function(ex) sum(ex$kij, na.rm = TRUE))
sums_kij


### Optimization: sensitivity to initial parameters ###########################

# True parameters used for simulation
if (use_wind_data) {
  true_param <- c(beta1, beta2, alpha1, alpha2, eta1, eta2)
  par_names <- c("beta1", "beta2", "alpha1", "alpha2", "eta1", "eta2")
  wind_df <- adv_df
} else {
  true_param <- c(beta1, beta2, alpha1, alpha2, adv[1], adv[2])
  par_names <- c("beta1", "beta2", "alpha1", "alpha2", "adv1", "adv2")
  wind_df <- NA
}

names(true_param) <- par_names

# Initial values to test
init_grid <- list(
  true       = true_param,
  plus_10pc  = true_param * 1.10,
  minus_10pc = true_param * 0.90,
  plus_30pc  = true_param * 1.30,
  minus_30pc = true_param * 0.70,
  high_beta  = true_param * c(1.5, 1.5, 1, 1, 1, 1),
  low_beta   = true_param * c(0.5, 0.5, 1, 1, 1, 1),
  high_alpha = true_param * c(1, 1, 1.3, 1.3, 1, 1),
  low_alpha  = true_param * c(1, 1, 0.7, 0.7, 1, 1),
  high_eta   = true_param * c(1, 1, 1, 1, 1.5, 1.5),
  low_eta    = true_param * c(1, 1, 1, 1, 0.5, 0.5)
)

# Safety: alpha and eta2 should stay positive
init_grid <- lapply(init_grid, function(x) {
  x <- pmax(x, 1e-4)
  names(x) <- par_names
  x
})

run_one_init <- function(init_name, init_params) {

  message("Running optimization for init: ", init_name)
  result_list_opt <- mclapply(
    1:M,
    function(i) {
      process_simulation(
        i = i,
        m = m,
        list_simu = list_rpar,
        u = u,
        list_lags = list_lags,
        list_excesses = list_excesses,
        init_params = init_params,
        distance = "euclidean",
        hmax = 7,
        wind_df = list_wind[[i]],
        fixed_eta1 = NA,
        fixed_eta2 = NA,
        normalize = FALSE
      )
    },
    mc.cores = num_cores
  )

  df <- do.call(rbind, result_list_opt)
  df <- as.data.frame(df)
  colnames(df) <- par_names

  df$init_type <- init_name
  df$sim_id <- seq_len(nrow(df))

  for (p in par_names) {
    df[[paste0("init_", p)]] <- init_params[p]
    df[[paste0("true_", p)]] <- true_param[p]
    df[[paste0("bias_", p)]] <- df[[p]] - true_param[p]
    df[[paste0("rel_bias_", p)]] <- (df[[p]] - true_param[p]) / true_param[p]
  }

  df$converged <- !is.na(df[[par_names[1]]])

  df
}

df_result_sensitivity <- do.call(
  rbind,
  Map(run_one_init, names(init_grid), init_grid)
)

print(head(df_result_sensitivity))

### Save results ###############################################################

foldername_result <- file.path(
  data_folder,
  "optim_results/rpar",
  paste0("rpar_", param_str),
  paste0("sim_", ngrid^2, "s_", length(temp), "t_s0_", s0_str, "_t0_", t0_str),
  s0_type,
  wind_type,
  distance_type,
  eta_type,
  "init_sensitivity"
)

if (!dir.exists(foldername_result)) {
  dir.create(foldername_result, recursive = TRUE)
}

name_file <- paste0(
  "optim_rpar_init_sensitivity_",
  M, "simu_", m, "rep_",
  ngrid^2, "s_", length(temp), "t_",
  param_str, ".csv"
)

write.csv(
  df_result_sensitivity,
  file.path(foldername_result, name_file),
  row.names = FALSE
)

message("Results saved in: ", file.path(foldername_result, name_file))

### Summary table ##############################################################

library(dplyr)
library(tidyr)

summary_sensitivity <- df_result_sensitivity %>%
  group_by(init_type) %>%
  summarise(
    n = n(),
    n_converged = sum(converged),
    convergence_rate = mean(converged),
    across(
      all_of(par_names),
      list(
        mean = ~ mean(.x, na.rm = TRUE),
        sd = ~ sd(.x, na.rm = TRUE),
        median = ~ median(.x, na.rm = TRUE)
      ),
      .names = "{.col}_{.fn}"
    ),
    .groups = "drop"
  )

write.csv(
  summary_sensitivity,
  file.path(foldername_result, "summary_init_sensitivity.csv"),
  row.names = FALSE
)

print(summary_sensitivity)



library(dplyr)
library(tidyr)
library(ggplot2)

plot_folder <- file.path(foldername_result, "plots")
if (!dir.exists(plot_folder)) dir.create(plot_folder, recursive = TRUE)

df_long <- df_result_sensitivity %>%
  pivot_longer(
    cols = all_of(par_names),
    names_to = "parameter",
    values_to = "estimate"
  ) %>%
  mutate(
    true_value = true_param[parameter],
    rel_bias = (estimate - true_value) / true_value
  )

# 1) Estimates by initialization
library(ggplot2)

p_est <- ggplot(df_long,
                aes(x = init_type,
                    y = estimate,
                    fill = init_type)) +

  geom_boxplot(width = 0.65,
               alpha = 0.85,
               outlier.size = 0.8,
               linewidth = 0.4) +

  geom_hline(aes(yintercept = true_value),
             colour = "firebrick",
             linetype = "22",
             linewidth = 0.8) +

  facet_wrap(~parameter,
             scales = "free_y",
             ncol = 2) +

  coord_flip() +

  labs(
    x = NULL,
    y = "Estimated value"
  ) +

  theme_minimal(base_size = 13) +
  theme(
    legend.position = "none",
    panel.grid.minor = element_blank(),
    strip.text = element_text(face = "bold", size = 13),
    axis.text.y = element_text(face = "bold"),
    axis.title.y = element_blank(),
    panel.spacing = unit(1.2, "lines")
  )
p_est

library(ggplot2)
library(latex2exp)

latex_params <- c(
  beta1  = TeX("$\\beta_1$"),
  beta2  = TeX("$\\beta_2$"),
  alpha1 = TeX("$\\alpha_1$"),
  alpha2 = TeX("$\\alpha_2$"),
  eta1   = TeX("$\\eta_1$"),
  eta2   = TeX("$\\eta_2$")
)

df_long$parameter_lab <- factor(
  df_long$parameter,
  levels = names(latex_params),
  labels = latex_params
)

ggplot(df_long,
       aes(x = init_type, y = estimate, fill = init_type)) +

  geom_violin(trim = FALSE,
              alpha = 0.35,
              colour = NA) +

  geom_boxplot(width = 0.18,
               outlier.size = 0.6,
               linewidth = 0.35) +

  geom_hline(aes(yintercept = true_value),
             colour = "red3",
             linetype = "dashed",
             linewidth = 0.7) +

  facet_wrap(~parameter_lab,
             scales = "free_y",
             labeller = label_parsed) +

  coord_flip() +


  labs(
    x = "Initial values",
    y = "Estimated value"
  ) +

  theme_minimal(base_size = 13) +
  theme(
    legend.position = "none",
    panel.grid.minor = element_blank(),
    strip.text = element_text(face = "bold", size = 13),
    axis.text.y = element_text(face = "bold"),
    panel.spacing = unit(1.1, "lines")
  )


ggsave(file.path(plot_folder, "boxplot_estimates_by_init.png"),
       p_est, width = 10, height = 7, dpi = 300)



library(ggplot2)
library(latex2exp)

latex_params <- c(
  beta1  = "$\\beta_1$",
  beta2  = "$\\beta_2$",
  alpha1 = "$\\alpha_1$",
  alpha2 = "$\\alpha_2$",
  eta1   = "$\\eta_1$",
  eta2   = "$\\eta_2$"
)
init_labels <- c(
  true        = "Reference",
  plus_10pc   = "+10%",
  minus_10pc  = "-10%",
  plus_30pc   = "+30%",
  minus_30pc  = "-30%",
  high_beta   = expression("High " * beta),
  low_beta    = expression("Low " * beta),
  high_alpha  = expression("High " * alpha),
  low_alpha   = expression("Low " * alpha),
  high_eta    = expression("High " * eta),
  low_eta     = expression("Low " * eta)
)

plot_one_param <- function(param_name, xlim = NULL) {
  
  df_p <- subset(df_long, parameter == param_name)
  
  p <- ggplot(df_p,
              aes(x = init_type, y = estimate, fill = init_type)) +
    
    # geom_violin(trim = FALSE,
    #             alpha = 0.35,
    #             colour = NA) +
    
    geom_boxplot(width = 0.18,
                 outlier.size = 0.6,
                 linewidth = 0.35) +
    
    geom_hline(aes(yintercept = true_value),
               colour = "red3",
               linetype = "dashed",
               linewidth = 0.7) +
    
    coord_flip() +
    
    # scale_fill_brewer(palette = "Set2") +
    scale_x_discrete(labels = init_labels) +
    labs(
      title = " ",
      x = "Initial values",
      y = "Estimated value"
    ) +
    
    theme_minimal(base_size = 13) +
    theme(
      legend.position = "none",
      panel.grid.minor = element_blank(),
      plot.title = element_text(hjust = 0.5, face = "bold", size = 15),
      axis.text.y = element_text(face = "bold")
    )
  
  if (!is.null(xlim)) {
    p <- p + coord_flip(ylim = xlim)
  }
  
  p
}

p_beta1  <- plot_one_param("beta1", xlim = c(0, 1))
p_beta2  <- plot_one_param("beta2", xlim = c(0, 1))
p_alpha1 <- plot_one_param("alpha1", xlim = c(0, 1))
p_alpha2 <- plot_one_param("alpha2", xlim = c(0, 1))
p_eta1   <- plot_one_param("eta1", xlim = c(0, 10))
p_eta2   <- plot_one_param("eta2", xlim = c(0, 10))

p_beta1
p_beta2
p_alpha1
p_alpha2
p_eta1
p_eta2

# SAve plots
plot_folder <- "../phd_extremes/thesis/resources/images/optim_results/rpar/rpar_03_07_01_07_4_2/"
ggsave(file.path(plot_folder, "boxplot_estimates_beta1.png"),
       p_beta1, width = 8, height = 7, dpi = 300)
ggsave(file.path(plot_folder, "boxplot_estimates_beta2.png"),
       p_beta2, width = 8, height = 7, dpi = 300)
ggsave(file.path(plot_folder, "boxplot_estimates_alpha1.png"),
       p_alpha1, width = 8, height = 7, dpi = 300)
ggsave(file.path(plot_folder, "boxplot_estimates_alpha2.png"),
       p_alpha2, width = 8, height = 7, dpi = 300)
ggsave(file.path(plot_folder, "boxplot_estimates_eta1.png"),
       p_eta1, width = 8, height = 7, dpi = 300)
ggsave(file.path(plot_folder, "boxplot_estimates_eta2.png"),
       p_eta2, width = 8, height = 7, dpi = 300)

