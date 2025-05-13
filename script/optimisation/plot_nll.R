# remove objects from the environment
rm(list = ls())

# clear console
cat("\014")

# LOAD LIBRARIES ###############################################################
library(parallel)
library(abind)
muse <- TRUE

if (muse) {
  # Get the muse folder
  folder_muse <- "/home/serrec/work_rainstsimu/rpareto/"
  setwd(folder_muse)
  # Load libraries and set theme
  source("load_libraries.R")
  source("config.R")
} else {
  # Load libraries and set theme
  source("./script/load_libraries.R")
  source("./script/optimisation/config.R")
}

# Get all files in the folder "R"
functions_folder <- "./R"
files <- list.files(functions_folder, full.names = TRUE)
# load all functions in files
invisible(lapply(files, function(f) source(f, echo = FALSE)))
library(latex2exp)

# SIMULATION ###################################################################
result_folder <- "./data/simulations_rpar/"
if (!dir.exists(result_folder)) {
  dir.create(result_folder, recursive = TRUE)
}

# Configuration
true_param <- c(params, adv)
beta1 <- params[1]
beta2 <- params[2]
alpha1 <- params[3]
alpha2 <- params[4]
spa <- 1:ngrid
nsites <- ngrid^2 # if the grid is squared

# Number of realizations
nres <- M * m

# Apply formatting
param_str <- format_value(true_param)
adv_str <- format_value(adv)
s0_str <- format_value(s0)
t0_str <- format_value(t0)

# Save the data
foldername <- paste0("./data/simulations_rpar/rpar_", param_str, "_", adv_str,
                   "/sim_", ngrid^2, "s_", length(temp), "t_s0_",
                    s0_str, "_t0_", t0_str, "/")

if (!dir.exists(foldername)) {
  dir.create(foldername, recursive = TRUE)
}

sites_coords <- generate_grid_coords(ngrid)

# Save the data
files <- list.files(foldername, full.names = TRUE)
length(files)
list_rpar <- list()
for (i in 1:M) {
  for (j in 1:m) {
    file_name <- paste0(foldername, "rpar_", ngrid^2, "s_",
                        length(temp), "t_simu", i, "_", j, ".csv")
    list_rpar[[((i-1) * m + j)]] <- read.csv(file_name)
  }
}

### Optimization ###############################################################
library(parallel)
num_cores <- detectCores() - 1  # Reserve 1 core for the OS

# Parallel execution
tau_vect <- 0:10
u <- 1
tmax <- max(tau_vect)
sites_coords <- as.data.frame(sites_coords)

# Compute lags and excesses for each episode inside the simulation
list_results <- mclapply(1:nres, function(i) {
  s0_x <- s0_list[[i]]$x
  s0_y <- s0_list[[i]]$y
  s0_coords <- sites_coords[sites_coords$Longitude == s0_x &
                            sites_coords$Latitude == s0_y, ]
  episode <- list_rpar[[i]]
  lags <- get_conditional_lag_vectors(sites_coords, s0_coords, t0,
                                      tau_vect, latlon = FALSE)
  excesses <- empirical_excesses(episode, u, lags, type = "rpareto",
                                 t0 = t0, threshold = TRUE)
  list(lags = lags, excesses = excesses)
}, mc.cores = num_cores)

# Extract lags and excesses from the list of results
list_lags <- lapply(list_results, `[[`, "lags")
list_excesses <- lapply(list_results, `[[`, "excesses")

# With wind data == adv to estimate eta1 and eta2 to 1
wind_df <- data.frame(vx = adv[1], vy = adv[2])
init_param <- c(params[1:4], 1.3, 1)

# plot the negative log-likelihood for different values of eta1 and eta2
eta1_seq <- seq(0.5, 2, length.out = 30)
eta2_seq <- seq(0.5, 2, length.out = 30)

nll_matrix <- matrix(NA, nrow = length(eta1_seq), ncol = length(eta2_seq))
for(i in 1:length(eta1_seq)) {
  for (j in 1:length(eta2_seq)) {
    eta1 <- eta1_seq[i]
    eta2 <- eta2_seq[j]
    params <- c(params[1:4], eta1, eta2)
    
    # Compute the negative log-likelihood
    nll <- neg_ll_composite( 
        params = params,
        list_episodes = list_rpar,
        list_lags = list_lags,
        list_excesses = list_excesses,
        hmax = 7, wind_df = wind_df,
        directional = TRUE,
        latlon = FALSE,,
        rpar = TRUE
        )
    
    nll_matrix[i, j] <- nll
  }
}

# plot
library(ggplot2)
library(reshape2)



# plot the negative log-likelihood for different values of eta1 and eta2
library(plotly)

plot_ly(
  x = eta1_seq, 
  y = eta2_seq, 
  z = ~nll_matrix, 
  type = "surface"
) %>%
  layout(
    title = "Negative Log-Likelihood Surface",
    scene = list(
      xaxis = list(title = "eta1"),
      yaxis = list(title = "eta2"),
      zaxis = list(title = "NLL")
    )
  )

# nll in fonction of eta1 only
eta1_seq <- seq(0.5, 2, length.out = 30)
nll_vector <- numeric(length(eta1_seq))

for(i in 1:length(eta1_seq)) {
  eta1 <- eta1_seq[i]
  params <- c(params[1:4], eta1, 1)
  
  # Compute the negative log-likelihood
  nll <- neg_ll_composite( 
    params = params,
    list_episodes = list_rpar,
    list_lags = list_lags,
    list_excesses = list_excesses,
    hmax = 7, wind_df = wind_df,
    directional = TRUE,
    latlon = FALSE,
    rpar = TRUE
  )
  
  nll_vector[i] <- nll
}
# plot
nll_df <- data.frame(eta1 = eta1_seq, nll = nll_vector)
ggplot(nll_df, aes(x = eta1, y = nll)) +
  geom_line() +
  labs(title = "Negative Log-Likelihood vs eta1",
       x = "eta1",
       y = "Negative Log-Likelihood") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))


# idem for eta2
eta2_seq <- seq(0.5, 2, length.out = 30)
nll_vector <- numeric(length(eta2_seq))

for(i in 1:length(eta2_seq)) {
  eta2 <- eta2_seq[i]
  params <- c(params[1:4], 1, eta2)
  
  # Compute the negative log-likelihood
  nll <- neg_ll_composite( 
    params = params,
    list_episodes = list_rpar,
    list_lags = list_lags,
    list_excesses = list_excesses,
    hmax = 7, wind_df = wind_df,
    directional = TRUE,
    latlon = FALSE,
    rpar = TRUE
  )
  
  nll_vector[i] <- nll
}
# plot
nll_df <- data.frame(eta2 = eta2_seq, nll = nll_vector)
ggplot(nll_df, aes(x = eta2, y = nll)) +
  geom_line() +
  labs(title = "Negative Log-Likelihood vs eta2",
       x = "eta2",
       y = "Negative Log-Likelihood") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))
