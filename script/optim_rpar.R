library(generain)
library(reshape2)
library(ggplot2)
library(dplyr)
library(latex2exp)

btf_theme <- theme_minimal() +
  theme(axis.text.x = element_text(size =  6, angle = 0),
        axis.text.y = element_text(size = 6),
        axis.title.x = element_text(size = 8),
        axis.title.y = element_text(size = 8),
        legend.title = element_text(size = 8),
        legend.text = element_text(size = 6),
        title = element_text(size = 10),
        axis.line = element_blank(),  # Remove axis lines
        panel.border = element_blank(),  # Remove plot border
        panel.background = element_rect(fill = "transparent", color = NA),
        # Remove plot background
        axis.text = element_blank(),  # Remove axis text
        axis.ticks = element_blank(),  # Remove axis ticks
        plot.background = element_rect(fill = "transparent", color = NA),
        panel.grid = element_line(color = "#5c595943"))

# my green color
btfgreen <- "#69b3a2"

# SIMULATION

adv <- c(0.3, 0.2) # advection
params <- c(0.4, 0.2, 1.5, 1) # ok verif sur simu
true_param <- c(params, adv)
beta1 <- params[1]
beta2 <- params[2]
alpha1 <- params[3]
alpha2 <- params[4]
ngrid <- 5
spa <- 1:ngrid
nsites <- ngrid^2 # if the grid is squared
temp <- 1:30

# Conditional point
s0 <- c(1, 1) # y, x ???
t0 <- 1

# Number of realizations
M <- 100
m <- 1000
nres <- M * m

# Simulate the process
set.seed(123)
simu <- sim_rpareto(beta1, beta2, alpha1, alpha2, spa, spa, temp, adv, s0,
                    t0, nres)

if (any(adv < 1 && adv >= 0.1)) {
  adv_int <- adv * 10
  adv_str <- sprintf("%02d_%02d", adv_int[1], adv_int[2])
} else if (adv < 0.1 && adv > 0) {
  adv_int <- adv * 100
  adv_str <- sprintf("%03d_%03d", adv_int[1], adv_int[2])
} else {
  adv_int <- adv
  adv_str <- sprintf("%02d_%02d", adv_int[1], adv_int[2])
}

param_str <- sprintf("%02d_%02d_%02d_%02d", true_param[1] * 10,
                    true_param[2] * 10, true_param[3] * 10, true_param[4] * 10)

s0_str <- sprintf("%01d_%01d", s0[1], s0[2])
# setwd("./script")
# Save the data
foldername <- paste0("./simulations_rpar/rpar_", param_str, "_", adv_str,
                   "/sim_", ngrid^2, "s_", length(temp), "t_s0_",
                    s0_str, "/")


if (!dir.exists(foldername)) {
  dir.create(foldername, recursive = TRUE)
}
save_simulations(simu, ngrid, nres, folder = foldername,
        file = paste0("rpar_", ngrid^2, "s_", length(temp), "t"))

list_simuM <- list()
for (i in 1:nres) {
  file_name <- paste0(foldername, "rpar_", ngrid^2, "s_",
                        length(temp), "t_", i, ".csv")
  list_simuM[[i]] <- read.csv(file_name)
}


# OPTIMISATION
library(parallel)
num_cores <- detectCores() - 1  # Reserve 1 core for the OS

# Parallel execution
sites_coords <- generate_grid_coords(ngrid)
df_lags <- get_conditional_lag_vectors(sites_coords, s0, t0, tau_max = 10)
u <- 1 # threshold corresponding to the r-pareto simulation
result_list <- mclapply(1:M, process_simulation, M = M, m = m,
                        list_simuM = list_simuM, u = u, df_lags = df_lags,
                        s0 = s0, t0 = t0, true_param = true_param,
                        mc.cores = num_cores)

# Combine results into a data frame
df_result_all <- do.call(rbind, result_list)
colnames(df_result_all) <- c("beta1", "beta2", "alpha1",
                              "alpha2", "adv1", "adv2")

df_bplot <- as.data.frame(df_result_all)
df_bplot <- stack(df_bplot)

ggplot(df_bplot, aes(x = ind, y = values)) +
  geom_boxplot() +
  labs(title = "",
    x = "Parameters", y = "Estimated values") +
  theme_minimal() +
  geom_point(aes(y = true_param[as.numeric(ind)]), color = "red", pch=4) 


# save plot
foldername <- "./optim/rpar/"
if (!dir.exists(foldername)) {
  dir.create(foldername, recursive = TRUE)
}
# setwd("./script")
name_file <- paste0("bp_optim_", M, "simu_", m, "rep_", ngrid^2,
                "s_", length(temp), "t_", param_str, "_", 
                adv_str, ".png")

ggsave(paste0(foldername, name_file), width = 5, height = 5)
