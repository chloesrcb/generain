# Load libraries and set theme
source("./script/load_libraries.R")

# Get all files in the folder "R"
functions_folder <- "./R"
files <- list.files(functions_folder, full.names = TRUE)
# # load all functions in files
invisible(lapply(files, function(f) source(f, echo = FALSE)))
library(latex2exp)

# 0.5454501 4.8048612 0.1253875 0.7045290
# Fixed variogram parameters
params_base <- c(beta1=0.5455, beta2=4.8049, alpha1=0.1254, alpha2=0.7045)

ep_id <- 1
df_lags <- list_lags[[ep_id]]

# --- Experiment ---
adv_x <- unique(V_episodes$vx)
adv_y <- unique(V_episodes$vy)
adv_grid <- expand.grid(advx = adv_x, advy = adv_y)

# Calculs
results <- apply(adv_grid, 1, function(row) {
  advx <- row["advx"]; advy <- row["advy"]
  params <- c(params_base, advx, advy)
  res <- theoretical_chi(params, df_lags, latlon=FALSE)
  res$advx <- advx
  res$advy <- advy
  return(res[, c("advx","advy","hx","hy","hnormV","chi")])
})

results_df <- do.call(rbind, results)

# plot
ggplot(results_df, aes(x=hnormV, y=chi)) +
  geom_point() +
  labs(title="Influence de advx et advy sur Chi",
       x="hnormV", y="chi") 

theoretical_chi_eta <- function(params, df_lags, V, latlon=FALSE, distance="euclidean") {
  # params: c(beta1, beta2, alpha1, alpha2, eta1, eta2)
  beta1 <- params[1]
  beta2 <- params[2]
  alpha1 <- params[3]
  alpha2 <- params[4]
  eta1  <- params[5]
  eta2  <- params[6]
  
  # V est un vecteur c(Vx, Vy)
  advx <- eta1 * abs(V[1])^eta2 * sign(V[1])
  advy <- eta1 * abs(V[2])^eta2 * sign(V[2])
  adv  <- c(advx, advy)
  
  # Appel à la fonction originale
  return(theoretical_chi(c(beta1, beta2, alpha1, alpha2, adv), df_lags, latlon, distance))
}


# Vecteur de vitesse connu (km/h)
ep_id <- 1
V <- c(W_episodes$wx[ep_id], W_episodes$wy[ep_id]) 
# convert m/s to km/h
V <- V * 3.6

# Base params
params_base <- c(beta1=0.5455, beta2=4.8049, alpha1=0.1254, alpha2=0.7045)

df_lags <- list_lags[[ep_id]]
tau_5min <- 1:10
eta1_vals <- seq(0, 10, 0.5)
eta2_vals <- seq(0, 3, 0.5)

for (tau in tau_5min) {
    tau_fixed <- unique(df_lags$tau)[tau]
    # keep only tau == tau_fixed
    df_lags_sub <- df_lags[df_lags$tau == tau_fixed, ]

    results <- expand.grid(eta1=eta1_vals, eta2=eta2_vals)
    results_out <- do.call(rbind, apply(results, 1, function(row) {
        eta1 <- row["eta1"]; eta2 <- row["eta2"]
        params <- c(params_base, eta1, eta2)
        res <- theoretical_chi_eta(params, df_lags_sub, V, latlon=FALSE)
        res$eta1 <- eta1
        res$eta2 <- eta2
        return(res[, c("eta1","eta2","hx","hy","hnormV","chi")])
    }))

    library(ggplot2)

    params_plot <- round(params_base, 2)
    ggplot(results_out, aes(x=eta1, y=chi, color=factor(eta2), group=eta2)) +
    geom_line() +
    geom_point() +
        labs(
            title = paste(
                "V = (", paste(round(V, 3), collapse = ","), ") m/s",
                "\nparams = (", params_plot["beta1"], ",", params_plot["beta2"],
                ",", params_plot["alpha1"], ",", params_plot["alpha2"], ")"
            ),
            x = expression(eta[1]),
            y = expression(chi),
            color = expression(eta[2])
        ) +
        btf_theme +
        ylim(0,1)

    # Save plot for one episode
    filename <- paste(im_folder, "optim/omsev/diagnose/chi_etas_q",
                    q * 100, "_delta", delta, "_dmin", min_spatial_dist,
                    "_episode", ep_id, "_tau", tau, "_wind.pdf", sep = "")
    ggsave(filename, width = 20, height = 15, units = "cm")
}





