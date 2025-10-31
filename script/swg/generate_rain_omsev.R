
# remove objects from the environment
rm(list = ls())

# clear console
cat("\014")

# Load libraries and set theme
source("./script/load_libraries.R")

# library(generain)
library(ggplot2)
library(reshape2)
library(animation)
library(RandomFields)
library(RandomFieldsUtils)

# Get all files in the folder "R"
functions_folder <- "./R"
files <- list.files(functions_folder, full.names = TRUE)
# load all functions in files
invisible(lapply(files, function(f) source(f, echo = FALSE)))
library(latex2exp)

# --- paramètres globaux ---
set.seed(123)  # pour reproductibilité

times <- 0:24  # 24 pas de temps (par ex. heures)
coords <- expand.grid(x = 1:10, y = 1:10)  # grille 10x10
colnames(coords) <- c("Longitude", "Latitude")
params_vario <- list(
  beta1 = 0.2, beta2 = 0.1,
  alpha1 = 1.2, alpha2 = 1.0
)

params_margins <- list(
  xi    = rep(0.1, nrow(coords)),
  sigma = rep(1.0, nrow(coords)),
  kappa = rep(0.5, nrow(coords))
)

# --- 10 advections différentes (vx, vy) ---
advections <- list(
  c(0.1,  0.0),
  c(0.1,  0.1),
  c(0.0,  0.1),
  c(-0.1, 0.1),
  c(-0.1, 0.0),
  c(-0.1, -0.1),
  c(0.0,  -0.1),
  c(0.1,  -0.1),
  c(0.2,  0.0),
  c(0.0,   0.2)
)

# --- fonction de simulation d'une série d'épisodes ---
n_episodes <- length(advections)

episodes <- vector("list", n_episodes)

for (i in seq_len(n_episodes)) {
  cat("Simulation épisode", i, "avec adv =", advections[[i]], "\n")
  episodes[[i]] <- sim_episode(
    params_vario = params_vario,
    params_margins = params_margins,
    coords = coords,
    times = times,
    adv = advections[[i]],
    t0 = 0,
    p_wet = 0.05,
    beta_occ = 0.1,
    alpha_occ = 1.2
  )
}

# --- Exemple d’accès au résultat du 1er épisode ---
X1 <- episodes[[1]]$X   # champ de pluie simulé
O1 <- episodes[[1]]$O   # occurrence binaire
