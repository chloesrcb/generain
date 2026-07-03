rm(list = ls())
cat("\014")

################################################################################
# LOAD LIBRARIES
################################################################################
source("./script/load_libraries.R")
library(sf)
library(dplyr)
library(ggplot2)
library(grid)
library(viridis)
source("./R/spg.R", echo = FALSE)
source("./R/distances.R", echo = FALSE)

################################################################################
# CONFIGURATION TO ADAPT
################################################################################
config <- list(
  study_geom_file = file.path(data_folder, "geometry/verdanson_basin.geojson"),
  cell_m = 100, # resolution of the grid in meters
  nT = 12, # number of time steps
  steps = 0:11, # time steps to simulate
  u_emp = 1, # empirical or random threshold
  s0_pixel_id = "pixel_2100", # conditioning pixel id
  adv = c(100, -200), # advection with correct units (here, in meters per 5 minutes)
  seed = 2026, # random seed for reproducibility
  episode_id = "verdanson_test", # identifier for the episode
  im_folder = im_folder, # folder to save images

  # parameters for the marginal distributions
  params_margins = list(
    p0 = 0.989, # occurrence probability
    xi = 0.244, # tail index EGPD
    sigma = 0.536, # scale parameter EGPD
    kappa = 0.308 # bulk parameter EGPD
  ),

  # parameters for the variogram (see model in README)
  params_vario = list(
    beta1 = 0.48,
    beta2 = 0.77,
    alpha1 = 0.12,
    alpha2 = 0.65
  )
)

################################################################################
# RUN THE SPG SIMULATION
################################################################################
res <- run_spg_episode(
  study_geom_file = config$study_geom_file,
  params_vario = config$params_vario,
  params_margins = config$params_margins,
  adv = config$adv,
  s0_pixel_id = config$s0_pixel_id,
  nT = config$nT,
  steps = config$steps,
  u_emp = config$u_emp,
  cell_m = config$cell_m,
  seed = config$seed,
  im_folder = config$im_folder,
  episode_id = config$episode_id,
  make_plots = FALSE # if TRUE, generates plots for each time step
)