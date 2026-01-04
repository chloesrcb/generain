compute_st_variogram_coords <- function(grid,
                                 gamma_space = NULL,
                                 gamma_space_x = NULL,
                                 gamma_space_y = NULL,
                                 gamma_temp,
                                 adv = c(0,0)) {
  sites  <- unique(grid[, c("shifted_x","shifted_y")])
  t_vals <- sort(unique(grid$t))
  
  key_sites <- paste(sites$x, sites$y, sep = "|")
  s_idx <- match(paste(grid$shifted_x, grid$shifted_y, sep = "|"), key_sites)
  t_idx <- match(grid$t, t_vals)

  # Spatial part
  if (!is.null(gamma_space)) {
    if (length(gamma_space) == nrow(grid)) {
      vario <- gamma_space
    } else {
      vario <- gamma_space[s_idx]
    }
  } else {
    if (length(gamma_space_x) == nrow(grid)) {
      vario <- gamma_space_x + gamma_space_y
    } else {
      vario <- (gamma_space_x + gamma_space_y)[s_idx]
    }
  }

  gamma_0 <- vario + gamma_temp[t_idx]
  return(gamma_0)
}


compute_st_gaussian_process_coords <- function(grid, W_s = NULL,
                                        W_s_x = NULL, W_s_y = NULL,
                                        W_t, adv = c(0,0)) {
  sites  <- unique(grid[, c("x","y")])
  t_vals <- sort(unique(grid$t))
  
  key_sites <- paste(sites$x, sites$y, sep = "|")
  s_idx <- match(paste(grid$shifted_x, grid$shifted_y, sep = "|"), key_sites)
  t_idx <- match(grid$t, t_vals)

  if (!is.null(W_s)) {
    W <- W_s[s_idx] + W_t[t_idx]
  } else {
    W <- (W_s_x[s_idx] + W_s_y[s_idx]) + W_t[t_idx]
  }
  
  return(W)
}

# sim_rpareto_coords <- function(beta1, beta2, alpha1, alpha2, coords, t,
#                                adv = c(0, 0), t0 = 0, nres = 1,
#                                random_s0 = FALSE, s0 = coords[1, ],
#                                s0_radius = Inf,
#                                distance = "euclidean", threshold = 1, 
#                                seed = NULL) {
  
#   # distance: "euclidean" or "lalpha" (lalpha = sum of 1D distances)
#   if (!(distance %in% c("euclidean", "lalpha"))) {
#     stop("Invalid distance type. Choose 'euclidean' or 'lalpha'.")
#   }

#   if (!is.null(seed)) {
#     set.seed(seed)
#   }

#   RandomFields::RFoptions(spConform = FALSE, allow_duplicated_locations = TRUE,
#                           install = "no")
  
#   n_sites <- nrow(coords)
#   lt <- length(t)
  
#   # Spatio-temporel grid
#   grid <- expand.grid(site = seq_len(n_sites), t = t)
#   grid$x <- rep(coords$Longitude, times = lt)
#   grid$y <- rep(coords$Latitude,  times = lt)
#   grid$shifted_x <- grid$x - grid$t * adv[1]
#   grid$shifted_y <- grid$y - grid$t * adv[2]
  
#   # Models
#   modelSpace <- RandomFields::RMfbm(alpha = alpha1, var = 2 * beta1)
#   modelTime  <- RandomFields::RMfbm(alpha = alpha2, var = 2 * beta2)
  
#   # Normalize s0 into a dataframe with consistent columns
#   if (is.null(dim(s0))) {
#     s0_df <- data.frame(x = s0[1], y = s0[2])
#   } else {
#     s0_df <- as.data.frame(s0)
#     names(s0_df) <- c("x", "y")
#   }
  
#   # Attach site name if present in coords
#   if ("Site" %in% names(coords)) {
#     site_match <- which(coords$Longitude == s0_df$x & coords$Latitude == s0_df$y)
#     if (length(site_match) == 1) {
#       s0_df$Site <- coords$Site[site_match]
#       s0_df <- s0_df[, c("Site", "x", "y")]
#     }
#   }
  
#   # If random_s0: candidates within radius
#   if (random_s0) {
#     s0_center <- as.numeric(s0_df[1, c("x", "y")])
#     grid_points <- coords[c("Longitude", "Latitude")]
#     colnames(grid_points) <- c("x", "y")
#     dist <- sqrt((grid_points$x - s0_center[1])^2 +
#                    (grid_points$y - s0_center[2])^2)
#     candidate_points <- grid_points[dist <= s0_radius, ]
#     if (nrow(candidate_points) == 0) {
#       stop("No coordinates found within specified radius.")
#     }
#     if ("Site" %in% names(coords)) {
#       candidate_points$Site <- coords$Site[dist <= s0_radius]
#       candidate_points <- candidate_points[, c("Site", "x", "y")]
#     }
#   }
  
#   Z <- array(NA, dim = c(n_sites, lt, nres))
#   s0_used <- vector("list", nres)
  
#   for (i in seq_len(nres)) {
#     # Select conditioning site
#     if (random_s0) {
#       selected_index <- sample(nrow(candidate_points), 1)
#       s0_curr <- candidate_points[selected_index, , drop = FALSE]
#     } else {
#       s0_curr <- s0_df[1, , drop = FALSE]
#     }
#     s0_used[[i]] <- s0_curr
    
#     # Index of conditioning point
#     ind_s0_t0 <- which(grid$x == s0_curr$x & 
#                          grid$y == s0_curr$y & 
#                          grid$t == t0)
    
#     # Temporal variogram
#     gamma_temp <- RandomFields::RFvariogram(modelTime, x = t - t0)
    
#     if (distance == "lalpha") {
#       gamma_space_x <- RandomFields::RFvariogram(modelSpace,
#                         x = grid$shifted_x - grid$shifted_x[ind_s0_t0])
#       gamma_space_y <- RandomFields::RFvariogram(modelSpace,
#                         x = grid$shifted_y - grid$shifted_y[ind_s0_t0])
      
#       gamma_0 <- compute_st_variogram_coords(
#         grid,
#         gamma_space_x = gamma_space_x,
#         gamma_space_y = gamma_space_y,
#         gamma_temp = gamma_temp,
#         adv = adv
#       )
      
#       W_s_x <- RandomFields::RFsimulate(modelSpace, grid$shifted_x, grid = FALSE)
#       W_s_y <- RandomFields::RFsimulate(modelSpace, grid$shifted_y,  grid = FALSE)
#       W_t <- RandomFields::RFsimulate(modelTime, t, n = 1, grid = TRUE)

#       W <- compute_st_gaussian_process_coords(
#         grid, W_s_x = W_s_x, W_s_y = W_s_y, W_t = W_t, adv = adv
#       )
      
#     } else {
#       gamma_space <- RandomFields::RFvariogram(
#         modelSpace,
#         x = grid$shifted_x - grid$shifted_x[ind_s0_t0],
#         y = grid$shifted_y - grid$shifted_y[ind_s0_t0]
#       )
      
#       gamma_0 <- compute_st_variogram_coords(
#         grid,
#         gamma_space = gamma_space,
#         gamma_temp = gamma_temp,
#         adv = adv
#       )
      
#       # W_s <- RandomFields::RFsimulate(modelSpace,
#       #                                 grid$shifted_x, grid$shifted_y,
#       #                                 grid = FALSE)
#       # W_t <- RandomFields::RFsimulate(modelTime, t, n = 1, grid = TRUE)
      
#       # W <- compute_st_gaussian_process_coords(
#       #   grid, W_s = W_s, W_t = W_t, adv = adv
#       # )

#       W_s_all <- RandomFields::RFsimulate(
#               modelSpace,
#               x = grid$shifted_x,
#               y = grid$shifted_y,
#               grid = FALSE
#             )

#       W_s_mat <- matrix(W_s_all, nrow = n_sites, ncol = lt, byrow = FALSE)
#       W_t_vec <- as.numeric(RandomFields::RFsimulate(modelTime, x = t, grid = TRUE))
#       W_t_mat <- matrix(rep(W_t_vec, each = n_sites), nrow = n_sites)

#       W <- W_s_mat + W_t_mat

#     }
    
#     # r-Pareto process
#     Y <- exp(W - W[ind_s0_t0] - gamma_0)
#     R <- evd::rgpd(n = 1, loc = 1, scale = 1, shape = 1)
#     Z[,, i] <- matrix(threshold * R * Y, n_sites, lt, byrow = FALSE)
#   }

#   dimnames(Z) <- list(Site = rownames(coords), Time = t, Realisation = seq_len(nres))
#   if (!is.null(coords$Site)) {
#     rownames(Z) <- coords$Site
#   }

#   return(list(Z = Z, s0_used = s0_used))

# }




save_simu <- function(simu, coords = NULL, prefix = "episode", folder = ".", 
                      simu_id = 1) {
  # simu     : array [n_sites, n_times, n_episodes]
  # coords   : optional data.frame with column 'Site' for site names
  # prefix   : prefix for CSV file names (default "episode")
  # folder   : folder path where to save CSV files (default current folder)
  # simu_id  : index of the simulation (default = 1)
  
  n_sites <- dim(simu)[1]
  n_times <- dim(simu)[2]
  n_episodes <- dim(simu)[3]
  
  res <- vector("list", n_episodes)  # list to store data.frames
  
  # Create folder if it does not exist
  if (!dir.exists(folder)) {
    dir.create(folder, recursive = TRUE)
  }
  
  # Site names: either from coords$Site or generic names
  if (!is.null(coords) && "Site" %in% names(coords)) {
    site_names <- coords$Site
  } else {
    site_names <- paste0("Site_", seq_len(n_sites))
  }
  
  for (m in seq_len(n_episodes)) {
    # Extract episode m
    mat <- simu[ , , m]
    
    # Transpose so rows = time and columns = sites
    df <- as.data.frame(t(mat))
    names(df) <- site_names
    
    # Add time index
    df$time <- seq_len(n_times)
    
    # Move time to first column
    df <- df[, c("time", site_names)]
    
    # Remove time column
    df <- df[, -1]

    # Store in list
    res[[m]] <- df
    
    # File name: include both episode and simulation
    file_name <- file.path(
      folder,
      paste0(prefix, "_", m, "_simu_", simu_id, ".csv")
    )
    write.csv(df, file_name, row.names = FALSE)
  }
  
  return(invisible(res))
}



save_simu <- function(simu, coords = NULL, prefix = "episode", folder = ".", 
                      simu_id = 1, save = FALSE) {
  # simu     : array [n_sites, n_times, n_episodes]
  # coords   : optional data.frame with column 'Site' for site names
  # prefix   : prefix for CSV file names (default "episode")
  # folder   : folder path where to save CSV files (default current folder)
  # simu_id  : index of the simulation (default = 1)
  # save     : if TRUE save CSVs, else return list of data.frames
  
  n_sites <- dim(simu)[1]
  n_times <- dim(simu)[2]
  n_episodes <- dim(simu)[3]
  
  res <- vector("list", n_episodes)  # list to store data.frames
  
  # Site names: either from coords$Site or generic names
  if (!is.null(coords) && "Site" %in% names(coords)) {
    site_names <- coords$Site
  } else {
    site_names <- paste0("S", seq_len(n_sites))
  }
  
  # Create folder only if save = TRUE
  if (save && !dir.exists(folder)) {
    dir.create(folder, recursive = TRUE)
  }
  
  for (m in seq_len(n_episodes)) {
    # Extract episode m
    mat <- simu[ , , m]
    
    # Transpose so rows = time and columns = sites
    df <- as.data.frame(t(mat))
    names(df) <- site_names
    
    # Add time index
    df$time <- seq_len(n_times)
    
    # Move time to first column
    df <- df[, c("time", site_names)]
    
    # Remove time column
    df <- df[, -1]

    # Store in list
    res[[m]] <- df
    
    if (save) {
      # File name: include both episode and simulation
      file_name <- file.path(
        folder,
        paste0(prefix, "_", m, "_simu_", simu_id, ".csv")
      )
      write.csv(df, file_name, row.names = FALSE)
    }
  }
  
  if (save) {
    message("CSV files saved in: ", normalizePath(folder))
    return(invisible(NULL))
  } else {
    return(res)
  }
}


