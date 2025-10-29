library(ggplot2)

V_episodes <- wind_df
# --- fonction générique ---
profile_param <- function(param_name, grid, base_params, index,
                          filename_prefix = "profile", wind_df = V_episodes,
                          foldername = im_folder) {
  
  negll_vals <- numeric(length(grid))
  
  for (i in seq_along(grid)) {
    params_test <- base_params
    params_test[index] <- grid[i]
    
    # Choisir si on passe wind_df ou non selon le cas
    if (is.null(wind_df) || is.na(wind_df)[1]) {
      negll_vals[i] <- neg_ll_composite(
        params_test,
        list_lags = list_lags,
        list_episodes = list_episodes,
        list_excesses = list_excesses,
        hmax = hmax,
        latlon = FALSE,
        distance = "lalpha"
      )
    } else {
      negll_vals[i] <- neg_ll_composite(
        params_test,
        list_lags = list_lags,
        list_episodes = list_episodes,
        list_excesses = list_excesses,
        hmax = hmax,
        wind_df = wind_df,
        latlon = FALSE,
        distance = "lalpha"
      )
    }
    
    cat(param_name, "=", grid[i], "-> NegLL =", negll_vals[i], "\n")
  }
  
  df <- data.frame(param = grid, negll = negll_vals)
  min_idx <- which.min(df$negll)
  min_val <- df$param[min_idx]
  
  # --- ggplot ---
  p <- ggplot(df, aes(x = param, y = negll)) +
    geom_line(linewidth = 0.9) +
    geom_point(size = 2) +
    geom_vline(xintercept = min_val, color = "red", linetype = "dashed") +
    labs(
      x = param_name,
      y = "-ll"
    ) +
    btf_theme
  
  # --- sauvegarde PNG ---
  filename <- paste0(foldername, "/", filename_prefix, "_", param_name, ".png")
  ggsave(
    filename = filename,
    plot = p, width = 6, height = 4, dpi = 300
  )
  
  return(list(df = df, plot = p, min = min_val))
}


base_params <- c(0.21, 0.76, 0.37, 0.70, 2, 5)

foldername <- paste0(im_folder, "optim/comephore/profiles/")

if(!dir.exists(foldername)) {
  dir.create(foldername, recursive = TRUE)
}

profiles <- list(
#   list(name = "beta1",  grid = seq(0.01, 1, length.out = 30), index = 1, wind_df = V_episodes),
#   list(name = "beta2",  grid = seq(0.01, 2, length.out = 30), index = 2, wind_df = V_episodes),
#   list(name = "alpha1", grid = seq(0.01, 1, length.out = 30), index = 3, wind_df = V_episodes),
#   list(name = "alpha2", grid = seq(0.01, 1.5, length.out = 30), index = 4, wind_df = V_episodes),
  list(name = "eta1",   grid = seq(0, 5, length.out = 30), index = 5, wind_df = V_episodes),
  list(name = "eta2",   grid = seq(0, 6, length.out = 30), index = 6, wind_df = V_episodes)
)


results <- lapply(profiles, function(p) {
  profile_param(p$name, p$grid, base_params, p$index,
                filename_prefix = "profile", wind_df = p$wind_df, foldername = foldername)
})


adv_profiles <- list(
  list(name = "advx", grid = seq(-2, 2, length.out = 30), index = 5, wind_df = NA),
  list(name = "advy", grid = seq(-1, 1, length.out = 30), index = 6, wind_df = NA)
)

adv_results <- lapply(adv_profiles, function(p) {
  profile_param(p$name, p$grid, base_params, p$index,
                filename_prefix = "profile", wind_df = p$wind_df, foldername = foldername)
})

# get min from adv 
advx_min <- adv_results[[1]]$min
advy_min <- adv_results[[2]]$min

V_episodes_adv <- V_episodes
head(V_episodes_adv)
V_episodes_adv$vx <- advx_min
V_episodes_adv$vy <- advy_min

profiles <- list(
  list(name = "alpha1", grid = seq(0.01, 1, length.out = 30), index = 3, wind_df = V_episodes_adv),
  list(name = "alpha2", grid = seq(0.01, 2, length.out = 30), index = 4, wind_df = V_episodes_adv),
  list(name = "beta1",  grid = seq(0.1, 5, length.out = 30), index = 1, wind_df = V_episodes_adv),
  list(name = "beta2",  grid = seq(0.1, 10, length.out = 30), index = 2, wind_df = V_episodes_adv),
  list(name = "eta1",   grid = seq(0.1, 10, length.out = 30), index = 5, wind_df = V_episodes_adv),
  list(name = "eta2",   grid = seq(0.1, 5, length.out = 30), index = 6, wind_df = V_episodes_adv)
)

foldername <- paste0(im_folder, "optim/omsev/profiles/V_is_advmin/")
if(!dir.exists(foldername)) {
  dir.create(foldername, recursive = TRUE)
}

results_advmin <- lapply(profiles, function(p) {
  profile_param(p$name, p$grid, base_params, p$index,
                filename_prefix = "profile", wind_df = p$wind_df, foldername = foldername)
})
