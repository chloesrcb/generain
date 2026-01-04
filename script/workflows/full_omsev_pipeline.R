#!/usr/bin/env Rscript

# ---------------------------------------------------------------------------
# Full OMSEV / COMEphore workflow:
# - Run COMEPHORE optimisation (script/optimisation/comephore_muse.R)
# - Run OMSEV optimisation + jackknife (script/optimisation/omsev_jk.R)
# - Compare theoretical vs empirical chi for both datasets
# - Run OMSEV-like simulations and validate dependence structure
# - Compare empirical chi from simulations vs real OMSEV chi
# ---------------------------------------------------------------------------

source("./script/load_libraries.R")

functions_folder <- "./R"
files <- list.files(functions_folder, full.names = TRUE)
invisible(lapply(files, function(f) source(f, echo = FALSE)))

# ---------------------------------------------------------------------------

workflow_config <- list(
  comephore = list(
    script = "./script/optimisation/comephore_pipeline.R",
    adv_group = "significant",
    tau_fixed = 0,
    h_breaks = seq(0, 20, by = 1)
  ),
  omsev = list(
    script = "./script/optimisation/omsev_pipeline.R",
    tau_fixed = 0
  ),
  simulation = list(
    tau_vect = 0:10,
    Nsim = 400,
    adv_vector = c(vx = 0.2, vy = 0.1)
  )
)

run_stage <- function(script_path, label) {
  message("============================================================")
  message("Running stage: ", label)
  env <- new.env()
  sys.source(script_path, envir = env)
  env
}

prepare_wind_df <- function(df_or_mat, cols = c("adv_x", "adv_y")) {
  if (all(cols %in% colnames(df_or_mat))) {
    wind <- as.data.frame(df_or_mat[, cols])
  } else {
    wind <- as.data.frame(df_or_mat)
  }
  colnames(wind) <- c("vx", "vy")
  wind
}

compute_group_chi <- function(list_lags,
                              list_excesses,
                              wind_df,
                              params,
                              tau_min = 0,
                              tau_fixed = 0,
                              h_breaks = seq(0, 10, by = 1),
                              latlon = FALSE,
                              adv_transform = TRUE) {
  plot_th_emp_chi(
    list_lags_filtered = list_lags,
    list_excesses_filtered = list_excesses,
    wind_df_filtered = wind_df,
    params_estimates = params,
    tau_min = tau_min,
    tau_fixed = tau_fixed,
    h_breaks = h_breaks,
    latlon = latlon,
    adv_transform = adv_transform
  )
}

summary_correlation <- function(res_cmp) {
  cor(res_cmp$chi_emp_bar, res_cmp$chi_theo_bar,
      use = "complete.obs", method = "spearman")
}

save_plot <- function(plot_obj, filename, width = 7, height = 6) {
  dir.create(dirname(filename), recursive = TRUE, showWarnings = FALSE)
  ggplot2::ggsave(filename, plot = plot_obj, width = width, height = height)
}

# ---------------------------------------------------------------------------
# 1. COMEPHORE optimisation and diagnostics
# ---------------------------------------------------------------------------
come_env <- run_stage(workflow_config$comephore$script, "COMEPHORE optimisation")

available_groups <- names(come_env$results_all_classes)
if (!length(available_groups)) {
  stop("COMEPHORE pipeline did not return any advection classes.")
}

adv_group <- workflow_config$comephore$adv_group
if (!adv_group %in% available_groups) {
  warning("Requested adv_group '", adv_group,
          "' not found. Using default group ", available_groups[1])
  adv_group <- available_groups[1]
}

class_obj <- come_env$results_all_classes[[adv_group]]
indices_come <- class_obj$indices

come_lags <- come_env$list_lags[indices_come]
come_excesses <- come_env$list_excesses[indices_come]
come_wind <- prepare_wind_df(
  come_env$selected_episodes[indices_come, c("adv_x", "adv_y")]
)
come_params <- class_obj$par

h_breaks <- seq(0, 20, by = 1)
chi_come <- compute_group_chi(
  list_lags = come_lags,
  list_excesses = come_excesses,
  wind_df = come_wind,
  params = come_params,
  tau_fixed = workflow_config$comephore$tau_fixed,
  h_breaks = h_breaks,
  adv_transform = TRUE
)

corr_come <- summary_correlation(chi_come$res_cmp)
message(sprintf("COMEPHORE chi Spearman correlation: %.3f", corr_come))

chi_come$res |>
  dplyr::filter(tau == 0, h == 0) |>
  dplyr::summarise(n = dplyr::n(),
                   emp = unique(chi_emp),
                   theo = unique(chi_theo))

res <- chi_come$res

# plot chi_come emp vs theo
print(chi_come$plots$all)
res_com_cmp <- chi_come$res_cmp
ggplot(res_com_cmp, aes(x = chi_theo_bar, y = chi_emp_bar, size = n_pairs)) +
  geom_point(alpha = 0.7, color = btfgreen) +
  geom_abline(linetype = "dashed", color = "red") +
  theme_bw() +
  scale_size_continuous(range = c(1, 4)) +
  labs(x = "Theoretical Chi", y = "Empirical Chi") +
  btf_theme

# save plot
foldername <- paste0(im_folder, "workflows/full_pipeline/")
if (!dir.exists(foldername)) {
  dir.create(foldername, recursive = TRUE)
}
filename <- paste0(foldername, "comephore_chi_th_emp.png")
ggsave(filename,
       width = 7,
       height = 7)


# remove NA hbin 
res_com_cmp <- res_com_cmp[!is.na(res_com_cmp$hbin), ]
res_com_cmp |>
  dplyr::mutate(diff = chi_emp_bar - chi_theo_bar) |>
  ggplot(aes(x = factor(hbin), y = diff)) +
  geom_boxplot(fill=btfgreen, alpha=0.5) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  theme_bw() +
  labs(
    x = "Spatial lag bin (km)",
    y = expression(hat(chi) - chi)
  ) +
  ylim(-0.1,0.1) +
  btf_theme

# save plot
filename <- paste0(foldername, "comephore_chi_diff_boxplot.pdf")
ggsave(filename,
       width = 7,
       height = 5)




# ---------------------------------------------------------------------------
# 2. OMSEV optimisation, jackknife, and diagnostics
# ---------------------------------------------------------------------------
omsev_env <- run_stage(workflow_config$omsev$script, "OMSEV optimisation & jackknife")
omsev_params <- c(omsev_env$result$par, omsev_env$eta1_class, omsev_env$eta2_class)
omsev_wind <- prepare_wind_df(omsev_env$V_episodes_filtered)


dist_mat <- get_dist_mat(location_gauges) / 1000
df_dist <- reshape_distances(dist_mat)
n_hbins <- 10
h_all <- df_dist$value

h_breaks_omsev <- quantile(
  h_all,
  probs = seq(0, 1, length.out = n_hbins + 1),
  na.rm = TRUE
)

h_breaks_omsev <- unique(as.numeric(h_breaks_omsev))
h_breaks_omsev[length(h_breaks_omsev)] <- 1.6

chi_omsev <- compute_group_chi(
  list_lags = omsev_env$list_lags_filtered,
  list_excesses = omsev_env$list_excesses_filtered,
  wind_df = omsev_wind,
  params = omsev_params,
  tau_fixed = workflow_config$omsev$tau_fixed,
  h_breaks = h_breaks_omsev,
  adv_transform = TRUE
)

corr_omsev <- summary_correlation(chi_omsev$res_cmp)
message(sprintf("OMSEV chi Spearman correlation: %.3f", corr_omsev))

jackknife_results <- omsev_env$jackknife_monthyear_results

# plot chi_come emp vs theo
print(chi_omsev$plots$all)
res_om_cmp <- chi_omsev$res_cmp
ggplot(res_om_cmp, aes(x = chi_theo_bar, y = chi_emp_bar, size = n_pairs)) +
  geom_point(alpha = 0.7, color = btfgreen) +
  geom_abline(linetype = "dashed", color = "red") +
  theme_bw() +
  scale_size_continuous(range = c(1, 4)) +
  labs(x = "Theoretical Chi", y = "Empirical Chi") +
  btf_theme

# save plot
foldername <- paste0(im_folder, "workflows/full_pipeline/")
if (!dir.exists(foldername)) {
  dir.create(foldername, recursive = TRUE)
}
filename <- paste0(foldername, "omsev_chi_th_emp_above02.png")
ggsave(filename,
       width = 7,
       height = 7)


# remove NA hbin 
res_om_cmp <- res_om_cmp[!is.na(res_om_cmp$hbin), ]
res_om_cmp |>
  dplyr::mutate(diff = chi_emp_bar - chi_theo_bar) |>
  ggplot(aes(x = factor(hbin), y = diff)) +
  geom_boxplot(fill=btfgreen, alpha=0.5) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  theme_bw() +
  labs(
    x = "Spatial lag bin (km)",
    y = expression(hat(chi) - chi)
  ) +
  btf_theme

# save plot
filename <- paste0(foldername, "omsev_chi_diff_boxplot_all.pdf")
ggsave(filename,
       width = 7,
       height = 5)



# ---------------------------------------------------------------------------
# 3. Simulation-based validation
# ---------------------------------------------------------------------------
sim_env <- new.env()
sys.source("./script/swg/sim_validate_chi.R", envir = sim_env)
sim_cfg <- utils::modifyList(sim_env$sim_config, workflow_config$simulation)
sim_inputs <- sim_env$prepare_validation_inputs(sim_cfg)
validation_sim <- sim_env$simulate_chi_validation(sim_inputs, sim_cfg)
attr(validation_sim$res_cmp, "correlation") <- summary_correlation(validation_sim$res_cmp)

# ---------------------------------------------------------------------------
# 4. Combined diagnostics and outputs
# ---------------------------------------------------------------------------


plot_sim_vs_obs <- ggplot(
  dplyr::bind_rows(
    dplyr::mutate(chi_omsev$res_cmp, dataset = "OMSEV"),
    dplyr::mutate(validation_sim$res_cmp, dataset = "Simulations")
  ),
  aes(x = chi_theo_bar, y = chi_emp_bar, color = dataset)
) +
  geom_point(alpha = 0.7) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  theme_minimal() +
  labs(
    title = "Empirical chi: OMSEV vs simulations",
    x = "Theoretical Chi",
    y = "Empirical Chi",
    color = ""
  )

output_folder <- file.path(im_folder, "workflows/full_pipeline")
dir.create(output_folder, recursive = TRUE, showWarnings = FALSE)

# save_plot(chi_come$plots$all, file.path(output_folder, "comephore_chi.png"))
# save_plot(chi_omsev$plots$all, file.path(output_folder, "omsev_chi.png"))
# save_plot(validation_sim$plots$all, file.path(output_folder, "simulated_chi.png"))
# save_plot(plot_combined, file.path(output_folder, "combined_chi.png"))
# save_plot(plot_sim_vs_obs, file.path(output_folder, "sim_vs_obs_chi.png"))

# pipeline_results <- list(
#   comephore = list(
#     env = come_env,
#     chi = chi_come,
#     correlation = corr_come
#   ),
#   omsev = list(
#     env = omsev_env,
#     chi = chi_omsev,
#     correlation = corr_omsev,
#     jackknife = jackknife_results
#   ),
#   simulations = list(
#     config = sim_cfg,
#     validation = validation_sim,
#     correlation = attr(validation_sim$res_cmp, "correlation")
#   ),
#   plots = list(
#     comephore = chi_come$plots$all,
#     omsev = chi_omsev$plots$all,
#     simulations = validation_sim$plots$all,
#     combined = plot_combined,
#     sim_vs_obs = plot_sim_vs_obs
#   )
# )

# message("Full workflow completed. Results saved to ", output_folder)
