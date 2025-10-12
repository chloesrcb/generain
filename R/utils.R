#' Formats p-values with asterisks
#'
#'This function formats coefficient estimates by adding asterisks
#' to indicate significance.
#'
#' @param estimate A numeric vector of coefficient estimates.
#' @param p_value A numeric vector of p-values corresponding to the
#'                 coefficient estimates.
#'
#' @return A character vector of formatted coefficient estimates
#'          with asterisks indicating significance.
#'
#' @export
format_estimate <- function(estimate, p_value) {
  asterisks <- ifelse(p_value < 0.001, "***",
                  ifelse(p_value < 0.01, "**",
                    ifelse(p_value < 0.05, "*", "")))
  paste0(estimate, asterisks)
}


#' Save a LaTeX table to a file
#'
#' This function saves a LaTeX table to a file.
#'
#' @param table_latex The LaTeX code representing the table.
#' @param filename The name of the file to save the table to. If not provided,
#'                 the table will not be saved to a file.
#' @param pval A logical value indicating whether to include p-values
#'             in the table. Default is FALSE.
#' @param tab_folder The folder where the table will be saved.
#'                   Default is '/tables/'.
#'
#' @return None
#'
#' @export
save_tablatex <- function(table_latex, filename = "", pval = FALSE,
                    tab_folder = "/tables/") {
    # Add custom LaTeX code for font size
    # font_size <- "\\small"
    # modified_latex <- paste0(font_size, "\n\n",  table_latex)
    # Modify the LaTeX code for the pvalues
    if (pval) {
        modified_latex <- gsub("\\*\\*\\* ", "$^{***}$", table_latex)
        modified_latex <- gsub("\\*\\* ", "$^{**}$", modified_latex)
        modified_latex <- gsub("\\* ", "$^{*}$", modified_latex)
        # Add pvalues information
        complete_latex <- sprintf(
                "\\begin{table}\n    \\centering\n%s\n    \\vspace{-0.5em}\n
                \\footnotesize{$^{*}$p-value$<0.05$; $^{**}$p-value$<0.001$; 
                $^{***}$p-value$<0.001$}\n\\end{table}",
                modified_latex
                )
    } else {
        complete_latex <- table_latex
    }

    # Write the LaTeX code to a file
    writeLines(complete_latex, paste0(tab_folder, filename))

}

#' wlse_tablatex function
#'
#' This function generates a LaTeX table from a summary dataframe of WLSE
#' results.
#'
#' @param summary_df The summary dataframe containing the data to be converted
#'                   into a table.
#' @param ind The index of the summary dataframe to be used for generating the
#'            table.
#' @param filename The name of the file to save the generated LaTeX table.
#'                 If not provided, the table will be printed to the console.
#'
#' @return The generated LaTeX table.
#'
#' @import kableExtra
#' @import knitr
#' @import dplyr
#'
#' @export
wlse_tablatex <- function(summary_df, ind, filename = "") {
    # Extract p-values
    p_values <- summary_df$coefficients[, "Pr(>|t|)"]

    df_wls_temp <- data.frame(summary_df$coefficients)
    df_table <- round(df_wls_temp, 3)

    # Apply formatting to coefficient estimates
    df_table$Estimate <- mapply(format_estimate, df_table$Estimate, p_values)

    table_latex <- kable(df_table[, 1:2], format = "latex",
                     caption = "Parameters estimation with WLSE",
                     col.names = c("Estimate", "Std. Error"),
                     align = "cccc", booktabs = TRUE, escape = FALSE) %>%
        # Use booktabs styling for a more professional appearance
        row_spec(0, bold = FALSE) %>%
        # Combine rows for the first column
        collapse_rows(columns = 1:2, valign = "top", latex_hline = "linespace")

    # Modify the LaTeX code for the row names
    modified_latex <- gsub("\\(Intercept\\)",
                            "$\\\\widehat\\{c\\}_i$",
                            table_latex)
    modified_latex <- gsub("lagtemp",
                        "$\\\\widehat\\{\\\\alpha\\}_i$",
                        modified_latex)
    modified_latex <- gsub("_i", paste0("_", ind), modified_latex)
    save_tablatex(modified_latex, pval = TRUE, filename = filename)
}


#' Create a GIF from a simulation
#'
#' This function creates a GIF from a simulation.
#'
#' @param simulation_data The data from the simulation.
#' @param sites_coords The coordinates of the sites.
#' @param params The parameters used for the simulation.
#' @param foldername The folder where the GIF will be saved.
#' @param type The type of simulation.
#' @param forcedtemp The number of time steps to include in the GIF.
#' @param s0 The conditional site.
#' @param threshold The threshold value to indicate on the plot. Default is 1.
#'
#' @return None, but saves the GIF to a file.
#'
#' @import ggplot2
#' @import reshape2
#' @import animation
#'
#' @export
create_simu_gif <- function(simulation_data, sites_coords, params,
                            foldername, type = "rpar", forcedtemp = NA, 
                            s0 = NULL, threshold = 1) {
  ngrid <- sqrt(ncol(simulation_data)) # Number of grid points in each dimension
  Tmax <- nrow(simulation_data) # Number of time steps
  if (is.na(forcedtemp)) {
    temp <- 1:Tmax
  } else {
    temp <- 1:forcedtemp
    Tmax <- forcedtemp
  }

  beta1 <- params[1]
  beta2 <- params[2]
  alpha1 <- params[3]
  alpha2 <- params[4]
  if (length(params) == 6) {
    adv <- params[5:6]
  } else {
    adv <- c(0, 0)
  }

  simulation_data$Time <- rownames(simulation_data) # Add a time column
  simulation_data_long <- melt(simulation_data) # Convert to long format
  simulation_data_long$Time <- as.numeric(simulation_data_long$Time)

  # Create a dataframe to represent grid points
  grid <- sites_coords
  if (dim(grid)[2] == 3) {
    grid <- grid[, 2:3] # Keep only x and y coordinates
  }
  colnames(grid) <- c("x", "y")

  plots <- list()
  cropped_data <- simulation_data_long[simulation_data_long$Time %in% temp, ]

  max_val <- max(cropped_data$value)
  min_val <- min(cropped_data$value)
  # for each time step
  for (i in unique(cropped_data$Time)) {
    # Add the simulated values to the grid dataframe
    grid$value <- cropped_data$value[cropped_data$Time == i]
    grid$above <- grid$value > threshold

    # Plot
    p <- ggplot(data = grid, aes(x = y, y = x, fill = value)) +
      geom_tile() +
      scale_fill_gradientn(
        colours = c("#70a7ae", "#b8967a", "#9d503d"),
        values = scales::rescale(c(min_val, threshold, max_val)),
        name = "Rainfall in mm",
        limits = c(min_val, max_val)
      ) +
      labs(title = paste0("t =", i - 1, " | Betas: ", beta1, ", ", beta2,
                          " | Alphas: ", alpha1, ", ", alpha2, " | Advection: ",
                          adv[1], ", ", adv[2])) +
      theme_minimal() +
      theme(plot.background = element_rect(fill = "#F9F8F6", color = "#F9F8F6"),
            panel.border = element_blank(),
            panel.grid = element_blank(),
            axis.text = element_blank(),
            axis.title = element_blank())

    if (!is.null(s0)) {
      p <- p +
        geom_text(data = data.frame(x = s0[1], y = s0[2]),
          mapping = aes(x = x, y = y),
          label = "s0", color = "black", size = 5, fontface = "bold", 
          inherit.aes = FALSE)
    }

    plots[[i]] <- p
  }

  param_str <- format_value(params)
  # Save the plots as a gif
  ani.options(interval = 0.5) # time between frames
  saveGIF({
    for (i in temp) {
      print(plots[[i]])
    }
  }, movie.name = file.path(foldername, paste0(type, "_", param_str, "_",
                                   ngrid^2, "s_", Tmax, "t", ".gif")),
  ani.width = 700, ani.height = 600, ani.units = "px", ani.type = "cairo",
  ani.dev = "png")

}

#' Generate data for a specific tau value for the variogram plot
#'
#' This function generates data for a specific tau value for the variogram plot.
#'
#' @param tau The value of tau to generate data for.
#' @param empirical_df The dataframe of empirical variogram values.
#' @param theorical_df The dataframe of theoretical variogram values.
#'
#' @return A dataframe containing the variogram values for the specified tau.
#'
#' @export
generate_data_for_tau <- function(tau, empirical_df, theorical_df) {
  empirical_data <- data.frame(
    h_lag = empirical_df$hnorm[empirical_df$tau == tau],
    variogram = empirical_df$vario[empirical_df$tau == tau],
    chi = empirical_df$chi[empirical_df$tau == tau],
    type = 'Empirical',
    tau = tau
  )

  theorical_data <- data.frame(
    h_lag = theorical_df$hnorm[theorical_df$tau == tau],
    variogram = theorical_df$vario[theorical_df$tau == tau],
    chi = theorical_df$chi[theorical_df$tau == tau],
    type = 'Theoretical',
    tau = tau
  )

  # Combiner les deux dataframes
  return(rbind(empirical_data, theorical_data))
}



#' Generate variogram plots
#'
#' This function generates theoretical and empirical variogram plots for
#' multiple tau values.
#'
#' @param result The result of the variogram estimation.
#' @param df_lags The dataframe of lags.
#' @param true_param The true parameters used to generate the data.
#' @param tau_values The values of tau to generate plots for.
#' @param chi A logical value indicating whether to plot the chi values.
#'            Default is FALSE.
#' @param latlon A logical value indicating whether to use latitude and
#'              longitude coordinates. Default is FALSE.
#' @param distance The type of spatial norm, "euclidean" or "lalpha".
#'                 Default is "euclidean".
#'
#' @return None, but generates plots.
#'
#' @import ggplot2
#'
#' @export
generate_variogram_plots <- function(result, df_lags, true_param, tau_values,
                                     chi = FALSE, latlon = FALSE,
                                     distance = "euclidean") {
  # Get the estimated variogram
  if(typeof(result) == "list" && "beta1" %in% colnames(result)){
    beta1_hat <- result$beta1
    beta2_hat <- result$beta2
    alpha1_hat <- result$alpha1
    alpha2_hat <- result$alpha2
    adv1_hat <- result$adv1
    adv2_hat <- result$adv2
  } else {
    beta1_hat <- result$par[1]
    beta2_hat <- result$par[2]
    alpha1_hat <- result$par[3]
    alpha2_hat <- result$par[4]
    adv1_hat <- result$par[5]
    adv2_hat <- result$par[6]
  }

  empirical_df <- theoretical_chi(c(beta1_hat, beta2_hat, alpha1_hat, 
                                  alpha2_hat, adv1_hat, adv2_hat), df_lags,
                                  latlon, distance)
  theorical_df <- theoretical_chi(true_param, df_lags, latlon, distance)

  combined_data <- do.call(rbind, lapply(tau_values, generate_data_for_tau,
                            empirical_df, theorical_df))

  if (chi) {
    ggplot(combined_data, aes(x = h_lag, y = chi, color = type,
                              linetype = type)) +
      geom_point() +
      geom_line() +
      facet_wrap(~ tau, scales = "free_x",
                  labeller = labeller(tau = label_both)) +
      labs(
        title = "Extremogram for Multiple Tau Values",
        x = "Spatial Lag",
        y = "Extremogram"
      ) +
      theme_minimal() +
      scale_color_manual(values = c("Empirical" = "#a1a1ba", 
                                    "Theoretical" = "#dc7c7c")) +
      theme(legend.position = "bottom")
  } else {
    ggplot(combined_data, aes(x = h_lag, y = variogram, color = type,
                                  linetype = type)) +
      geom_point() +
      geom_line() +
      facet_wrap(~ tau, scales = "free_x",
                labeller = labeller(tau = label_both)) +
      labs(
        title = " ",
        x = "Spatial Lag",
        y = "Variogram"
      ) +
      theme_minimal() +
      scale_color_manual(values = c("Empirical" = "#a1a1ba",
                                    "Theoretical" = "#dc7c7c")) +
      theme(legend.position = "bottom")
  }
}


#' Generate a variogram plot for rpareto
#' 
#' This function generates a variogram plot for the rpareto model.
#' 
#' @param result The result of the variogram estimation.
#' @param df_lags The dataframe of lags.
#' @param wind_df The dataframe of wind data.
#' @param tau_values The values of tau to generate plots for.
#' @param chi A logical value indicating whether to plot the chi values.
#'           Default is FALSE.
#' 
#' @return None, but generates plots.
#' 
#' @import ggplot2
#' 
#' @export
generate_variogram_plots_rpareto <- function(result, lags_r, wind_r) {
  # Get the estimated variogram
  beta1 <- result$beta1
  beta2 <- result$beta2
  alpha1 <- result$alpha1
  alpha2 <- result$alpha2
  eta1 <- result$eta1
  eta2 <- result$eta2

  # Compute advection for each row
  adv_df <- wind_r %>%
    mutate(
      adv_x = (abs(vx)^eta1) * sign(vx) * eta2,
      adv_y = (abs(vy)^eta1) * sign(vy) * eta2
    )

  # Apply theorical_chi for each row's advection
  empirical_chi <- adv_df %>%
    rowwise() %>%
    mutate(chi = list(theoretical_chi(c(beta1, beta2, alpha1, alpha2,
                              adv_x, adv_y), lags_r, latlon = TRUE))) %>%
    unnest(chi)

  # keep chi for hnormv > 0
  empirical_chi <- empirical_chi %>%
    filter(hnormV > 0)

  df_chi <- empirical_chi[, c("hnormV", "vario", "tau")]

  # Plot chi$vario in function on hnorm for each tau
  ggplot(df_chi, aes(x = hnormV, y = vario)) +
    geom_line(color=btfgreen) +
    facet_wrap(~ tau, scales = "free_x",
        labeller = labeller(tau = function(x) paste0("tau = ", x))) +
    labs(
      title = "",
      x = "Spatial Lag",
      y = "Variogram"
    ) +
    theme_minimal() +
    theme(legend.position = "bottom")

}


#' Format Numeric Values into Custom String Representation
#'
#' This function takes a numeric vector and formats each value into a custom
#' string representation. The formatted values are concatenated with
#' underscores.
#'
#' @param x A numeric vector containing the values to be formatted.
#'
#' @return A single string where each formatted value is concatenated with
#' underscores. The formatting rules are as follows:
#' - If the value is negative, the string "neg" is prefixed to the formatted
#'   value.
#' - If the value is an integer, it is formatted as is.
#' - If the value is a decimal:
#'   - For values greater than or equal to 1:
#'     - If there is 1 decimal place, the value is multiplied by 10 and
#'       formatted with 2 digits.
#'     - If there are more than 1 decimal places, the value is multiplied by 100
#'       and formatted with 3 digits.
#'   - For values less than 1:
#'     - If there is 1 decimal place, the value is multiplied by 10 and
#'       formatted with 2 digits.
#'     - If there are more than 1 decimal places, the value is multiplied by 100
#'       and formatted with 3 digits.
#'
#' @examples
#' format_value(c(-1.23, 0.5, 2, -0.01))
#' # Returns: "neg123_05_200_neg001"
#'
#' @export
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

    # Add "neg" prefix if the number was negative
    if (is_negative) {
      formatted_value <- paste0("neg", formatted_value)
    }

    return(formatted_value)
  })

  # Concatenate values
  return(paste(formatted_values, collapse = "_"))
}
