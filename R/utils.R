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
#' @param type The type of simulation.
#' @param forcedtemp The number of time steps to include in the GIF.
#'
#' @return None, but saves the GIF to a file.
#'
#' @import ggplot2
#' @import reshape2
#' @import animation
#'
#' @export
create_simu_gif <- function(simulation_data, sites_coords, params,
                            type = "rpar", forcedtemp = NA) {
  ngrid <- sqrt(ncol(simulation_data)) # Number of grid points in each dimension
  Tmax <- nrow(simulation_data) # Number of time steps
  if (is.na(forcedtemp)) {
    temp <- 1:Tmax
  } else {
    temp <- 1:forcedtemp
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

  if (length(params) == 6) {
    adv <- params[5:6]
  } else {
    adv <- c(0, 0)
  }

  simulation_data$Time <- rownames(simulation_data) # Add a time column
  simulation_data_long <- melt(simulation_data) # Convert to long format
  simulation_data_long$Time <- as.numeric(simulation_data_long$Time)
  # Create a dataframe to represent grid points
  grid <- expand.grid(x = 1:ngrid, y = 1:ngrid)

  plots <- list()
  cropped_data <- simulation_data_long[simulation_data_long$Time %in% temp, ]
  # for each time step
  for (i in unique(cropped_data$Time)) {
    # Add the simulated values to the grid dataframe
    grid$value <- cropped_data$value[cropped_data$Time == i]

    # Plot
    p <-  ggplot(data = grid, aes(x = x, y = y, fill = value)) +
      geom_tile() +
      scale_fill_gradient(low = "#70a7ae", high = "#9d503d",
                          name = "Rainfall in mm",
                          limits = c(min(cropped_data$value),
                                     max(cropped_data$value))) +
      labs(title = paste0("t =", i, " | Betas: ", beta1, ", ", beta2,
                          " | Alphas: ",
                          alpha1, ", ", alpha2, " | Advection: ", adv[1],
                          ", ", adv[2])) +
      theme_minimal() +
      theme(plot.background = element_rect(fill = "#F9F8F6",
                                         color = "#F9F8F6"),
            panel.border = element_blank(),
            panel.grid = element_blank(),
            axis.text.x = element_blank(),
            axis.text.y = element_blank(),
            axis.title.x = element_blank(),
            axis.title.y = element_blank())

    plots[[i]] <- p
  }

  # Save the plots as a gif
  ani.options(interval = 0.5) # time between frames
  saveGIF({
    for (i in temp) {
      print(plots[[i]])
    }
  }, movie.name = paste0("/user/cserreco/home/Documents/These/generain/images",
                         "/simu_gif/simu_", type, "/", type, "_", ngrid^2, "s_",
                         Tmax, "t.gif"),
  ani.width = 700, ani.height = 600, ani.units = "px", ani.type = "cairo")

}
