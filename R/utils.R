library(kableExtra)
library(knitr)

# Format coefficient estimates with asterisks
format_estimate <- function(estimate, p_value) {
  asterisks <- ifelse(p_value < 0.001, "***",
                  ifelse(p_value < 0.01, "**",
                    ifelse(p_value < 0.05, "*", "")))
  paste0(estimate, asterisks)
}


save_tablatex <- function(table_latex, filename = "", pval = FALSE,
                    tab_folder = "../thesis/resources/tables/") {
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

