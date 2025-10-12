library(xtable)

# Load libraries and set theme
source("./script/load_libraries.R")

# Get all files in the folder "R"
functions_folder <- "./R"
files <- list.files(functions_folder, full.names = TRUE)
# # load all functions in files
invisible(lapply(files, function(f) source(f, echo = FALSE)))
library(latex2exp)

foldername <- paste0(data_folder,
                    "/comephore/optim_results/lalpha/free_eta/")

filenames <- list.files(foldername, full.names = TRUE, pattern = "\\.csv$")
# remove combined_optim_results.csv if exists
filenames <- filenames[!grepl("combined_optim_results.csv", filenames)]

all_results <- list()

for(file in filenames){
    data <- read.csv(file)
    
    row <- data[1, ]
    file_base <- basename(file)
    parts <- unlist(strsplit(file_base, "_"))
    q_value <- as.numeric(sub("q", "", parts[2]))
    delta_value <- as.numeric(sub("delta", "", parts[3]))
    dmin_value <- as.numeric(sub("dmin", "", sub("\\.csv$", "", parts[4])))
    row$q <- q_value
    row$delta <- delta_value
    row$dmin <- dmin_value
    all_results[[length(all_results) + 1]] <- row
}

# reorder columns to have q, delta, dmin at the beginning
final_df <- do.call(rbind, all_results)
final_df <- final_df[, c("q", "delta", "dmin", "beta1", "beta2", "alpha1", "alpha2", "eta1", "eta2")]
final_df <- final_df[order(final_df$q, final_df$delta, final_df$dmin), ]

# save df as csv
write.csv(final_df, file = paste0(foldername, "combined_optim_results.csv"),
          row.names = FALSE)

# colnames in latex format
sanitize_names <- function(x){
  x <- gsub("beta1", "$\\\\widehat{\\\\beta}_1$", x)
  x <- gsub("beta2", "$\\\\widehat{\\\\beta}_2$", x)
  x <- gsub("alpha1", "$\\\\widehat{\\\\alpha}_1$", x)
  x <- gsub("alpha2", "$\\\\widehat{\\\\alpha}_2$", x)
  x <- gsub("eta1", "$\\\\widehat{\\\\eta}_1$", x)
  x <- gsub("eta2", "$\\\\widehat{\\\\eta}_2$", x)
  x <- gsub("q", "$q$", x)
  x <- gsub("delta", "$\\\\delta$", x)
  x <- gsub("dmin", "$d_{\\\\min}$", x)
  return(x)
}

folder_tab <- "../phd_extremes/thesis/resources/tables/"
output_file <- paste0(folder_tab, "comephore_optim_results.tex")
ncol <- ncol(final_df)
align_string <- paste0("c", paste(rep("c", ncol), collapse=""))
latex_table <- xtable(final_df, align = align_string)

print(latex_table,
  include.rownames = FALSE,
  sanitize.colnames.function = sanitize_names,
  floating = FALSE,
  file = output_file,
  booktabs = TRUE)


foldername <- paste0(data_folder,
                    "/comephore/optim_results/jackknife_estimates/")

filenames <- list.files(foldername, full.names = TRUE, pattern = "\\.csv$")
# keep only those that begin with "results_jk" not "all_results_jk"
filenames <- filenames[grepl("results_jk_by_seasonyear", filenames) & !grepl("all_results_jk", filenames)]

all_results <- list()

file <- filenames[1]
file_base <- basename(file)
parts <- unlist(strsplit(file_base, "_"))
q_value <- as.numeric(sub("q", "", parts[6]))
delta_value <- as.numeric(sub("delta", "", parts[7]))
dmin_value <- as.numeric(sub("dmin", "", sub("\\.csv$", "", parts[8])))


library(dplyr)
library(knitr)
library(kableExtra)

# Read the data
data <- read.csv(file)

# Format CI as [lower, upper]
data <- data %>%
  mutate(
    JK = sprintf("$%.3f ~ [%.3f, %.3f]$", Estimate_jk, CI_lower, CI_upper)  # note the ~ before [
  ) %>%
  select(Parameter, Estimate_full, JK)

# Replace parameter names with LaTeX equivalents
data$Parameter <- sanitize_names(data$Parameter)

# round estimates to 3 decimal places
data$Estimate_full <- round(data$Estimate_full, 3)
# Make the LaTeX table
caption_text <- paste0("Estimates of variogram parameters and Jackknife estimates with 95\\% CI for $q = ", q_value / 100,
                       "$, $\\delta = ", delta_value,
                       "$ hours, $d_{\\min} = ", dmin_value, "$ km.",
                       " Jackknife estimates are computed by removing one season-year at a time.")
results <- kable(
  data,
  format = "latex",
  booktabs = TRUE,
  col.names = c("Parameter", "Estimate", "Jackknife estimate"),
  escape = FALSE,
  align = "c",
  caption = caption_text,
  linesep = ""   # <- empêche l'insertion de \addlinespace
) %>%
  kable_styling(latex_options = c("H", "scale_down"))

# colnames(results) <- sanitize_names(colnames(results))

folder_tab <- "../phd_extremes/thesis/resources/tables/"
filename <- paste0("comephore_optim_results_jk_ci_seasons_q", q_value, "_delta", delta_value, "_dmin", dmin_value, ".tex")
output_file <- paste0(folder_tab, filename)

# Save the LaTeX table to a file
writeLines(results, con = output_file)

