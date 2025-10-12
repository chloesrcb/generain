# ================================
# R Package Installation Script
# ================================

# List of required packages
packages <- c(
    # General
    "plyr", "dplyr", "tidyr", "data.table",
    
    # Spatial
    "sf", "geodist", "geosphere",
    
    # Temporal
    "datetime",
    
    # Plot
    "reshape2", "Rcpp", "ggplot2", "gridExtra", "latex2exp",
    "ggpubr", "factoextra", "tidyverse", "tibble",
    "hrbrthemes", "reshape", "matrixStats",
    
    # Extremes
    "fExtremes", "extRemes", "ismev", "evd", "mev", "POT",
    
    # Models
    "fields", "lmtest", "bbmle",
    
    # Tables and Latex
    "kableExtra",
    
    # Parallel
    "parallel"
)

# Check and install missing packages
install_if_missing <- function(pkg){
    if(!requireNamespace(pkg, quietly = TRUE)){
        install.packages(pkg, repos = "https://cloud.r-project.org")
    }
}

invisible(lapply(packages, install_if_missing))

# Load all packages
lapply(packages, library, character.only = TRUE)

cat("\nAll packages are installed and loaded.\n")
