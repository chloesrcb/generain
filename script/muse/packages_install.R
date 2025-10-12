
# ================================

# Liste des packages nécessaires
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

# Vérifie et installe les packages manquants
install_if_missing <- function(pkg){
  if(!requireNamespace(pkg, quietly = TRUE)){
    install.packages(pkg, repos = "https://cloud.r-project.org")
  }
}

invisible(lapply(packages, install_if_missing))

# Charge tous les packages
lapply(packages, library, character.only = TRUE)

cat("\n✅ Tous les packages sont installés et chargés.\n")

