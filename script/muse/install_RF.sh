#!/bin/bash

# Exit on any error
set -e

# Check if package name(s) provided (optional)
# PKG="$1"

# Create/append Makevars file for Linux Werror/format-security fix
MAKEVARS_FILE="$HOME/.R/Makevars"
if ! grep -q "Wno-error=format-security" "$MAKEVARS_FILE" 2>/dev/null; then
  echo "CXXFLAGS=-O2 -Wall -Wno-error=format-security" >> "$MAKEVARS_FILE"
  echo "Updated ~/.R/Makevars with format-security fix."
fi

# Uncomment and adjust if needed:
# module purge
module load R
module load gcc/11.3.0

# Install devtools if not installed
Rscript -e "if(!requireNamespace('devtools', quietly=TRUE)) install.packages('devtools', repos='https://cloud.r-project.org')"

# Install RandomFieldsUtils from GitHub
Rscript -e "devtools::install_github('iflint1/RandomFieldsUtils')"

# Install RandomFields from GitHub
Rscript -e "devtools::install_github('iflint1/RandomFields')"

echo "Installation complete! You can now load the library with: library(RandomFields)"
