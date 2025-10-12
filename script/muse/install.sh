#!/bin/bash
# Script pour installer un package R proprement

if [ -z "$1" ]; then
  echo "Usage: $0 package_name"
  exit 1
fi

PKG=$1

# Nettoyer l'environnement
module load cmake
module load R
module load gcc/11.3.0
module load gdal/3.4.2
module load geos/3.7.2
module load proj/7.2.1

# Forcer l'utilisation de GCC moderne et C99
export CC=gcc
export CXX=g++
export CFLAGS="-std=gnu99"
export CXXFLAGS="-std=gnu++11"

# Dossier perso pour les packages R
export R_LIBS_USER=~/R/x86_64-pc-linux-gnu-library/4.4
mkdir -p "$R_LIBS_USER"

echo ">>> Installation du package R: $PKG"
Rscript -e "install.packages('$PKG', repos='https://cloud.r-project.org')"

