#!/bin/bash
# Script pour installer terra 1.6.x et sf sur un cluster
# Compatible avec GDAL 3.4.2, PROJ 7.2.1, GEOS 3.10.6
# ====================================================

# Définir répertoire local pour R
R_LIBS_USER=$HOME/R/library
mkdir -p $R_LIBS_USER

# Charger les modules disponibles sur le cluster
module load R
module load cmake
module load gdal/3.4.2
module load proj/7.2.1
module load geos/3.10.6

# Exporter variables pour compilation R
export LD_LIBRARY_PATH=$HOME/local/lib:$LD_LIBRARY_PATH
export CPPFLAGS="-I$HOME/local/include"
export LDFLAGS="-L$HOME/local/lib"
export PKG_CONFIG_PATH=$HOME/local/lib/pkgconfig:$PKG_CONFIG_PATH

# Installer sf (pré-requis pour terra)
Rscript -e "install.packages('sf', lib='$R_LIBS_USER', repos='https://cloud.r-project.org')"

# Installer terra version 1.6.28 compatible cluster
Rscript -e "install.packages('terra', repos='https://rspatial.r-universe.dev')"

echo "Installation terminée !"
echo "Pour utiliser les packages :"
echo "R_LIBS_USER=$R_LIBS_USER"
echo "Dans R :"
echo ".libPaths(c('$R_LIBS_USER', .libPaths()))"
echo "library(sf)"
echo "library(terra)"

