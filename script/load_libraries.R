# Package
library(generain)

# General
library(plyr)
library(dplyr)
library(tidyr)
library(reshape2)
library(data.table)

# spatial
# library(raster)
# library(rasterVis)
# library(elevatr)
# library(rgeos)
# library(leaflet)
# library(RColorBrewer)
# library(rgl)
library(sf)
# library(ggspatial)
# library(stars)
library(geodist)
library(geosphere)
# library(igraph)

# temporal
library(datetime)

# Plot
library(reshape2)
library(Rcpp)
library(ggplot2)
library(gridExtra)
library(latex2exp)
library(ggpubr)
library(factoextra)
library(tidyverse)
library(tibble)
library(hrbrthemes)
library(reshape)
library(matrixStats)
library(gridExtra)

# Extremes
library(fExtremes)
library(extRemes)
library(ismev)
library(evd)
library(mev)
library(POT)

# Models
library(fields)
library(lmtest)
library(bbmle)

# Tables and Latex
library(kableExtra)
library(data.table)

# Parallel
library(parallel)

# install.packages(c("reshape2", "Rcpp",
#                     "gridExtra", "latex2exp", "ggpubr", "factoextra",
#                     "tidyverse", "tibble", "hrbrthemes", "reshape",
#                     "matrixStats", "fExtremes", "extRemes", "ismev", "evd",
#                     "mev", "POT", "fields", "lmtest", "bbmle", "parallel"))


# PERSONNALIZATION ############################################################
# ggplot theme
# personnalized theme
btf_theme <- theme_minimal() +
  theme(axis.text.x = element_text(size = 15, angle = 0),
        axis.text.y = element_text(size = 15),
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20),
        axis.line = element_blank(),  # Remove axis lines
        panel.border = element_blank(),  # Remove plot border
        panel.background = element_rect(fill = "transparent", color = NA),
        # Remove plot background
        axis.text = element_blank(),  # Remove axis text
        axis.ticks = element_blank(),  # Remove axis ticks
        plot.background = element_rect(fill = "transparent", color = NA),
        panel.grid = element_line(color = "#5c595943"))

btf_boxplot_theme <- btf_theme +
  theme(
    axis.line = element_blank(),  # Remove axis lines
    panel.border = element_blank(),  # Remove plot border
    panel.background = element_rect(fill = "transparent",
                                    color = NA),  # Remove plot background
    axis.text = element_blank(),  # Remove axis text
    axis.text.x = element_text(angle = 0),
    axis.ticks = element_blank(),  # Remove axis ticks
    panel.grid = element_line(color = "#5c595943")
  )
# my green color
btfgreen <- "#69b3a2"

# images folder
im_folder <- "../phd_extremes/thesis/resources/images/"
data_folder <- "../phd_extremes/data/"
