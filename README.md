# generain package

Model for extreme and moderate precipitation


# Installation

```devtools::install_github("chloesrcb/generain")```

or 

```remotes::install_github('chloesrcb/generain', force=TRUE)```

# Loading

```library(generain)```

# On the cluster

### Needed modules

```
module load R
module load gdal
module load proj
```

### Installation

Previously upload [RandomFields archive](https://cran.r-project.org/src/contrib/Archive/RandomFields/) and [RandomFieldsUtils archive](https://cran.r-project.org/src/contrib/Archive/RandomFieldsUtils/).

```
Rscript -e "install.packages('remotes')"
Rscript -e "install.packages(pkgs = './inst/archives/RandomFields_3.3.14.tar.gz', type='source', repos=NULL)"
Rscript -e "install.packages(pkgs = './inst/archives/RandomFieldsUtils_1.2.5.tar.gz', type='source', repos=NULL)"
Rscript -e "remotes::install_github('chloesrcb/generain', force=TRUE)"
```
