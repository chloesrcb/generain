# `generain` package

This package implements a spatio-temporal stochastic model for extreme rainfall using
EGPD marginals and r-Pareto processes, while incorporating advection effects.

<p align="center">
  <img src="https://raw.githubusercontent.com/chloesrcb/extreme-rainfall-generator/main/images/simulations/simulated_episode_1.gif" width="300"><br>
  <em>Simulated extreme rainfall episode using the proposed spatio-temporal model.</em>
</p>


## Installation

```r
# install.packages("devtools") # if not installed
devtools::install_github("chloesrcb/generain")
```

or 

```r
# install.packages("remotes") # if not installed
remotes::install_github('chloesrcb/generain', force=TRUE)
```

## Loading

```r
library(generain)
```

## On the cluster

If you are working on a cluster, load the required modules before installing.

### Needed modules

```
module load R
module load gdal
module load proj
```

### Installation

```
 Rscript -e "remotes::install_github('chloesrcb/generain', force=TRUE)"
```

```
 Rscript -e 'install.packages("reticulate", repos="https://cloud.r-project.org")'
```

### Installation on R versions < 4.4

Previously upload [RandomFields archive](https://cran.r-project.org/src/contrib/Archive/RandomFields/) and [RandomFieldsUtils archive](https://cran.r-project.org/src/contrib/Archive/RandomFieldsUtils/).

```
Rscript -e "install.packages('remotes',repos = 'https://cran.r-project.org')"
Rscript -e "install.packages(pkgs = './inst/archives/RandomFields_3.3.14.tar.gz', type='source', repos=NULL)"
Rscript -e "install.packages(pkgs = './inst/archives/RandomFieldsUtils_1.2.5.tar.gz', type='source', repos=NULL)"
Rscript -e "remotes::install_github('chloesrcb/generain', force=TRUE)"
```

## Installation issue

### RandomFields Installation Guide for R 4.4.3

The `RandomFieldsUtils` and `RandomFields` packages are no longer actively maintained, but fixed versions are proposed by [Ian Flint](https://github.com/iflint1). The fixed versions are available via the GitHub repositories [iflint1/RandomFieldsUtils](https://github.com/iflint1/RandomFieldsUtils) and [iflint1/RandomFields](https://github.com/iflint1/RandomFields).


#### Clean Previous Installations

Before installing, make sure you've removed any old or broken installations of `RandomFields`.

Check your R library paths by running in R:

```r
.libPaths()
```

Then manually delete any existing `RandomFields` folder from those directories if present, for example:

```
rm -rf ~/R/x86_64-pc-linux-gnu-library/4.4/RandomFields
```

#### Step 1: Install the Packages in R

Open R and run:

```r
# install.packages("devtools") # if not installed
devtools::install_github("iflint1/RandomFieldsUtils")
devtools::install_github("iflint1/RandomFields")
```

If this completes without error, you're good and can upload the library with `library(RandomFields)`.

#### Step 2: Fix for `Werror=format-security` Error (Linux only)

##### 1. Create the `.R` configuration directory

In your terminal:

```
mkdir ~/.R
```

##### 2. Create or edit the `Makevars` file

Open the file with a text editor:

```
nano ~/.R/Makevars
```

Then add the following line and save the file:

```make
CXXFLAGS=-O2 -Wall -Wno-error=format-security
```

##### 3. Retry Installation

Go back into R and rerun the installation (Step 1).

