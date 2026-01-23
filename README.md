# `generain` package

## Purpose and model

This R package implements a spatio-temporal stochastic model for extreme rainfall using
EGPD marginals and r-Pareto processes, while incorporating advection effects.

<p align="center">
  <img src="https://raw.githubusercontent.com/chloesrcb/extreme-rainfall-generator/main/images/simulations/simulated_episode_1.gif" width="300"><br>
  <em>Simulated extreme rainfall episode using the proposed spatio-temporal model.</em>
</p>

The proposed framework combines:
- **EGPD marginals** for extreme rainfall intensities,
- **r-Pareto processes** to model extreme spatio-temporal dependence structure,
- **advection effects** to account for rain storm displacement

## Installation

Before installing `generain`, make sure the dependencies `RandomFieldsUtils` and `RandomFields` are installed. These packages are no longer actively maintained, but fixed versions are proposed by [Ian Flint](https://github.com/iflint1). The fixed versions are available via the GitHub repositories [iflint1/RandomFieldsUtils](https://github.com/iflint1/RandomFieldsUtils) and [iflint1/RandomFields](https://github.com/iflint1/RandomFields).

You can install `generain` directly from GitHub using either `devtools` or `remotes`:

```r
devtools::install_github("chloesrcb/generain")
```

or 

```r
remotes::install_github('chloesrcb/generain', force=TRUE)
```

## Loading

Once installed, load the package with:

```r
library(generain)
```
