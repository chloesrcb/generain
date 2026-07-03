# `generain` package

## Purpose and model

This R package implements a spatio-temporal stochastic model for extreme rainfall using
EGPD marginals and r-Pareto processes, while incorporating advection effects.

The proposed framework combines:

- **EGPD marginals** for extreme rainfall intensities,
- **r-Pareto processes** to model extreme spatio-temporal dependence structure,
- **advection effects** to account for rain storm displacement


This package can be used to simulate high-resolution spatio-temporal extreme rainfall episodes over the Verdanson water catchment in Montpellier, France. This application and the methodology is described in:

> **Spatio-temporal modeling of urban extreme rainfall events at high resolution**  
> Serre-Combe, Meyer, Opitz, and Toulemonde

Available at [arXiv.2602.19774](https://arxiv.org/abs/2602.19774)

<p align="center">
  <img src="https://raw.githubusercontent.com/chloesrcb/extreme-rainfall-generator/main/images/simulations/output.gif" width="360"><br>
  <em>Simulated extreme rainfall episode over the Verdanson catchment using the proposed model.</em>
</p>


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


## SPG simulations on the Verdanson water catchment

The script `script/spg/spg_verdanson_simu.R` simulates spatio-temporal rainfall episodes over a regular fine-scale grid covering the Verdanson catchment.

Simulations are based on a spatial r-Pareto process combined with an extended generalized Pareto distribution (EGPD) for the marginal rainfall distribution. The spatial and temporal dependence structure is controlled by fitted variogram parameters and a chosen advection vector.

### Generalization to other study areas

Although this script is currently configured for the Verdanson catchment, it can be adapted to any other spatial domain. This requires:

- replacing the catchment geometry file (`verdanson_basin.geojson`) by the target domain geometry;
- ensuring that all spatial inputs (rain gauges, grids, and coordinates) are expressed in a consistent projected coordinate system (e.g. metres);
- updating the grid resolution (`cell_m`) according to the desired spatial scale;
- checking that the advection vector `adv_m5` is expressed in the same spatial units per time step as the coordinate system and variogram parameters;
- updating file paths (`data_folder`, `im_folder`, and input data locations) accordingly.

Care must be taken to maintain consistency between spatial units, temporal resolution, and model parameters, particularly the variogram and advection components.

### Require inputs

Before running the script, make sure that the following objects and files are available:

- `data_folder` and `im_folder` defined in the project configuration or in `script/load_libraries.R`;
- `omsev/loc_rain_gauges.csv`, containing the longitude, latitude, and station name of the rain gauges (optional);
- `geometry/verdanson_basin.geojson`, defining the Verdanson catchment boundary;
- the `generain` package, or the local source files `spg.R` and `distances.R`.

### 1. Build the simulation grid

A regular square grid is created over the Verdanson catchment. The grid resolution is controlled by:

```{r}
cell_m <- 100 # in meters
```

Only grid cells whose centroids fall inside the catchment are retained. Each cell is assigned a unique identifier of the form `pixel_1`, `pixel_2`, etc. The conditioning location is selected through `s0_pixel_id`. This pixel should be chosen among the grid cells located inside the catchment.

### 2. Define the simulation settings

The main simulation settings are:

- `u_emp`: the rainfall threshold used in the r-Pareto simulation;
- `nT`: the number of simulated time steps;
- `steps`: the corresponding temporal indices.

For 5-minute rainfall data, `nT = 12` corresponds to a one-hour rainfall episode, , with `steps` ranging from 0 to 11

### 3. Specify marginal parameters

The marginal distribution is described by common EGPD parameters with zero point mass in `params_margins_common` with:

- `p0`: the occurrence probability;
- `xi`: the tail parameter of the EGPD;
- `sigma`: the scale parameter of the EGPD;
- `kappa`: the bulk parameter of the EGPD.

They can be replaced by estimated parameters from another dataset, time resolution, or study area.

### 4. Specify dependence and advection parameters

The spatio-temporal dependence model is defined through variogram parameters $(\beta_1, \beta_2,\alpha_1, \alpha_2)$, stored in `params_m5min` at the desired resolution scale (here, meters per 5 minutes).
The variogram is defined as:
\[
\gamma(h, \tau) = \beta_1 ||h - \tau V||^{\alpha_1} + \beta_2 |\tau|^{\alpha_2}
\]
where $h \in \mathbb{R}^2$ is the spatial lag, $\tau \in \mathbb{R}$ is the temporal lag and $V \in \mathbb{R}^2$ is the advection vector.

The rainfall field is advected using a constant displacement vector, specified in metres per time step through `adv_m5`.
For example:
```{r}
adv_m5 <- c(100, 200)
```

Here, the simulated rainfall pattern moves by 100 m eastward and 200 m northward at each 5-minute time step.

### 5. Run a simulation

A single rainfall episode is generated with:

```{r}
sim_episode <- sim_episode_grid_m5(
  params_vario = params_m5min,
  params_margins_common = params_margins_common,
  coords = grid_df_l93,
  steps = steps,
  adv_m5 = adv_5m,
  t0 = 0,
  s0_pixel_id = s0_pixel_id,
  u_emp = u_emp,
  seed = 2026,
  cell_m = 100
)
```

The output `sim_episode` is a matrix where:

- rows correspond to grid pixels;
- columns correspond to time steps;
- entries contain simulated rainfall intensities in mm per 5 minutes.

The row names match the grid-cell identifiers used in the spatial grid.
