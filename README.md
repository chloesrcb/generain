# `generain` package

## Purpose and model

This R package implements a spatio-temporal stochastic model for extreme rainfall using
EGPD marginals and r-Pareto processes, while incorporating advection effects.

The package provides a generic simulation engine for high-resolution spatio-temporal extreme rainfall based on:

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


#### Installation

Before installing `generain`, make sure the dependencies `RandomFieldsUtils` and `RandomFields` are installed. These packages are no longer actively maintained, but fixed versions are proposed by [Ian Flint](https://github.com/iflint1). The fixed versions are available via the GitHub repositories [iflint1/RandomFieldsUtils](https://github.com/iflint1/RandomFieldsUtils) and [iflint1/RandomFields](https://github.com/iflint1/RandomFields).

You can install `generain` directly from GitHub using either `devtools` or `remotes`:

```r
devtools::install_github("chloesrcb/generain")
```

or 

```r
remotes::install_github('chloesrcb/generain', force=TRUE)
```

#### Loading

Once installed, load the package with:

```r
library(generain)
```

## SPG simulations (spatio-temporal extreme rainfall framework)

The framework is implemented in a fully modular way and can be applied to any spatial domain with appropriate inputs.

---

## General simulation engine

The core simulation is implemented in:

`R/spg_episode_simulation.R`

It is domain-independent and can be used for any study area with correct units for each parameter.

### Model inputs

The simulation requires three main components:

#### 1. Marginal model (EGPD and occurrence probability)

The marginal rainfall distribution is defined by:

- `p0`: the occurrence probability;
- `xi`: the tail parameter of the EGPD;
- `sigma`: the scale parameter of the EGPD;
- `kappa`: the bulk parameter of the EGPD.

These parameters can be estimated from observed data and adapted.

#### 2. Extreme spatio-temporal dependence

Dependence is modeled via a non-separable variogram:

$$
\gamma(h, \tau) = \beta_1 ||h - \tau V||^{\alpha_1} + \beta_2 |\tau|^{\alpha_2}
$$

where $h \in \mathbb{R}^2$ is the spatial lag, $\tau \in \mathbb{R}$ is the temporal lag and $V \in \mathbb{R}^2$ is the advection vector.

The parameters $(\beta_1, \beta_2,\alpha_1, \alpha_2)$ can be estimated using the current `generain` package.

#### 3. Advection

Storm displacement is modeled by a constant velocity vector.

For example:
```{r}
adv <- c(100, -200) # here, in m/5min
```

Here, the simulated rainfall pattern moves by 100 m eastward and 200 m northward at each 5-minute time step.

### Simulation settings


The main simulation settings are:

- `u_emp`: the rainfall threshold used in the r-Pareto simulation;
- `nT`: the number of simulated time steps;
- `steps`: the corresponding temporal indices.

For 5-minute rainfall data, `nT = 12` corresponds to a one-hour rainfall episode, , with `steps` ranging from 0 to 11



## Example: Verdanson case study

The script `script/spg/spg_verdanson_simu.R` simulates spatio-temporal rainfall episodes over a regular fine-scale grid covering the Verdanson urban catchment (Montpellier, France).

It includes:

- construction of a regular spatial grid over the catchment,
- simulation of an extreme rainfall episode,
- visualization of spatio-temporal rainfall fields with the advection over a map.


## Generalization to other study areas

The Verdanson example is only a **case study configuration**.

To apply the model to another region, the following elements must be adapted:

- **Geometry file**: replace the geometry file with the target domain
- **Coordinate system**: all spatial data must be in a projected CRS
- **Grid resolution**: adjust `cell_m` depending on spatial scale
- **Advection vector**: ensure consistency with spatial units
- **Model parameters**: recalibrate variogram and EGPD parameters for the new region
- **File paths**: update `data_folder` and output directories (`im_folder`)

Care must be taken to maintain consistency between spatial scale, temporal resolution, and model parameters.