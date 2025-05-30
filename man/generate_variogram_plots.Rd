% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{generate_variogram_plots}
\alias{generate_variogram_plots}
\title{Generate variogram plots}
\usage{
generate_variogram_plots(
  result,
  df_lags,
  true_param,
  tau_values,
  chi = FALSE,
  latlon = FALSE,
  directional = FALSE
)
}
\arguments{
\item{result}{The result of the variogram estimation.}

\item{df_lags}{The dataframe of lags.}

\item{true_param}{The true parameters used to generate the data.}

\item{tau_values}{The values of tau to generate plots for.}

\item{chi}{A logical value indicating whether to plot the chi values.
Default is FALSE.}

\item{latlon}{A logical value indicating whether to use latitude and
longitude coordinates. Default is FALSE.}

\item{directional}{A logical value indicating whether to use directional
variogram. Default is FALSE.}
}
\value{
None, but generates plots.
}
\description{
This function generates theoretical and empirical variogram plots for
multiple tau values.
}
