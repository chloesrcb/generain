% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/optimization.R
\name{theorical_chi_ind}
\alias{theorical_chi_ind}
\title{Calculate the theoretical chi value.}
\usage{
theorical_chi_ind(params, h, tau)
}
\arguments{
\item{params}{The variogram parameter values (beta1, beta2, alpha1, alpha2).}

\item{h}{The spatial lag value.}

\item{tau}{The temporal lag value.}
}
\value{
The theoretical chi value.
}
\description{
This function calculates the individual theoretical spatio-temporal chi value
based on a the variogram parameters, a spatial lag and a temporal lag.
}
