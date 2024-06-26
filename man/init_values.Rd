% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/egpd.R
\name{init_values}
\alias{init_values}
\title{Function to initialize values for EGPD fitting}
\usage{
init_values(y, u)
}
\arguments{
\item{y}{The data to be fitted.}

\item{u}{The threshold value.}
}
\value{
GPD estimates for sigma and xi.
}
\description{
This function is used to initialize values for the EGPD fitting with
the estimates of the classical GPD parameters.
It performs the necessary setup and initialization steps given a classical
GPD distribution with a threshold equals to u. For EGPD, classical threshold
will be u=0.
}
