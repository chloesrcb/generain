% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/simulations.R
\name{sim_BR}
\alias{sim_BR}
\title{sim_BR function}
\usage{
sim_BR(beta1, beta2, alpha1, alpha2, x, y, t, adv = c(0, 0), nres = 1)
}
\arguments{
\item{beta1}{The value of beta1.}

\item{beta2}{The value of beta2.}

\item{alpha1}{The value of alpha1.}

\item{alpha2}{The value of alpha2.}

\item{x}{Vector for the first dimension (spatial x in our case).}

\item{y}{Vector for the second dimension (spatial y in our case)}

\item{t}{Vector for the third dimension (time in our case).}

\item{adv}{The advection coordinates vector. Default is c(0, 0).}

\item{nres}{The number of simulations to perform. Default is 1.}
}
\value{
The result of the simulation.
}
\description{
This function performs a simulation of a spatio-temporal Brown-Resnick
process using a fractionnal Brownian motion model and based on the
David Leber code with advection.
}
