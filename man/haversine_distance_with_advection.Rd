% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/distances.R
\name{haversine_distance_with_advection}
\alias{haversine_distance_with_advection}
\title{haversine_distance_with_advection function}
\usage{
haversine_distance_with_advection(lat1, lon1, lat2, lon2, adv, tau)
}
\arguments{
\item{lat1}{Latitude of the first point.}

\item{lon1}{Longitude of the first point.}

\item{lat2}{Latitude of the second point.}

\item{lon2}{Longitude of the second point.}

\item{adv}{A vector of advection values (in m/s).}

\item{tau}{Temporal lag (in seconds).}
}
\value{
The Haversine distance between the two points with advection applied
and the direction of the vector in radians.
}
\description{
This function calculates the Haversine distance between two points
with advection applied and gives the direction of the vector in radians
with meteorological orientation (0 degrees is North).
}
