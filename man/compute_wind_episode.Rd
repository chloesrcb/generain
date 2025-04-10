% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/optimization.R
\name{compute_wind_episode}
\alias{compute_wind_episode}
\title{compute_wind_episode function}
\usage{
compute_wind_episode(episode, s0, u, wind_df, delta)
}
\arguments{
\item{episode}{The extreme episode dataframe.}

\item{s0}{The starting location name.}

\item{u}{The quantile threshold.}

\item{wind_df}{The wind dataframe.}

\item{delta}{The temporal window size.}
}
\value{
The wind episode dataframe.
}
\description{
Compute the wind speed and direction for a given extreme episode.
}
