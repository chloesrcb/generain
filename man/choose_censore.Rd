% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/egpd.R
\name{choose_censore}
\alias{choose_censore}
\title{choose_censore function}
\usage{
choose_censore(rain_df, censores, n_samples = 100)
}
\arguments{
\item{rain_df}{The input dataframe containing the rain data.}

\item{censores}{The censored vector.}

\item{n_samples}{The number of samples to generate. Default is 100.}
}
\value{
The selected censored data in a datafram according to the NRMSE and
        RMSE. # TODO: add CRPS
}
\description{
This function selects the censored data from a rain_df dataframe for each
site according to the NRMSE and RMSE.
}
