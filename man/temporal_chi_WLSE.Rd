% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dependence.R
\name{temporal_chi_WLSE}
\alias{temporal_chi_WLSE}
\title{Calculate temporal chi using WLSE method}
\usage{
temporal_chi_WLSE(dftemp, weights)
}
\arguments{
\item{dftemp}{The input data frame containing the temporal chi values.}

\item{weights}{The weights to be used in the calculation: "residuals", "exp",
or "none"}
}
\value{
The calculated temporal chi value.
}
\description{
This function calculates the temporal chi using the Weighted Least Squares
Estimation (WLSE) method.
}
\examples{
data <- data.frame(time = c(1, 2, 3), value = c(10, 20, 30))
weights <- c(0.5, 0.3, 0.2)
result <- temporal_chi_WLSE(data, weights)
print(result)

}
