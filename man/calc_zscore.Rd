% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/calc.R
\name{calc_zscore}
\alias{calc_zscore}
\title{Calculate a pseudo-z-score of a difference distribution}
\usage{
calc_zscore(diff, disease.code)
}
\arguments{
\item{diff}{A "long" dataframe containing Metabs and diff values}

\item{disease.code}{The 3-letter disease code for the current disease}
}
\value{
A zscore dataframe containing zscores for each metabolite
}
\description{
This function uses a calculated diff (calc_diff) to calculate a pseudo-z-score for the
distribution.
}
