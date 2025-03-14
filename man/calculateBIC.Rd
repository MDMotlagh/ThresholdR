% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/calculateBIC.R
\name{calculateBIC}
\alias{calculateBIC}
\title{Claculate BIC values}
\usage{
calculateBIC(data, k = 3)
}
\arguments{
\item{data}{A data frame of CLR normalized ADT expression values. Cells are in the rows and ADT probes (surface markers) are in columns.}

\item{k}{The maximimum expected number of underlying components across all markers. Default is set to 3.}
}
\value{
A list of numerical vectors of BIC values corresponding to the markers. If BIC value calculation for a marker encounters an error, the function returns NA for that marker in the list.
}
\description{
Claculate BIC values
}
\examples{
bic_values <- calculateBIC(data = adt, k = 3)
}
