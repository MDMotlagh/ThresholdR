% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fit_k.R
\name{fit_k}
\alias{fit_k}
\title{Fit for all Ks from k_list}
\usage{
fit_k(
  data,
  k_list,
  margin.den = 0.1,
  epsilon = 0.01,
  maxit = 10000,
  maxrestarts = 5,
  seed
)
}
\arguments{
\item{data}{dataframe of normalized ADT values: Cells as rows and Protein as columns}

\item{k_list}{a list of vectors, each vector contains values for possible components fitted best between k=1 to k=k for normalized distribution of each protein}

\item{margin.den}{a margin to evaluate the density of the first component at its mean value (x=mean(C1)) whether it is at most 0.1 times higher than the density of the whole population at x=mean(C1).
if density of C1 is greater than 1.1*density of the whole data at the point of x=mean(C1), a correction step is needed for the standard deviation of the first component (C1), because in some cases that left side of the distribution is truncated, SD is underestimated which leads to overshooting of the density of the first component.}

\item{epsilon}{inherent to normalmixEM function: the convergence criterion.}

\item{maxit}{The maximum number of iterations (Default=10000)}

\item{maxrestarts}{The maximum number of restarts (Default=5)}
}
\value{
Returns a list, `fittings.corrected` including a list for each k fittings with values:

* `x` - The normalized data.
* `lambda` - mixing weights or probabilities.
* `mu` - The mean parameters
* `sigma` - The standard deviation parameters.
* `loglik` - log-likelihood.
* `posterior` - A matrix of posterior probabilities.
* `all.loglik` - A vector of each iteration's log-likelihood
* `restarts` - The number of restarts.
* `ft` - function name, a character vector.
}
\description{
This function fits for all possible Ks extracted through GetAllPossibleFits function.
Fittig is done through normalmixEm function from mixtools R package
}
\examples{
data("adt.df")
fittings.test <- fit_k(data = adt.df[,1:3], k_list = k_list.test, margin.den = 0.1, seed = 42)
}
