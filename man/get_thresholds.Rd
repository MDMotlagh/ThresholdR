% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_thresholds.R
\name{get_thresholds}
\alias{get_thresholds}
\title{Calculate different thresholds for bi- and tri-modal fittings}
\usage{
get_thresholds(fittings, k_list)
}
\arguments{
\item{fittings}{A list of fitted parameters for each distribution (component) under each k across all markers.}

\item{k_list}{A list with set of possible k values for each marker as an element.}
}
\value{
A data frame with calculated thresholds: bi_thr1 defines the mean+3SD of the first component and bi_cut refers to the cutpoint on x axis where the two components meet in a bimodal (k=2) distribution. tri_thr1 and tri_thr2 indicate mean+3SD of the first and second components in a trimodal (k=3) distribution. Subsequently, tri_cut1 and tri_cut2 represent the cutpoints between first and second component, and second and third components, repectively, in a trimodal distribution. The user should choose which thereshold to use the recomended value is mean+3SD of the first component in k=2.
}
\description{
Calculate different thresholds for bi- and tri-modal fittings
}
\examples{
thresholds <- get_thresholds(fittings = fittings, k_list = k.list)
}
