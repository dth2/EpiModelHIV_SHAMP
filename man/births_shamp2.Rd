% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mod.births.shamp2.R
\name{births_shamp2}
\alias{births_shamp2}
\title{Births Module for up to 5 race groups heterosexuals and MSM for preserving demographic distribution.}
\usage{
births_shamp2(dat, at)
}
\arguments{
\item{dat}{Master data list object of class \code{dat} containing networks,
individual-level attributes, and summary statistics.}

\item{at}{Current time step.}
}
\value{
This function updates the \code{attr} list with new attributes for each new
population member, and the \code{nw} objects with new vertices.
}
\description{
Module function for births or entries into the sexually active
             population.
}
\details{
New population members are added based on expected numbers of entries among
all five race/immigrant groups and two sexes, stochastically determined with draws from Poisson
distributions. The proportion of men who are MSM are determined by \code{msm.frac}. 
The proportion of men who are MSMF are determined by \code{msmf.frac}. 
For each new entry, a set of attributes is added for that node,
and the nodes are added onto the network objects. Only attributes that are
a part of the network model formula are updated as vertex attributes on the
network objects.
}
\keyword{het}
\keyword{module}
\keyword{msm}
