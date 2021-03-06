% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/conattractionmodel.R
\name{settle.con}
\alias{settle.con}
\title{settle.con}
\usage{
settle.con(radius, ID.mat, hab.mat, n.individ, A.coef, ...)
}
\arguments{
\item{radius}{Sets how many cells away from the occupied cells should be included in the "neighborhood".
Passed to `neighbors`.}

\item{ID.mat}{A matrix of the same size as the simulated landscape but with the value of each cell corresponding
to it's cell ID. Passed to `neighbors`.}

\item{hab.mat}{The matrix containing the "grown" habitats enumerated by the strength of preference
(e.g., the output from pref.strength, not convert.cell.) Note that which hab.mat to use (i.e., which
preference strength to use) is arbitrary for the null model, but the choice must match the `A.coef` supplied (below).}

\item{n.individ}{The number of individuals included in the ABM.}

\item{A.coef}{The habitat preference parameter that matches the hab.mat supplied.  This is the "key" that
allows the function to interpret each cell of the matrix as either "Habitat A" or "Habitat B".}

\item{...}{Additional arguments as necessary.}
}
\value{
The proportion of simulated settled locations located in Habitat A as well as the
raw locations for each simulation run
}
\description{
Agents settle landscape based on conspecific attraction
This function acts iteratively to allow agents to settle the landscape selecting randomly from cells defined as
available based on proximity to previously settled cells.  The first agent settles the landscape completely
randomly. Instanciation as a function allows the entire ABM to be run an arbitray number of times (see 'reps' below)
 in order to generate a sampling distribution of apparent selections generated by the model.
}
