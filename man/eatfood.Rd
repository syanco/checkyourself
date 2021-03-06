% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/eatfood.R
\name{eatfood}
\alias{eatfood}
\title{eatfood}
\usage{
eatfood(num_individ, steps, num_sources, diet_prop, food)
}
\arguments{
\item{num_individ}{Integer, the number of individuals to simulate}

\item{steps}{Integer, the number of time steps to simulate. Passed to
`eatfoodsteps`.}

\item{num_sources}{integer, the number of food sources to choose from.
Can not exceed the actual number of available food sources. Passed to
`eatfoodsteps` and then to `eatfoodonce`.}

\item{diet_prop}{vector, the relative probability of selecting a source.
Must have an entry for each source. Passed to `eatfoodsteps` and then to
`eatfoodonce`.}

\item{food}{list, the simulated food sources from which to select. Passed to
`eatfoodsteps` and then to `eatfoodonce`.}
}
\value{
A list, each element of which represents an individual consumer and
 contains the time step matrix supplied by `eatfoodsteps`.
}
\description{
Wrapper function to iterate `eatfoodsteps` across multiple individuals.
}
