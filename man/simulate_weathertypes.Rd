% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Weather_types.R
\name{simulate_weathertypes}
\alias{simulate_weathertypes}
\title{Simulate Weather Types}
\usage{
simulate_weathertypes(first_state, dates_sim, dates, tm, K)
}
\arguments{
\item{first_state}{The initial weather type from which to start the simulation.}

\item{dates_sim}{A vector of dates for which the weather types are to be simulated.}

\item{dates}{A vector of dates corresponding to the transition matrices in 'tm'.}

\item{tm}{A list of transition matrices, each representing the transition probabilities between weather types for a specific date.}

\item{K}{The total number of unique weather types or states.}
}
\value{
A vector of simulated weather types or states for each date in 'dates_sim'.
}
\description{
Function to simulate weather types over specified dates using transition matrices.
}
\keyword{internal}
