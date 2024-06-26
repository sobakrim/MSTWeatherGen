% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Simulation_helper_functions.R
\name{list_to_array}
\alias{list_to_array}
\title{Convert List to 3D Array}
\usage{
list_to_array(Y, names, dates)
}
\arguments{
\item{Y}{A list of matrices, where each matrix represents simulated data for a specific time step.}

\item{names}{Vector of names representing the variables in the simulated data (columns of the matrices).}

\item{dates}{Vector of dates representing the time steps.}
}
\value{
A 3D array where the first dimension corresponds to time steps, the second dimension corresponds to spatial locations, and the third dimension corresponds to variables.
}
\description{
Converts a list of matrices into a 3D array.
}
\keyword{internal}
