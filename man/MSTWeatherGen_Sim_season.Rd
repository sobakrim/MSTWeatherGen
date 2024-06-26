% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/MSTWeatherGen.R
\name{MSTWeatherGen_Sim_season}
\alias{MSTWeatherGen_Sim_season}
\title{Seasonal Simulation in MSTWeatherGen}
\usage{
MSTWeatherGen_Sim_season(
  dates,
  names,
  first_state = NULL,
  Z_initial = NULL,
  parm,
  AR_lag = 1,
  bk
)
}
\arguments{
\item{dates}{Vector of dates for which to simulate weather data.}

\item{names}{Names of weather variables to be simulated.}

\item{first_state}{Optional initial state (weather type) for the simulation. If not provided, it is determined based on state frequencies.}

\item{Z_initial}{Optional initial conditions for the autoregressive (AR) process. If not provided, conditions are generated based on the initial state's covariance structure.}

\item{parm}{Parameters object containing the stochastic weather generator settings, including weather types, transition probabilities, and scaling parameters.}

\item{AR_lag}{Lag order of the AR model to be used in the simulation.}

\item{bk}{A list of Bk matrices for the AR model, one for each weather type, essential for generating correlated weather variables across time.}
}
\value{
A list containing:
\itemize{
  \item \code{sim}: A 3D array of simulated weather data for the specified dates.
  \item \code{wt}: A vector indicating the sequence of simulated weather types.
  \item \code{Z_initial_next}: Initial conditions for subsequent simulations, facilitating sequential simulation processes.
}
}
\description{
Simulates seasonal weather data using the Multivariate Space-Time Stochastic Weather Generator (MSTWeatherGen).
This function is designed to generate synthetic weather data for a specified season, leveraging an autoregressive model 
and predefined parameters to accurately reflect weather type transitions and spatial-temporal correlations.
}
\keyword{internal}
