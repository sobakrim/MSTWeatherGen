#' Seasonal Estimation for Space-Time Stochastic Weather Generator
#'
#' Performs seasonal estimation for a space-time stochastic weather generator (MSTSWG). This function processes
#' input weather data to identify weather types for a specified season, estimates scaling parameters, computes
#' transition probabilities between weather types, and estimates parameters for the Gaussian field representing
#' spatial dependencies. It can optionally focus on precipitation as a primary variable for analysis.
#'
#' @param data Array of weather data with dimensions [time, location, variable], where each entry represents
#' a weather observation at a specific time and location.
#' @param dates Vector of dates corresponding to the time dimension of the data, used to filter data within the specified season.
#' @param precipitation Logical flag indicating if precipitation is considered as a primary variable for analysis. Defaults to TRUE.
#' @param names Optionally, names of the variables in the data array to be used for analysis. If not provided
#' and `precipitation` is TRUE, "Precipitation" is assumed to be the first variable, with others following numerically named.
#' If `precipitation` is FALSE, variables are numerically named starting from the first dimension.
#' @param names_weather_types Subset of `names`, specifying the variables to be used for weather type classification.
#' If not provided, all variables specified in `names` are used.
#' @param coordinates Matrix with columns for coordinates of each location, essential for spatial analysis.
#' @param season Vector of integers representing months to define the season for analysis, enabling focused
#' seasonal analysis.
#' @param max_it Maximum number of iterations for optimization procedures within the analysis process.
#' @param tmax Maximum time lag for temporal analysis, relevant for estimating parameters of the Gaussian field.
#' @param dmax Maximum distance for spatial analysis, influencing the estimation of spatial dependencies.
#' @param n1, n2 Parameters defining spatial window size for analysis, crucial for detailed spatial analysis within the Gaussian field model.
#' @return A list containing essential components for weather simulation: filtered dates for the specified season (`dates`), 
#' scaling parameters for each variable (`scale_parm`), identified weather types (`wt`), transition probabilities between weather types (`transitions`),
#' parameters for lambda transformations for each weather type and variable (`lmbd`), and parameters for the Gaussian field model (`gf_par`).
#' @export


MSTSWG_estimation_season = function(data, dates, precipitation = T, names = NULL, 
                                    names_weather_types = NULL, coordinates,
                                 season, max_it, tmax, dmax, n1, n2) {
  # Function to perform seasonal estimation for a space-time stochastic weather generator.
  # It processes input weather data, identifies weather types for a given season,
  # estimates scaling parameters, computes transition probabilities between weather types,
  # and estimates parameters for the Gaussian field representing spatial dependencies.
  #
  # Arguments:
  #   data: Array of weather data with dimensions [time, location, variable].
  #   dates: Vector of dates corresponding to the time dimension of the data.
  #   names: Names of the variables in the data array to be used for analysis.
  #   names_weather_types: Subset of 'names', the variables to be used for weather type classification.
  #   coordinates: Matrix with columns for coordinates of each location.
  #   season: Vector of integers representing months to define the season for analysis.
  #   max_it: Maximum number of iterations for certain optimization procedures.
  #   tmax: Maximum time lag for temporal analysis.
  #   n1, n2: Parameters defining spatial window size for analysis.
  #   return_plots: Logical, indicating if plots should be generated and returned.
  #
  # Returns:
  #   A list containing:
  #     - dates: Filtered dates for the specified season.
  #     - scale_parm: Scaling parameters for each variable.
  #     - wt: Identified weather types for the season.
  #     - transitions: Transition probabilities between weather types.
  #     - lmbd: Parameters for lambda transformations for each weather type and variable.
  #     - gf_par: Parameters for the Gaussian field model.
  
  # Load necessary libraries for statistical operations and date handling
  require(MASS)
  require(lubridate)
  
  # Check names
  if(is.null(names)&precipitation){
    names = c("Precipitation", paste0("v",2:dim(data)[3]))
  }else{
    if(is.null(names)&!precipitation){
      names = paste0("v",1:dim(data)[3])
    }else{
      if(!is.null(names)&precipitation){
        names[1] = "Precipitation"
      }
    }
  }
  if(is.null(names_weather_types)){
    names_weather_types = names
  }
  # Step 1: Filter the input data and dates for the specified season
  filtered = filter_season_data(data, dates, season, names)
  data = filtered$data_filtered
  dates = filtered$dates_filtered
  rm(filtered)
  
  # Step 2: Scale the data
  #scale = scale_data(data, names, dates, window_size = 10)
  #scale_parm = scale$scale_parm  # Scaling parameters for reverting the scaling if needed
  #data = scale$data  # Scaled data
  #rm(scale)
  
  # Step 3-1: Identify weather types for the season
  wt = weather_types(data = data, variables = names_weather_types, dates = dates,coordinates =  coordinates,
                     max_number_wt = 6, return_plots = F)
  wt = wt$cluster # extract weather types 
  
  # Step 3-2: Estimate transition probabilities between weather types
  transitions = estimate_transitions(cluster = wt, dates = dates, nb = 30, K = length(unique(wt)))
  
  # Step 4: Transformations for each variable in each weather type
  lmbd = estimate_lambda_transformations(data = data, wt = wt, names = names, coordinates = coordinates)
  threshold_precip = lmbd$threshold_precip
  lmbd = lmbd$lambda_transformations
  data = transformations(data = data,wt = wt,names = names, coordinates = coordinates,lmbd = lmbd)
  
  # Step 5: Estimate parameters for the Gaussian field model
  gf_par = estimate_gaussian_field_params(data = data, wt = wt, names = names, coordinates = coordinates, 
                                          tmax = tmax, max_it = max_it, n1 = n1, n2 = n2, 
                                          dates = dates, threshold_precip = threshold_precip)
  
  # Return the analysis results
  return(list(dates = dates, wt = wt, 
              transitions = transitions, lmbd = lmbd, gf_par = gf_par))
}

#' Estimate Parameters for Multivariate Space-Time Stochastic Weather Generator
#'
#' This function estimates parameters for a Multivariate Space-Time Stochastic Weather Generator (MSTSWG),
#' allowing for detailed analysis and modeling of weather data. Estimation can be performed on a seasonal basis
#' or annually, depending on the provided data and specified parameters. It handles multiple weather variables,
#' with special consideration for precipitation if indicated.
#'
#' @param data A multi-dimensional array of weather data, encompassing time, location, and various weather variables.
#' @param dates A vector of dates corresponding to the time dimension in the data array, used for temporal analysis.
#' @param by_season Logical flag indicating whether to perform the estimation seasonally (`TRUE`) or annually (`FALSE`).
#' @param seasons A list defining the seasons, each with start and end days and months, required if `by_season` is `TRUE`.
#' @param precipitation Logical, indicating if precipitation should be considered as a primary variable for analysis. Defaults to `TRUE`.
#' @param names Optionally, names of the variables in the data array to be used for analysis. If `precipitation` is `TRUE`
#' and `names` is not provided, "Precipitation" is assumed to be the first variable, with other variables numerically named.
#' @param names_weather_types Specific variables from `names` to be used for classifying weather types. If not provided, it defaults to using all variables specified in `names`.
#' @param coordinates A matrix or data frame containing the geographical coordinates for each location in the data.
#' @param max_it The maximum number of iterations for optimization procedures within the estimation process.
#' @param tmax The maximum temporal lag to be considered in the analysis.
#' @param n1, n2 Parameters defining spatial window sizes for the analysis, influencing the granularity of spatial analysis.
#' @return A list containing the results of the `MSTSWG_estimation_season` function for each season (or for the entire year if `by_season` is `FALSE`),
#' including estimated parameters and other outputs relevant to weather generation, such as weather type classifications and spatial dependencies.
#' @export


MSTSWG_estimation = function(data, dates, by_season = TRUE, seasons, precipitation = T, names = NULL, 
                             names_weather_types = NULL, coordinates,
                             max_it, tmax, n1, n2) {
  # Function for estimating parameters for a Multivariate Space-Time Stochastic Weather Generator (MSTSWG).
  # This can be done either on a seasonal basis or annually, based on the provided data and parameters.
  #
  # Arguments:
  #   data: A multi-dimensional array of weather data, typically including time, location, and various weather variables.
  #   dates: A vector of dates corresponding to the time dimension in the data array.
  #   by_season: Logical flag indicating whether the estimation should be done seasonally (TRUE) or annually (FALSE).
  #   seasons: A list defining the seasons, each with start and end days and months. Necessary if by_season is TRUE.
  #   names: Names of the variables in the data array that will be used in the analysis.
  #   names_weather_types: Specific variables from 'names' that will be used to classify weather types.
  #   coordinates: A matrix or data frame containing the coordinates coordinates for each location in the data.
  #   max_it: The maximum number of iterations for the optimization procedures within the estimation process.
  #   tmax: The maximum temporal lag to be considered in the analysis.
  #   n1, n2: Parameters defining spatial window sizes for the analysis.
  #
  # Returns:
  #   A list of results from the MSTSWG_estimation_season function for each season (or annually, if by_season is FALSE),
  #   including estimated parameters and potentially other analytical outputs relevant to weather generation.
  
  # Validate input parameters, especially the requirement for 'seasons' when 'by_season' is TRUE
  if(by_season & missing(seasons)){
    stop("Please provide a lsit of seasons that cover all the year in following format: seasons <- list(
          s1 = list(min_day = 1, max_day = 29, min_month = 12, max_month = 2),
          s2 = list(min_day = 1, max_day = 31, min_month = 3, max_month = 5),
          s3 = list(min_day = 1, max_day = 31, min_month = 6, max_month = 8),
          s4 = list(min_day = 1, max_day = 30, min_month = 9, max_month = 11)
          )")
  }
  
  # If not analyzing by season, create a default season that spans the entire year
  if (!by_season) {
    seasons <- list(s1 = list(min_day = 1, max_day = 31, min_month = 1, max_month = 12))
  }
  # Check names
  if(is.null(names)&precipitation){
    names = c("Precipitation", paste0("v",2:dim(data)[3]))
  }else{
    if(is.null(names)&!precipitation){
      names = paste0("v",1:dim(data)[3])
    }else{
      if(!is.null(names)&precipitation){
        names[1] = "Precipitation"
      }
    }
  }
  if(is.null(names_weather_types)){
    names_weather_types = names
  }
  # Perform the estimation for each season or for the entire year, based on 'by_season' flag
  swg <- lapply(seasons, function(season) {
    MSTSWG_estimation_season(data = data, dates = dates, precipitation = precipitation, 
                             names = names, names_weather_types = names_weather_types, 
                          coordinates = coordinates, season = season, max_it = max_it, tmax = tmax, 
                          dmax = dmax, n1 = n1, n2 = n2)
  })
  
  return(list(swg = swg, by_season = by_season, names = names, names_weather_types = names_weather_types))
}

#' Simulate Seasonal Weather Data Using MSTSWG
#'
#' Simulates seasonal weather data using the Multivariate Space-Time Stochastic Weather Generator (MSTSWG).
#' This function accounts for specified weather types, their transitions, and employs an autoregressive (AR) model
#' for generating synthetic weather variables across a sequence of dates.
#'
#' @param dates Vector of dates for which to simulate weather data, specifying the period over which weather data is to be simulated.
#' @param first_state Optional initial weather type for starting the simulation. If not provided, the initial state is sampled based on the relative frequencies of different weather types observed in the parameter data.
#' @param Z_initial Optional initial conditions for the AR process. If missing, these are autonomously generated based on the covariance structure of the initial state, ensuring a realistic starting point for the simulation.
#' @param parm Parameters object for the stochastic weather generator, containing essential information such as weather type sequences, transition probabilities, and other relevant parameters for the simulation.
#' @param AR_lag The lag order of the AR model to be utilized in the simulation, defining the dependency of current weather observations on past data.
#' @param bk A list of Bk matrices corresponding to the AR model, one for each weather type. These matrices are pivotal for modeling the temporal evolution of weather variables under different weather conditions.
#' @return A list containing three elements: 'sim', the simulated weather data as a 3D array; 'wt', the sequence of simulated weather types; and 'Z_initial_next', the initial conditions for subsequent simulation steps. This structure provides a comprehensive output for further analysis or continuous simulation processes.
#' @export


MSTSWG_simulation_season = function(dates, first_state = NULL, Z_initial = NULL, parm, AR_lag = 1, bk) {
  # Function to simulate seasonal weather data using the multivariate space-time stochastic weather generator.
  #
  # Arguments:
  #   dates: Vector of dates for which to simulate weather data.
  #   first_state: The initial state (weather type) for the simulation. If missing, it is sampled based on the frequency of states.
  #   Z_initial: Initial conditions for the AR process. If missing, they are generated based on the covariance structure of the initial state.
  #   parm: Parameters of the stochastic weather generator.
  #   AR_lag: The lag order of the AR model to be used in the simulation.
  #   bk: A list of Bk matrices for the AR model, one for each weather type.
  #
  # Returns:
  #   A list containing the simulated weather data ('sim'), the sequence of simulated weather types ('wt'), and the initial conditions for the next simulation ('Z_initial_next').
  
  # Load required libraries
  library(lubridate)
  
  names <- parm$names
  # Determine the number of weather types (K) and their relative frequencies (fr)
  K <- length(unique(parm$wt))
  fr <- sapply(1:K, function(k) length(which(parm$wt == k)) / length(parm$wt))
  parm$wt <- factor(parm$wt, levels = 1:K)
  
  # Estimate transition probabilities if not provided
  #parm$transitions <- estimate_transitions(parm$wt, parm$dates, nb = 20, K)
  
  # Simulate weather types for the given dates
  wt <- if (is.null(first_state)) {
    first_state <- sample(1:K, 1, prob = fr)
    simulate_weathertypes(first_state, dates, parm$dates, parm$transitions, K)
  } else {
    simulate_weathertypes(first_state, dates, parm$dates, parm$transitions, K)
  }
  wt <- as.numeric(wt)
  
  # Generate initial conditions for the AR process if not provided
  if (is.null(Z_initial)) {
    Z_initial <- generate_initial_conditions(AR_lag, bk, wt)
  }
  
  # Simulate the AR process to generate synthetic weather data
  Y <- simulate_Z(bk, AR_lag, length(dates), Z_initial, wt)
  
  # Convert the list of matrices Y to a 3D array 'sim'
  sim <- list_to_array(Y, names, dates)
  
  # Apply inverse transformations to the simulated data

  sim <- apply_inverse_transformations(sim = sim, wt = wt, parm = parm, names = names)

  # Rescale the simulated data for variables other than "Precipitation"
  #sim <- rescale_data(sim, parm, names, dates)
  
  # Prepare initial conditions for the next simulation
  Z_initial_next <- sapply(1:AR_lag, function(m) Y[[length(Y)-m+1]])
  
  return(list(sim = sim, wt = wt, Z_initial_next = Z_initial_next))
}
#' Simulate Weather Data Using MSTSWG
#'
#' This function simulates weather data over specified dates using the Multivariate Space-Time Stochastic Weather Generator (MSTSWG). The simulation can be conducted on either a seasonal basis or for the entire period, depending on the provided parameters and the structure of the historical data. It utilizes autoregressive models and specified parameters to generate realistic weather variables.
#'
#' @param dates_sim Dates for which to simulate weather data, specifying the target period for simulation.
#' @param dates_original The original dates corresponding to the historical weather data, used for aligning and deriving simulation parameters.
#' @param data Historical weather data array, serving as the basis for deriving parameters for simulation models and ensuring realistic weather patterns.
#' @param seasons Optional definitions of seasons, which allows for applying different simulation parameters or models according to season-specific dynamics.
#' @param parm A comprehensive parameters object containing essential data and model parameters for conducting the weather simulation.
#' @param AR_lag The lag order of the autoregressive (AR) model used in the simulation, dictating how previous weather data influences future simulations.
#' @param bk Coefficients matrices for the AR model, which may vary by season or weather type, playing a key role in the temporal modeling of weather variables.
#' @return A 3D array (`sim`) of simulated weather data for the specified dates, encapsulating simulated values for each weather variable across spatial locations and over time.
#' @export

MSTSWG_simulation = function(dates_sim, dates_original, data, seasons = NULL, parm, AR_lag=1, bk) {
  # Function to simulate weather data using the Multivariate Space-Time Stochastic Weather Generator (MSTSWG).
  #
  # Arguments:
  #   dates_sim: Dates for which to simulate weather data.
  #   dates_original: The original dates corresponding to the input data.
  #   data: Historical weather data used for deriving simulation parameters.
  #   seasons: Definitions of seasons to potentially apply different simulation parameters seasonally.
  #   names: Names of weather variables to be simulated.
  #   names_weather_types: Names of weather variables used for weather type classification.
  #   parm: Parameters object containing the parameters.
  #   AR_lag: The lag order of the autoregressive model to be used in the simulation.
  #   bk: Coefficients matrices for the AR model, potentially varying by season or weather type.
  #
  # Returns:
  #   sim: A 3D array of simulated weather data for the specified dates.
  
  by_season = parm$by_season  # Flag indicating whether to simulate data by season.
  parm = parm$swg  # Extract simulation parameters.
  names = parm$names
  names_weather_types = parm$names_weather_types
  if (by_season) {
    # If simulation is to be performed seasonally, assign seasons to the simulation dates.
    seasons_assigned <- assign_seasons(dates_sim, seasons)
    wt_seasons <- lapply(1:length(seasons), function(s) parm[[s]]$wt)
    
    # Identify indices where the season changes.
    change_season_indices <- c(which(diff(seasons_assigned) != 0) + 1, length(seasons_assigned) + 1)
    
    # Find centroids for weather types in each season.
    centroids <- find_centroids(data, dates_original, seasons, wt_seasons, names_weather_types)
    
    # Simulate weather data for the first season.
    sm <- MSTSWG_simulation_season(dates = dates_sim[1:(change_season_indices[1]-1)], names = names, 
                                parm = parm[[seasons_assigned[1]]], AR_lag = AR_lag, bk = bk[[seasons_assigned[1]]])
    sim <- sm$sim
    Z_initial <- sm$Z_initial_next
    
    # Iterate over each season change and simulate weather data.
    for (s in 1:(length(change_season_indices)-1)) {
      # Determine the most probable initial weather type for the next season.
      first_state <- most_probable_weather_type(sim = sim[dim(sim)[1], , ],
                                                centroids = centroids[[seasons_assigned[change_season_indices[s]]]], 
                                                transitions = parm[[seasons_assigned[change_season_indices[s]]]]$transitions, 
                                                names_weather_types = names_weather_types)
      
      # Simulate weather data for the current season segment.
      sm <- MSTSWG_simulation_season(dates = dates_sim[change_season_indices[s]:(change_season_indices[s+1]-1)], 
                                  names = names, Z_initial = Z_initial, first_state = first_state,
                                  parm = parm[[seasons_assigned[change_season_indices[s]]]], 
                                  AR_lag = AR_lag, bk = bk[[seasons_assigned[change_season_indices[s]]]])
      Z_initial <- sm$Z_initial_next
      sim <- abind(sim, sm$sim, along = 1)  # Append the new simulation data to the overall simulation.
    }
  } else {
    # For non-seasonal simulation, simply simulate weather data for the entire period.
    sim <- MSTSWG_simulation_season(dates = dates_sim, names = names, 
                                   parm = parm, AR_lag = AR_lag, bk = bk)$sim
  }
  
  return(sim)
}
