
#' Estimate Parameters for Multivariate Space-Time Stochastic Weather Generator
#'
#' This function estimates parameters for a Multivariate Space-Time Stochastic Weather Generator (MSTWeatherGen),
#' allowing for detailed analysis and modeling of weather data. Estimation can be performed on a seasonal basis
#' or annually, depending on the provided data and specified parameters. It handles multiple weather variables,
#' with special consideration for precipitation if indicated.
#'
#' @param data A multi-dimensional array of weather data, encompassing time, location, and various weather variables.
#' @param dates A vector of dates corresponding to the time dimension in the data array, used for temporal analysis.
#' @param by_season Logical flag indicating whether to perform the estimation seasonally (`TRUE`) or annually (`FALSE`).
#' @param seasons A list defining the seasons, each with start and end days and months, required if `by_season` is `TRUE`.
#' @param scale Logical, indicating if the data needs to be standardized (`TRUE`) or not (`FALSE`). If `scale` is `TRUE` each meteorological variable
#' is standardized at each location using smoothed mean and standard deviation. 
#' @param precipitation Logical, indicating if precipitation should be considered as a primary variable for analysis. Defaults to `TRUE`.
#' @param names Optionally, names of the variables in the data array to be used for analysis. If `precipitation` is `TRUE`
#' and `names` is not provided, "Precipitation" is assumed to be the first variable, with other variables numerically named.
#' @param names_weather_types Specific variables from `names` to be used for classifying weather types. If not provided, it defaults to using all variables specified in `names`.
#' @param coordinates A matrix or data frame containing the geographical coordinates for each location in the data.
#' @param max_it The maximum number of iterations for optimization procedures within the estimation process.
#' @param tmax The maximum temporal lag to be considered in the analysis.
#' @param n1 First parameter defining spatial window size for analysis, crucial for detailed spatial analysis within the Gaussian field model.
#' @param n2 Second parameter defining spatial window size.
#' @return A list containing the results of the `MSTWeatherGen_Estim_season` function for each season (or for the entire year if `by_season` is `FALSE`),
#' including estimated parameters and other outputs relevant to weather generation, such as weather type classifications and spatial dependencies.
#' @details
#' The model considers a \(p\)-dimensional stochastic process \(Y(s,t) = [Y_i(s,t)]_{i=1}^{p}\), defined in \(R^2 x R\), 
#' representing \(p\) meteorological variables at space-time coordinates \((s,t)\). 
#' Additionally, a \(1\)-dimensional stochastic process \(X(t)\) characterizes the weather state at time \(t\), 
#' taking values in a discrete state space \(S = \{1,...,K\}\). We assume the existence of \(K\) independent 
#' \(p\)-dimensional latent Gaussian random fields \(Z_k(s,t)=[Z_{k,i}(s,t)]_{i=1}^{p}\), with \(k=1,...,K\), 
#' such that when \(X(t)\) is in state \(k\), the stochastic vector \(Y(s,t)\) is related to \(Z_k(s,t)\) 
#' through a non-linear transformation function \(Psi_{k,i,s}\).
#'
#' Each component of the random field \(Z_k(s,t)\) is assumed to have zero mean, and the field is 
#' considered second-order stationary, with its covariance function depending only on the space-time lag \((h,u)\).
#' @export


MSTWeatherGen_Estim = function(data, dates, by_season = TRUE,  seasons, scale = FALSE, precipitation = T, names = NULL, 
                               names_weather_types = NULL, coordinates,
                               max_it, tmax, n1, n2){
  # Function for estimating parameters for a Multivariate Space-Time Stochastic Weather Generator (MSTWeatherGen).
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
  #   A list of results from the MSTWeatherGen_Estim_season function for each season (or annually, if by_season is FALSE),
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
    MSTWeatherGen_Estim_season(data = data, dates = dates, scale = scale, precipitation = precipitation, 
                               names = names, names_weather_types = names_weather_types, 
                               coordinates = coordinates, season = season, max_it = max_it, tmax = tmax, 
                               n1 = n1, n2 = n2)
  })
  
  return(list(swg = swg, by_season = by_season, names = names, names_weather_types = names_weather_types))
}

#' Simulate Weather Data Using MSTWeatherGen
#'
#' This function simulates weather data over specified dates using the Multivariate Space-Time Stochastic Weather Generator (MSTWeatherGen). The simulation can be conducted on either a seasonal basis or for the entire period, depending on the provided parameters and the structure of the historical data. It utilizes autoregressive models and specified parameters to generate realistic weather variables.
#'
#' @param dates_sim Dates for which to simulate weather data, specifying the target period for simulation.
#' @param dates_original The original dates corresponding to the historical weather data, used for aligning and deriving simulation parameters.
#' @param data Historical weather data array, serving as the basis for deriving parameters for simulation models and ensuring realistic weather patterns.
#' @param seasons Optional definitions of seasons, which allows for applying different simulation parameters or models according to season-specific dynamics.
#' @param parm A comprehensive parameters object containing essential data and model parameters for conducting the weather simulation.
#' @param AR_lag The lag order of the autoregressive (AR) model used in the simulation, dictating how previous weather data influences future simulations.
#' @param bk Coefficients matrices for the AR model, which may vary by season or weather type, playing a key role in the temporal modeling of weather variables.
#' @return A 3D array (`sim`) of simulated weather data for the specified dates, encapsulating simulated values for each weather variable across spatial locations and over time.
#' @importFrom abind abind
#' @export

MSTWeatherGen_Sim = function(dates_sim, dates_original, data, seasons = NULL, parm, AR_lag=1, bk) {
  # Function to simulate weather data using the Multivariate Space-Time Stochastic Weather Generator (MSTWeatherGen).
  #
  # Arguments:
  #   dates_sim: Dates for which to simulate weather data.
  #   dates_original: The original dates corresponding to the input data.
  #   data: Historical weather data used for deriving simulation parameters.
  #   seasons: Definitions of seasons to potentially apply different simulation parameters seasonally.
  #   parm: Parameters object containing the parameters.
  #   AR_lag: The lag order of the autoregressive model to be used in the simulation.
  #   bk: Coefficients matrices for the AR model, potentially varying by season or weather type.
  #
  # Returns:
  #   sim: A 3D array of simulated weather data for the specified dates.
  
  by_season = parm$by_season # Flag indicating whether to simulate data by season.
  names = parm$names
  names_weather_types = parm$names_weather_types
  parm = parm$swg  # Extract simulation parameters.
  if (by_season) {
    # If simulation is to be performed seasonally, assign seasons to the simulation dates.
    seasons_assigned <- assign_seasons(dates_sim, seasons)
    wt_seasons <- lapply(1:length(seasons), function(s) parm[[s]]$wt)
    
    # Identify indices where the season changes.
    change_season_indices <- c(which(diff(seasons_assigned) != 0) + 1, length(seasons_assigned) + 1)
    
    # Find centroids for weather types in each season.
    centroids <- find_centroids(data, dates_original, seasons, wt_seasons, names_weather_types)
    
    # Simulate weather data for the first season.
    sm <- MSTWeatherGen_Sim_season(dates = dates_sim[1:(change_season_indices[1]-1)], names = names, 
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
      sm <- MSTWeatherGen_Sim_season(dates = dates_sim[change_season_indices[s]:(change_season_indices[s+1]-1)], 
                                     Z_initial = Z_initial, first_state = first_state, names = names,
                                     parm = parm[[seasons_assigned[change_season_indices[s]]]], 
                                     AR_lag = AR_lag, bk = bk[[seasons_assigned[change_season_indices[s]]]])
      Z_initial <- sm$Z_initial_next
      sim <- abind::abind(sim, sm$sim, along = 1)  # Append the new simulation data to the overall simulation.
    }
  } else {
    # For non-seasonal simulation, simply simulate weather data for the entire period.
    sim <- MSTWeatherGen_Sim_season(dates = dates_sim, names = names, 
                                    parm = parm, AR_lag = AR_lag, bk = bk)$sim
  }
  
  return(sim)
}
#' Seasonal Estimation for MSTWeatherGen
#'
#' Performs seasonal estimation for the MSTWeatherGen package, which includes identification of weather types,
#' estimation of scaling parameters, computation of transition probabilities between weather types, and estimation
#' of parameters for the Gaussian field model.
#'
#' @param data Array of weather data with dimensions [time, location, variable].
#' @param dates Vector of dates corresponding to the time dimension of the data.
#' @param precipitation Logical indicating if precipitation is to be considered as a primary variable.
#' @param scale Logical, indicating if the data needs to be standardized (`TRUE`) or not (`FALSE`). If `scale` is `TRUE` each meteorological variable
#' is standardized at each location using smoothed mean and standard deviation. 
#' @param names Names of the variables in the data array to be used for analysis.
#' @param names_weather_types Subset of 'names', variables to be used for weather type classification.
#' @param coordinates Matrix with columns for the coordinates of each location.
#' @param season Vector of integers representing months to define the season for analysis.
#' @param max_it Maximum number of iterations for optimization procedures.
#' @param tmax Maximum time lag for temporal analysis.
#' @param n1, n2 Parameters defining spatial window size for analysis.
#'
#' @return A list containing dates, scale parameters, weather types, transition probabilities, lambda transformations,
#' and parameters for the Gaussian field model for the specified season.
#'
#' @keywords internal
MSTWeatherGen_Estim_season = function(data, dates, precipitation = T, scale = FALSE, names = NULL, 
                                      names_weather_types = NULL, coordinates,
                                      season, max_it, tmax, n1, n2) {
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
  
  if(scale){
    # Step 2: Scale the data
    scale = scale_data(data, names, dates, window_size = 10)
    scale_parm = scale$scale_parm  # Scaling parameters for reverting the scaling if needed
    data = scale$data  # Scaled data
    rm(scale)
  }else{
    scale_parm = NULL 
  }
  
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
  return(list(dates = dates, wt = wt, scale_parm = scale_parm,
              transitions = transitions, lmbd = lmbd, gf_par = gf_par))
}
#' Seasonal Simulation in MSTWeatherGen
#'
#' Simulates seasonal weather data using the Multivariate Space-Time Stochastic Weather Generator (MSTWeatherGen).
#' This function is designed to generate synthetic weather data for a specified season, leveraging an autoregressive model 
#' and predefined parameters to accurately reflect weather type transitions and spatial-temporal correlations.
#'
#' @param dates Vector of dates for which to simulate weather data.
#' @param names Names of weather variables to be simulated.
#' @param first_state Optional initial state (weather type) for the simulation. If not provided, it is determined based on state frequencies.
#' @param Z_initial Optional initial conditions for the autoregressive (AR) process. If not provided, conditions are generated based on the initial state's covariance structure.
#' @param parm Parameters object containing the stochastic weather generator settings, including weather types, transition probabilities, and scaling parameters.
#' @param AR_lag Lag order of the AR model to be used in the simulation.
#' @param bk A list of Bk matrices for the AR model, one for each weather type, essential for generating correlated weather variables across time.
#'
#' @return A list containing:
#' \itemize{
#'   \item \code{sim}: A 3D array of simulated weather data for the specified dates.
#'   \item \code{wt}: A vector indicating the sequence of simulated weather types.
#'   \item \code{Z_initial_next}: Initial conditions for subsequent simulations, facilitating sequential simulation processes.
#' }
#' @keywords internal

MSTWeatherGen_Sim_season = function(dates, names, first_state = NULL, Z_initial = NULL, parm, AR_lag = 1, bk) {
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
  
  
  if(!is.null(parm$scale_parm)){
    # Rescale the simulated data for variables other than "Precipitation" if required
    sim <- rescale_data(sim, parm, names, dates)
  }
  
  # Prepare initial conditions for the next simulation
  Z_initial_next <- sapply(1:AR_lag, function(m) Y[[length(Y)-m+1]])
  
  return(list(sim = sim, wt = wt, Z_initial_next = Z_initial_next))
}
