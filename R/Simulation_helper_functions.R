#' Calculate Bk Matrices for Autoregressive (AR) Model
#'
#' Computes the coefficient matrices (Bk) for an autoregressive (AR) model using lagged covariance matrices.
#'
#' @param C_k_matrices A list of covariance matrices, where each matrix represents the covariance at a different time lag. 
#'                     C_k_matrices[[1]] is the covariance matrix at lag 0.
#'
#' @return A list of Bk matrices (coefficient matrices) for the AR model, including Bk_0, which is derived separately from the rest of the Bk matrices.
#'
#' @importFrom Matrix chol2inv
#'
#' @keywords internal
calculate_Bk_matrices <- function(C_k_matrices, Bk_0) {
  # Function for calculating coefficients matrices (Bk) for an autoregressive (AR) model
  # using lagged covariance matrices (C_k_matrices).
  #
  # Arguments:
  #   C_k_matrices: A list of covariance matrices, where each matrix represents the covariance
  #                 at a different time lag. C_k_matrices[[1]] is the covariance matrix at lag 0.
  #
  # Returns:
  #   A list of Bk matrices (coefficient matrices) for the AR model, including Bk_0,
  #   which is derived separately from the rest of the Bk matrices.
  
  M <- length(C_k_matrices) - 1  # Determine AR maximum lag
  
  matrix_size <- nrow(C_k_matrices[[1]])  
  
  # Step 1: Construct the matrix A as a block matrix composed of lagged covariance matrices
  A_blocks <- list()
  for (i in 1:M) {
    row_blocks <- list()
    for (j in 1:M) {
      row_blocks[[j]] <- C_k_matrices[[abs(i - j) + 1]]  # Use absolute difference in lag indices to find corresponding C_k matrix
    }
    A_blocks[[i]] <- do.call(cbind, row_blocks)  # Combine row blocks horizontally
  }
  A <- do.call(rbind, A_blocks)  # Combine all row blocks vertically to form the matrix A
  
  # Step 2: Construct the matrix v as a vertical concatenation of covariance matrices from C_k_matrices, excluding the first one
  v <- do.call(rbind, C_k_matrices[2:(M + 1)])
  
  # Step 3: Solve for Bk matrices (Bk_1 to Bk_M) 
  Bk_matrices_solution <- t(v) %*% Matrix::chol2inv(chol(A)) 
  
  # Step 4: Reshape the solution vector into individual Bk matrices
  Bk_matrices_list <- list()
  for (i in 1:M) {
    start_col <- (i - 1) * matrix_size + 1
    end_col <- i * matrix_size
    Bk_matrices_list[[i]] <- Bk_matrices_solution[, start_col:end_col]
  }
  
  # Step 5: Calculate Bk_0 by adjusting the lag 0 covariance matrix with contributions from higher-order lags
  Bk_0_rhs <- C_k_matrices[[1]]  # Start with the lag 0 covariance matrix
  for (i in 1:M) {
    Bk_0_rhs <- Bk_0_rhs - Bk_matrices_list[[i]] %*% C_k_matrices[[i + 1]]  # Adjust for contributions from Bk_i matrices
  }
  Bk_0 <- try(t(chol(Bk_0_rhs)), silent = T) # Perform Cholesky decomposition to obtain Bk_0
  if(is.character(Bk_0)) Bk_0 = Bk_0
  Bk_matrices_list <- list(Bk_0 = Bk_0, bk = Bk_matrices_list)
  
  return(Bk_matrices_list)
}
#' Calculate AR Coefficient Matrices
#'
#' Computes autoregressive (AR) coefficient matrices for each season and weather state based on the parameters of a fitted spatial weather model and the spatial coordinates of observation locations. This function is part of the process for preparing data for spatial weather simulations, allowing the AR process to consider spatial correlations.
#'
#' @param parm A list containing the fitted parameters of the spatial weather model. This list should include elements for each season, and within each season, parameters (`gf_par`) for each weather state.
#' @param coordinates A matrix or data frame of geographic coordinates for the locations represented in the model. The coordinates are used to calculate spatial covariance matrices.
#' @param AR_lag The lag order of the autoregressive model to be used in calculating the coefficient matrices. Specifies how many past time points should be considered in the AR process.
#'
#' @return A list of AR coefficient matrices (`bk`) for each season and weather state, along with the base covariance matrix (`cov0`) used in their calculation. The structure of the return list mirrors that of the `parm` input, with a set of coefficient matrices for each season and weather state.
#'
#' @examples
#' \dontrun{
#'   # Assuming `fitted_parms` contains the fitted model parameters and
#'   # `site_coordinates` contains the geographic coordinates for your locations:
#'   AR_coeffs = calculate_AR_coefficients_matrices(parm = fitted_parms,
#'                                                  coordinates = site_coordinates,
#'                                                  AR_lag = 1)
#' }
#'
#' @export
calculate_AR_coefficients_matrices <- function(parm, coordinates, AR_lag){
  names = parm$names
  bk = lapply(1:length(parm$swg), function(s){
    K = length(parm$swg[[s]]$gf_par)
    bk = lapply(1:K, function(k){
      cova = cov_matrices(par = parm$swg[[s]]$gf_par[[k]], coordinates = coordinates, 
                          names = names, M = AR_lag)
      bk = try(calculate_Bk_matrices(cova), silent = T)
      return(list(bk = bk, cov0 = cova[[1]]))
    })
    return(bk)
  })
  for(s in 1:length(parm$swg)){
    K = length(parm$swg[[s]]$gf_par)
    for(k in 1:K){
      j = k-1
      while (is.character(bk[[s]][[k]]$bk$Bk_0)) {
        bk[[s]][[k]]$bk$Bk_0 = try(bk[[s]][[j]]$bk$Bk_0, silent = T)
        bk[[s]][[k]]$bk$bk = list(try(bk[[s]][[j]]$bk$Bk, silent = T))
        bk[[s]][[k]]$cov0 = try(bk[[s]][[k]]$cov0, silent = T)
        j = k+1
      }
    }
  }
  return(bk)
}
#' Simulate AR Process
#'
#' Simulates the autoregressive (AR) process to generate synthetic weather data.
#'
#' @param Bk A list of Bk matrices for the AR model.
#' @param M The lag order of the AR model.
#' @param num_steps The number of time steps to simulate.
#' @param Z_initial Initial conditions for the AR process.
#' @param wt A vector of weather types corresponding to each time step.
#'
#' @return A list containing the simulated values of Z for each time step.
#'
#' @keywords internal
simulate_Z <- function(Bk, M, num_steps, Z_initial,wt) {
  n <- nrow(Bk[[1]]$bk$Bk_0)  # Assuming Bk0 is square and represents the dimension of Z
  Z_list <- vector("list", num_steps)  # List to store the simulated values of Z
  
  # Fill in initial values for Z
  for (t in 1:min(M, num_steps)) {
    Z_list[[t]] <- Z_initial[,t]
  }
  
  # Simulating the process
  for (t in (M + 1):num_steps) {
    
    BkZ_product <- matrix(0, n, 1)
    for (i in 1:M) {
      BkZ_product <- BkZ_product + Bk[[wt[t]]]$bk$bk[[i]] %*%Z_list[[t-i]]
    }
    
    Z_list[[t]] <- Bk[[wt[t]]]$bk$Bk_0 %*%  rnorm(n) + BkZ_product
  }
  
  return(Z_list)
}
#' Generate Initial Conditions for AR Process
#'
#' Generates initial conditions for the autoregressive (AR) process based on Bk matrices and initial weather type.
#'
#' @param AR_lag The lag order of the AR model.
#' @param bk A list of Bk matrices for the AR model, one for each weather type.
#' @param wt The initial weather type (state) for the simulation.
#'
#' @return A matrix representing the initial conditions for the AR process simulation.
#'
#' @keywords internal
generate_initial_conditions <- function(AR_lag, bk, wt) {
  # Generate initial conditions for the AR process based on Bk matrices and initial weather type.
  #
  # Arguments:
  #   AR_lag: The lag order of the AR model.
  #   bk: A list of Bk matrices for the AR model, one for each weather type.
  #   wt: The initial weather type (state) for the simulation.
  #
  # Returns:
  #   Z_initial: A matrix representing the initial conditions for the AR process simulation.
  
  Z_initial <- sapply(1:AR_lag, function(m) {
    cov0 <- bk[[wt[m]]]$cov0
    L <- chol(cov0)
    rnorm_vals <- rnorm(ncol(cov0))
    Zm <- L %*% rnorm_vals
    return(Zm)
  })
  num_steps <- 20
  Z_initial <- simulate_Z(Bk = bk, M = AR_lag, num_steps = num_steps, Z_initial = Z_initial, wt = wt)
  Z_initial <- as.matrix(do.call(cbind, Z_initial)[,(length(Z_initial)-AR_lag+1):length(Z_initial)])
  return(Z_initial)
}
#' Convert List to 3D Array
#'
#' Converts a list of matrices into a 3D array.
#'
#' @param Y A list of matrices, where each matrix represents simulated data for a specific time step.
#' @param names Vector of names representing the variables in the simulated data (columns of the matrices).
#' @param dates Vector of dates representing the time steps.
#'
#' @return A 3D array where the first dimension corresponds to time steps, the second dimension corresponds to spatial locations, and the third dimension corresponds to variables.
#'
#' @keywords internal
list_to_array <- function(Y, names, dates) {
  # Converts a list of matrices into a 3D array.
  #
  # Arguments:
  #   Y: A list of matrices, where each matrix represents simulated data for a specific time step.
  #   names: Vector of names representing the variables in the simulated data (columns of the matrices).
  #
  # Returns:
  #   A 3D array where the first dimension corresponds to time steps, the second dimension corresponds to
  #   spatial locations, and the third dimension corresponds to variables.
  
  sim = lapply(1:length(Y), function(i) matrix(Y[[i]], ncol = length(names)))
  sim = array(unlist(sim), dim = c(nrow(sim[[1]]), ncol(sim[[1]]),length(sim)))  
  sim = aperm(sim , c(3,1,2))
  
  sim_array = array(sim, dim = dim(sim) ,
                    dimnames = list(dates,
                                    1:dim(sim)[2], names))
  
  return(sim_array)
}
#' Predict Most Probable Weather Type
#'
#' Predicts the most likely weather type for the next time step in a simulated weather series.
#'
#' @param sim A matrix containing simulated weather data for a single time step.
#' @param centroids A list of centroids, where each centroid is a vector representing the mean values of weather variables for a specific weather type.
#' @param transitions A list (or matrix) of transition probabilities between weather types. Each element transitions[[i]][j, k] represents the probability of transitioning from weather type j to k.
#' @param names_weather_types Names of the weather variables used to determine the weather type.
#'
#' @return The predicted weather type for the next time step, selected based on the calculated probabilities.
#'
#' @keywords internal
most_probable_weather_type = function(sim, centroids, transitions, names_weather_types) {
  # Function to predict the most likely weather type for the next time step in a simulated weather series.
  #
  # Arguments:
  #   sim: A matrix containing simulated weather data for a single time step.
  #   centroids: A list of centroids, where each centroid is a vector representing the mean values of weather variables for a specific weather type.
  #   transitions: A list (or matrix) of transition probabilities between weather types. 
  #                Each element transitions[[i]][j, k] represents the probability of transitioning from weather type j to k.
  #   names_weather_types: Names of the weather variables used to determine the weather type.
  #
  # Returns:
  #   new_wt: The predicted weather type for the next time step, selected based on the calculated probabilities.
  
  # Determine the number of weather types from the length of the centroids list
  K = length(centroids)
  
  # Find the previous weather type by identifying the centroid closest to the current simulated weather data
  previous_wt = which.min(sapply(1:K, function(k) {
    mean((c(sim[, names_weather_types]) - centroids[[k]])^2)
  }))
  
  # Predict the next weather type by sampling from the distribution defined by transition probabilities from the previous weather type
  new_wt = sample(1:K, 1, prob = transitions[[1]][previous_wt, ])
  
  return(new_wt)
}
#' Find Centroids
#'
#' Calculates centroids of weather types for each season.
#'
#' @param data A multi-dimensional array of weather data, with dimensions [time, locations, variables].
#' @param dates A vector of dates corresponding to the time dimension in the data array.
#' @param seasons A list defining the start and end dates of each season.
#' @param wt_seasons A list containing vectors of weather type classifications for each season.
#' @param names_weather_types Names of the weather variables used to determine the weather type centroids.
#'
#' @return A list of lists, where each inner list contains the centroids for each weather type within a specific season.
#'
#' @keywords internal
find_centroids = function(data, dates, seasons, wt_seasons, names_weather_types) {
  # Calculates centroids of weather types for each season.
  #
  # Arguments:
  #   data: A multi-dimensional array of weather data, with dimensions [time, locations, variables].
  #   dates: A vector of dates corresponding to the time dimension in the data array.
  #   seasons: A list defining the start and end dates of each season.
  #   wt_seasons: A list containing vectors of weather type classifications for each season.
  #   names_weather_types: Names of the weather variables used to determine the weather type centroids.
  #
  # Returns:
  #   centroids: A list of lists, where each inner list contains the centroids for each weather type within a specific season.
  
  centroids = lapply(1:length(seasons), function(s) {
    # Identify the indices for the current season
    season_indices = season_indices(dates, seasons[[s]])
    # Subset the data for the current season and relevant weather variables
    season_data = data[season_indices, , names_weather_types]
    
    # Determine the number of unique weather types in the current season
    K = length(unique(wt_seasons[[s]]))
    
    # Calculate centroids for each weather type within the season
    season_centroids = lapply(1:K, function(k) {
      # Calculate the mean of each weather variable for the current weather type
      centroid = colMeans(season_data[wt_seasons[[s]] == k, , ])
      return(centroid)
    })
    
    return(season_centroids)
  })
  
  return(centroids)
}
#' Assign Seasons
#'
#' Assigns a season to each date based on predefined season boundaries.
#'
#' @param dates A vector of dates to classify into seasons.
#' @param seasons A list where each element represents a season with its start and end dates and months.
#'
#' @return A numeric vector where each element represents the season assigned to the corresponding date in the input vector.
#'
#' @keywords internal
assign_seasons <- function(dates, seasons) {
  # Assigns a season to each date based on predefined season boundaries.
  #
  # Arguments:
  #   dates: A vector of dates to classify into seasons.
  #   seasons: A list where each element represents a season with its start and end dates and months.
  #
  # Returns:
  #   A numeric vector where each element represents the season assigned to the corresponding date in the input vector.
  
  # Initialize a vector to store the assigned season for each date
  seasons_assigned <- numeric(length(dates))
  
  # Loop through each date to assign a season
  for (i in seq_along(dates)) {
    # Extract the month and day from the current date
    month <- as.numeric(format(dates[i], "%m"))
    day <- as.numeric(format(dates[i], "%d"))
    
    # Determine which season the date falls into
    for (s in names(seasons)) {
      season <- seasons[[s]]
      # Handle seasons within a single year
      if (season$max_month >= season$min_month) {
        if ((month > season$min_month && month < season$max_month) ||
            (month == season$min_month && day >= season$min_day) ||
            (month == season$max_month && day <= season$max_day)) {
          seasons_assigned[i] <- as.numeric(substr(s, 2, 2))
          break
        }
      } else {
        # Handle seasons that cross the year boundary
        if ((month > season$min_month || month < season$max_month) ||
            (month == season$min_month && day >= season$min_day) ||
            (month == season$max_month && day <= season$max_day)) {
          seasons_assigned[i] <- as.numeric(substr(s, 2, 2))
          break
        }
      }
    }
  }
  
  return(seasons_assigned)
}