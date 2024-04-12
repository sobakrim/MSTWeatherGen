#' Generate All Combinations of Variable Index Pairs
#'
#' Constructs a matrix containing all possible combinations of variable names, including self-pairs and cross-pairs. This facilitates the analysis and modeling of interactions between different variables in a multivariate dataset.
#'
#' @param names Vector of variable names from which pairs are to be generated.
#'
#' @return A data frame where each row represents a pair of variables, with columns 'v1' and 'v2' indicating the variable names in the pair. This includes both self-pairs (where 'v1' and 'v2' are the same) and cross-pairs (where 'v1' and 'v2' are different).
#'
#' @keywords internal
#' @noRd

generate_variable_index_pairs <- function(names) {
  # This function creates a matrix with all combinations of variable names
  u1 = sapply(names, function(v1) sapply(names, function(v2) v1))
  u2 = sapply(names, function(v1) sapply(names, function(v2) v2))
  ep = data.frame(cbind(u1[!upper.tri(u1)],u2[!upper.tri(u2)]))
  names(ep) = c("v1","v2")
  ep = rbind(ep[ep$v1==ep$v2,], ep[!ep$v1==ep$v2,])
  
  return(ep)
}
#' Initialize Model Parameters If Missing
#'
#' Sets up the `par_all` vector with default values or based on provided parameters if it hasn't been initialized. This function ensures that all necessary model parameters are prepared for the modeling process.
#'
#' @param par_all Existing vector of all model parameters; if NULL, it will be initialized.
#' @param names Vector of variable names involved in the model.
#' @param pairs Generated pairs of variables for which parameters are set.
#' @param par_s Initial scaling parameters for the covariance function.
#' @param ax Correction term parameters to be updated in `par_all`.
#' @param cr Initial correlation matrix used for beta computation.
#'
#' @return Updated `par_all` vector with all model parameters, including default and specified values.
#'
#' @keywords internal
#' @noRd

initialize_par_all_if_missing <- function(par_all, names, pairs, par_s, ax, cr) {
  # Initialize the `par_all` vector if it is missing, with default values or using `par_s`
  if (is.null(par_all)) {
    names_par_all <- c(paste(pairs, "dij", sep = ":"), "a1", "d1", "g1", "a2", "d2", "g2",
                       "b1", "e1", "l1", "b2", "e2", "l2", "c", "f", "m",
                       paste(names, "ai", sep = ":"), paste(names, "bi", sep = ":"),
                       paste(names, "ci", sep = ":"),
                       paste(pairs, "rij", sep = ":"), paste(pairs, "vij", sep = ":"), 
                       paste(pairs, "ax", sep = ":"))
    
    par_all <- setNames(rep(0.1, length(names_par_all)), names_par_all)
    
    par_all[paste(pairs, "dij", sep = ":")] = 1
    par_all[paste(pairs[1:length(names)], "rij", sep = ":")] <- par_s[1,] 
    par_all[paste(pairs[1:length(names)], "vij", sep = ":")] <- par_s[2,] 
    par_all[paste(pairs, "ax", sep = ":")] <- 0
    parms <- c("a1", "a2", "d1", "d2", "g1", "g2")
    par_all[parms] <-rep(1, length(parms))
  }
  
  # Update ax parameters based on covariance information
  par_all <- update_ax_parameters(par_all, names, ax)
  
  parm <- param(par_all, names)
  beta <- try(compute_beta(parm, names, cr), silent = T)
  ch <- try(chol(beta), silent = T)
  if(is.character(ch)){
    par_s <- matrix(rep(1, length(names)^2), ncol = length(names), nrow= length(names))
    par_all[paste(pairs[1:length(names)], "rij", sep = ":")] <- 1
    par_all[paste(pairs[1:length(names)], "vij", sep = ":")] <- 1
    par_all <- update_ax_parameters(par_all, names, ax)
  } 
  return(par_all)
}
#' Update Ax Parameters in Model Parameters
#'
#' Modifies the 'ax' parameters within the complete set of model parameters (`par_all`) using the covariance information provided by the 'ax' matrix. This adjustment is crucial for ensuring accurate covariance structures in the model.
#'
#' @param par_all The complete set of model parameters, including 'ax' values to be updated.
#' @param names Vector of variable names, indicating the variables for which 'ax' adjustments are applied.
#' @param ax Matrix or data frame containing the updated covariance information to adjust 'ax' parameters in `par_all`. If `ax` is not a matrix, it will be transformed to ensure positive definiteness before updating.
#'
#' @return The modified `par_all` vector with updated 'ax' parameters reflecting the provided covariance information.
#'
#' @keywords internal
#' @noRd
#' @importFrom Matrix nearPD
update_ax_parameters <- function(par_all, names, ax) {
  # Update the `ax` parameters in `par_all` based on the covariance information in `ax`
  if(!is.matrix(ax)){
    for (v1 in names) {
      for(v2 in names){
        par_all[paste(paste(v1,v2,sep = "-"), "ax", sep = ":")] = ax$cov[ax$v1==v1&ax$v2==v2|ax$v2==v1&ax$v1==v2]
      }
    }
    a = sapply(names, function(v1){
      sapply(names, function(v2){
        ax = par_all[paste(paste(v1,v2,sep = "-"), "ax", sep = ":")]
        if(is.na(ax)) ax = par_all[paste(paste(v2,v1,sep = "-"), "ax", sep = ":")]
        return(ax)
      })
    })
    ax = Matrix::nearPD(a)$mat
  }else{
    ax = Matrix::nearPD(ax)$mat
  }
  colnames(ax) = rownames(ax) = names
  for (v1 in names) {
    for(v2 in names){
      par_all[paste(paste(v1,v2,sep = "-"), "ax", sep = ":")] = ax[v1,v2]
    }
  }
  return(par_all)
}
#' Initialize Spatial Parameters for Variables
#'
#' Optimizes initial spatial parameters for each variable in a dataset. This function employs
#' a log-likelihood based optimization approach to determine initial values of spatial parameters
#' that best fit the spatial structure of the data for each variable.
#'
#' @param data The dataset for which spatial parameters are being initialized, typically
#'        a 3D array or a list of spatial observations for multiple variables.
#' @param names A vector of variable names for which spatial parameters are to be initialized.
#' @param h Vector of spatial distances used in the likelihood function, representing
#'        the spatial relationship between observation points.
#' @param uh Matrix containing additional data required by the likelihood function, 
#'        typically combining spatial and temporal distances with other relevant information.
#' @param max_it Maximum number of iterations for the optimization process, defaulting to 2000.
#'        This parameter controls the depth of the optimization search.
#'
#' @return A list containing optimized spatial parameters for each variable specified in 'names'.
#'         Each element of the list corresponds to a set of parameters (e.g., range and smoothness)
#'         optimized for the spatial structure of the corresponding variable.
#' @keywords internal
#' @importFrom parallel mclapply
init_space_par <- function(data, names, h, uh, max_it = 2000) {
  # Initializes spatial parameters for each variable
  # by optimizing an initial log-likelihood function (loglik0)
  
  # Arguments:
  #   data: The dataset for which spatial parameters are being initialized. This could be
  #         a matrix or data frame of observed values across locations.
  #   names: A vector of variable names for which spatial parameters are to be initialized.
  #   D: A matrix or data frame containing distances.
  #   h: Vector of spatial distances used in the likelihood function.
  #   uh: Matrix containing combined spatial and temporal distances, along with additional data,
  #       used in the likelihood function.
  #   mpar: A list or vector of model parameters that are held fixed during the optimization process.
  #   max_it: Maximum number of iterations for the optimization process. Default is 2000.
  
  # Returns:
  #   A list of optimized parameters for each variable. Each element in the list corresponds to
  #   the set of parameters optimized for one of the variables specified in 'names'.
  
  # Perform parallel optimization for each variable using mclapply
  par = parallel::mclapply(names, function(v) {
    optim(
      par = c(1, 1),  # Initial parameter guesses
      fn = loglik_spatial,     # Objective function to minimize (negative log-likelihood)
      data = data,
      v = v,            # Current variable being optimized
      h = h,
      uh = uh,
      control = list(maxit = max_it, trace = 2)  # Optimization control settings
    )$par  # Extract the optimized parameters
  })
  
  return(par)
}
#' Optimize Spatial Parameters for Variable Pairs
#'
#' Performs optimization of model parameters for each pair of variables based on spatial data.
#' This function iteratively optimizes spatial parameters to maximize the log-likelihood
#' of observing the data given the model, focusing on the spatial interaction between pairs of variables.
#'
#' @param par_all Initial set of all model parameters before optimization.
#' @param data The dataset used for optimization, typically a 3D array or a list
#'        of spatial observations for multiple variables over time.
#' @param names Vector of variable names involved in the optimization process.
#' @param Vi Matrix indicating the pairs of variables for which parameters are optimized.
#' @param uh Matrix containing additional data required by the likelihood function,
#'        combining spatial and temporal distances with other relevant information.
#' @param cr Initial correlation matrix used in the optimization process to provide
#'        a starting point for parameter estimation.
#' @param max_it Maximum number of iterations allowed for the optimization algorithm,
#'        controlling the depth of the search process.
#' @param ep Data frame defining pairs of variables for analysis and optimization,
#'        used to guide the optimization process by specifying which variable interactions to consider.
#'
#' @return A vector `par_all` updated with optimized model parameters for each variable pair.
#'         This vector consolidates the model parameters post-optimization, ready for model application or further analysis.
#'
#' @keywords internal
#' @importFrom stringr str_split
optimize_spatial_parameters <- function(par_all, data, names, Vi, uh, cr, max_it, ep) {
  pairs <- paste(ep[,1],ep[,2], sep = "-")
  parms <- c(paste(pairs, "ax", sep = ":"), paste(names, "ci", sep = ":"),
             paste(pairs[1:length(names)], "rij", sep = ":"), 
             paste(pairs[1:length(names)], "vij", sep = ":"))
  optimized_par <- optim(par_all[parms], fn = loglik, data = data, parms = parms,
                         par_all = par_all, ep = ep, names = names, 
                         Vi = Vi, uh = uh, cr = cr, 
                         control = list(maxit = max_it))$par
  par_all[parms] <- optimized_par
  return(update_ax_parameters(par_all, names, compute_ax(param(par_all, names), names)))
}
#' Optimize Spatio-Temporal Parameters for Variable Pairs
#'
#' Optimizes spatio-temporal model parameters for each pair of variables to enhance the
#' log-likelihood of the observed data under the model. This function focuses on both spatial
#' and temporal interactions between pairs of variables, refining the model's ability to capture
#' complex spatio-temporal dependencies.
#'
#' @param par_all Initial comprehensive set of model parameters to be refined through optimization.
#' @param data Dataset containing spatio-temporal observations, typically a 3D array or list
#'        with dimensions representing locations, times, and variables.
#' @param names Vector containing the names of the variables considered in the model,
#'        indicating the scope of the optimization.
#' @param Vi Matrix specifying variable pairs, directing the optimization towards relevant
#'        variable interactions.
#' @param uh Matrix that combines spatial and temporal distances with additional identifiers,
#'        essential for defining the context of each data point in spatio-temporal space.
#' @param cr Initial correlation matrix offering a preliminary estimate of relationships
#'        between variables, serving as a foundation for optimization.
#' @param max_it Specifies the maximum number of iterations for the optimization process,
#'        setting a limit to ensure computational feasibility.
#' @param ep Data frame defining pairs of variables for detailed analysis, guiding the optimization
#'        process by highlighting specific interactions of interest.
#'
#' @return Updated `par_all` vector containing optimized parameters for each variable pair,
#'         representing an enhanced set of model parameters post-optimization.
#'         This updated parameter set is ready for subsequent model application or further analysis.
#'
#' @keywords internal
#' @importFrom stringr str_split
optimize_pairs_spatiotemporal <- function(par_all, data, names, Vi, uh, cr, max_it, ep) {
  pairs <- paste(ep[,1],ep[,2], sep = "-")
  # Optimize model parameters for each pair of variables using the log-likelihood function
  for (i in seq(nrow(ep))) {
    pair <- pairs[i]
    sp <- unlist(stringr::str_split(pair, "-"))
    if (sp[1] == sp[2]) {
      parms <- c(paste(sp[1], "ci", sep = ":"),paste(sp[2], "ci", sep = ":"),
                 paste(sp[1], "ai", sep = ":"),paste(sp[2], "ai", sep = ":"),
                 paste(pair, "ax", sep = ":"),
                paste(pair, "rij", sep = ":"),paste(pair, "vij", sep = ":"))
      par_all[parms] <- optim(par_all[parms], fn = loglik_pair, data = data, pair = pair, parms = parms,
                              par_all = par_all, ep = ep, names = names, 
                              Vi = Vi, uh = uh, cr = cr,
                              control = list(maxit = max_it))$par
    }else{
      pair <- paste(ep[i,1],ep[i,1], sep = "-")
      parms <- c(paste(sp[1], "ci", sep = ":"),paste(sp[2], "ci", sep = ":"),
                 paste(sp[1], "ai", sep = ":"),paste(sp[2], "ai", sep = ":"))
      #pair <- paste(ep[i,2],ep[i,2], sep = "-")
      #parms <- c(parms, paste(pair, "aij", sep = ":"))
      par_all[parms] <- optim(par_all[parms], fn = loglik_pair, data = data, pair = pairs[i], parms = parms,
                              par_all = par_all, ep = ep, names = names, 
                              Vi = Vi, uh = uh, cr = cr,
                              control = list(maxit = max_it))$par
      
    }
  }
  return(par_all)
}
#' Optimize Temporal Parameters Across All Variable Pairs
#'
#' Performs a final optimization step to refine the temporal parameters of the model, 
#' considering interactions across all variable pairs. This function aims to enhance the 
#' model's temporal dynamics by optimizing a subset of parameters that influence temporal 
#' relationships.
#'
#' @param par_all Comprehensive set of model parameters, including both spatial and temporal
#'        parameters, to be optimized in this step.
#' @param data Dataset containing spatio-temporal observations, used to evaluate the model's
#'        performance and guide the optimization process.
#' @param names Vector of variable names involved in the model, indicating the scope of the
#'        optimization across variable interactions.
#' @param Vi Matrix indicating the pairs of variables considered in the model, focusing the
#'        optimization on relevant interactions.
#' @param uh Matrix combining spatial and temporal distances along with additional identifiers,
#'        critical for contextualizing each observation in spatio-temporal analysis.
#' @param cr Initial correlation matrix providing a baseline of variable relationships, serving
#'        as a starting point for optimization.
#' @param max_it Maximum number of iterations allowed in the optimization process, setting
#'        computational bounds to ensure completion.
#' @param ep Data frame specifying pairs of variables to be analyzed, guiding the optimization
#'        towards significant interactions.
#'
#' @return Updated `par_all` vector with optimized temporal parameters, representing an
#'         improved set of model parameters after the optimization process. This vector is
#'         essential for applying the model to new data or further analysis.
#'
#' @keywords internal
optimize_temporal_parameters <- function(par_all, data, names, Vi, uh, cr, max_it, ep) {
  # Final optimization step for the subset of parameters across all variable pairs
  parms <- c("a1", "d1", "g1", "a2", "d2", "g2",
             "b1", "e1", "l1", "b2", "e2", "l2", "c", "f", "m",
             paste(names, "ai", sep = ":"), paste(names, "bi", sep = ":"),
             paste(names, "ci", sep = ":"))
  optimized_par <- optim(par_all[parms], fn = loglik, data = data, parms = parms,
                         par_all = par_all, ep = ep, names = names, 
                         Vi = Vi, uh = uh, cr = cr, 
                         control = list(maxit = max_it))$par
  par_all[parms] <- optimized_par
  return(par_all)
}
#' Estimate Geostatistical Parameters for Multivariate Spatio-Temporal Data
#'
#' This function estimates the geostatistical parameters for a multivariate space-time model
#' based on observed spatio-temporal data. It integrates various preprocessing and optimization
#' steps to derive the optimal model parameters that best fit the observed data.
#'
#' @param data A 3D array containing the observed data values across time, space, and variables.
#' @param wt_id Vector of identifiers indicating specific weather types or variable categorizations within the dataset.
#' @param max_it Maximum number of iterations allowed during the optimization process.
#' @param dates Vector of dates corresponding to the temporal observations within the dataset.
#' @param tmax The maximum temporal lag considered for the model, influencing temporal dependency estimation.
#' @param names Vector containing the names of the variables included in the analysis.
#' @param par_all (Optional) Initial or current complete set of model parameters. If not provided, parameters are initialized within the function.
#' @param coordinates Matrix containing the geographical coordinates of the spatial locations in the dataset.
#' @param n1, n2 Parameters that define the granularity for generating spatial index pairs, affecting the spatial resolution of the model.
#' @param ax Matrix of precomputed correction terms used to adjust the covariance matrix, aiding in model stabilization.
#' @param cr Initial correlation matrix representing the base relationships between variables, used as a starting point for optimization.
#' @param threshold_precip Threshold values for precipitation, used in preprocessing to distinguish between different precipitation intensities.
#'
#' @return A list containing two elements: `parm`, the parameter matrix formatted for interpretation and use in subsequent model applications, and
#' `par_all`, the complete vector of optimized model parameters.
#'
#' @keywords internal
#' @importFrom parallel mclapply

estimation_gf <- function(data, wt_id, max_it, dates, tmax, names, par_all = NULL,
                          coordinates, n1, n2, ax, cr, threshold_precip) {
  # Estimate geostatistical parameters for spatio-temporal data 
  
  # Arguments:
  #   data: A 3D array of observed values over time and space for multiple variables.
  #   wt_id: Identifiers for weather types or variables within the dataset.
  #   max_it: Maximum number of iterations for the optimization process.
  #   dates: Dates corresponding to the temporal observations in the dataset.
  #   tmax: Maximum temporal lag to consider in the model.
  #   names: Names of the variables included in the model.
  #   fixed_par: Parameters that are held fixed during the optimization process.
  #   fixed_variable: Variable name that should be fixed during parameter estimation.
  #   par_all: Initial or current set of all model parameters.
  #   coordinates: Coordinates for the spatial locations in the dataset.
  #   n1, n2: Parameters defining the granularity of spatial index pairs.
  #   ax: Precomputed correction terms for the covariance matrix.
  #   cr: Initial correlation matrix for the variables.
  #   threshold_precip: Threshold values for precipitation to be used in preprocessing.
  
  # Returns:
  # The optimal parameters 
  
  # Dimensions of the data
  Nt <- dim(data)[1]  # Number of time points
  Ns <- dim(data)[2]  # Number of spatial locations
  Nv <- dim(data)[3]  # Number of variables
  
  # Generate spatial, temporal, and variable index pairs
  Si <- generate_spatial_index_pairs(coordinates, n1=n1, n2=n2)
  Ti <- generate_temporal_index_pairs(wt_id, dates, tmax)
  Vi <- generate_variable_index_pairs(names)
  
  # Preprocess data to adjust for thresholds and compute distances
  preprocessed_data <- preprocess_data(Ti, Si, coordinates)
  uh <- preprocessed_data$uh
  uh <- cbind(uh, threshold_precip[uh[,5]], threshold_precip[uh[,6]])
  u <- preprocessed_data$u
  h <- preprocessed_data$h  
  
  
  # Initialize spatial parameters 
  par_s <- init_space_par(data = data, names = names, h = h[u == 0], uh = uh[u == 0,], max_it = max_it)
  par_s <- do.call(cbind, par_s)
  
  # Construct parameter matrix for covariance model
  ep <- generate_variable_index_pairs(names)
  pairs <- paste(ep[,1],ep[,2], sep = "-")
  
  # Check and initialize par_all if missing
  par_all <- initialize_par_all_if_missing(par_all, names, pairs, par_s, ax, cr = cr)
  
  par_all <- optimize_spatial_parameters(par_all, data, names, Vi, uh[uh[,1]==0,], cr, max_it, ep)
  
  for (v in 1:2) {
    # Optimize temporal parameters
    par_all <- optimize_temporal_parameters(par_all, data, names, Vi, uh, cr, max_it, ep)
    # Optimize spatial parameters
    par_all <- optimize_spatial_parameters(par_all, data, names, Vi, uh, cr, max_it, ep)
  }
  
  # Construct parameter and beta matrices
  par_all <- update_ax_parameters(par_all, names, compute_ax(param(par_all, names), names))
  parm <- param(par_all, names)
  beta <- compute_beta(parm, names, cr)
  beta <- sapply(1:nrow(ep), function(i) beta[ep[i,1], ep[i,2]])
  par_all[1:length(beta)] <- beta
  parm <- param(par_all, names)
  
  return(list(parm = parm, par_all = par_all))
}

#' Select Uniformly Distributed Points Indices
#'
#' This function selects indices of points from a given set of coordinates that are uniformly distributed across the spatial domain.
#' The selection process involves performing k-means clustering on the coordinates and then selecting the closest data point to each cluster center.
#'
#' @param coordinates A matrix or data frame of geographical coordinates (longitude and latitude).
#' @param N The number of uniformly distributed points to select.
#'
#' @return A vector of indices corresponding to the selected uniformly distributed points.
#'
#' @keywords internal
#' @importFrom stats kmeans
selectUniformPointsIndices <- function(coordinates, N) {

  colnames(coordinates) <- c("lon", "lat")
  # K-means clustering
  clusters <- stats::kmeans(coordinates, centers = N)
  
  # For each cluster center, find the index of the closest original data point
  indices <- sapply(1:N, function(cluster_num) {
    subset_indices <- which(clusters$cluster == cluster_num)
    distances <- sqrt((coordinates$lon[subset_indices] - clusters$centers[cluster_num, "lon"])^2 + 
                        (coordinates$lat[subset_indices] - clusters$centers[cluster_num, "lat"])^2)
    return(subset_indices[which.min(distances)])
  })
  
  return(indices)
}

#' Select Points Based on Probabilities
#'
#' Selects a set of points from a dataset based on probabilities calculated from their spatial distances to a specific reference point.
#' The function calculates the probability for each point as inversely proportional to its squared distance from the reference point, ensuring a higher chance of selecting closer points.
#'
#' @param coordinates A matrix or data frame containing the coordinates of points.
#' @param betaIndex The index of the reference point in the coordinates matrix/data frame.
#' @param v The number of points to select, including the reference point.
#'
#' @return A vector of selected indices from the coordinates matrix/data frame, including the index of the reference point.
#'
#' @keywords internal


selectPoints <- function(coordinates, betaIndex, v) {
  if (v > nrow(coordinates) || betaIndex > nrow(coordinates)) {
    stop("Invalid v or betaIndex")
  }
  
  calculateProbability <- function(alphaIndex) {
    1 / (sum((coordinates[alphaIndex,] - coordinates[betaIndex,])^2) + 1)
  }
  
  probabilities <- sapply(1:nrow(coordinates), calculateProbability)
  
  selectedIndices <- unique(c(betaIndex, sample(1:nrow(coordinates), size = v, prob = probabilities)))
  return(selectedIndices)
}
#' Generate Spatial Index Pairs
#'
#' Generates pairs of spatial indices based on a subset of uniformly distributed points and additional points selected based on proximity.
#' This function first selects a set of n1 points uniformly distributed across the spatial domain. For each of these points,
#' it then selects n2 points based on their proximity, ensuring a diverse yet focused selection of spatial pairs for further analysis.
#'
#' @param coordinates A matrix or data frame containing the coordinates of points.
#' @param n1 The number of points to uniformly distribute across the spatial domain.
#' @param n2 The number of points to select based on proximity to each of the n1 points.
#'
#' @return A matrix where each row represents a pair of indices corresponding to the spatial location pairs selected for analysis.
#'
#' @keywords internal
#' @importFrom stringr str_split
generate_spatial_index_pairs <- function(coordinates,n1, n2) {
  D = as.matrix(dist(coordinates))
  Ns = nrow(coordinates)
  rs <- selectUniformPointsIndices(coordinates, n1)
  Si <- sapply(rs, function(p1) {
    
    rs1 = selectPoints(coordinates, p1, v=n2)
    
    sapply(rs1, function(p2) {
      return(paste(min(c(p1, p2)), max(c(p1, p2)), sep = "-"))
    })
  })
  Si <- unique_elements(Si)
  Si <- matrix(unlist(lapply(stringr::str_split(Si, pattern = "-"), as.numeric)), byrow = T, nrow = length(Si))
  return(Si)
}
#' Generate Temporal Index Pairs
#'
#' Creates pairs of temporal indices based on specified time lags up to a maximum time lag (tmax). This function is useful for
#' analyzing temporal relationships and covariances at different lags in spatio-temporal data.
#'
#' @param wt_id Identifiers (indices) corresponding to specific time points in the dataset.
#' @param dates Vector of dates corresponding to the time dimension in the dataset.
#' @param tmax Maximum temporal lag for which pairs are to be generated.
#'
#' @return A matrix with each row representing a pair of temporal indices and their corresponding time lag. The columns 't1' and 't2'
#' represent the indices of the paired time points, and 'u' represents the time lag between them.
#'
#' @keywords internal

generate_temporal_index_pairs <- function(wt_id,dates, tmax) {
  Ti = lapply(0:tmax, function(i){
    Ti = cbind(wt_id-i,wt_id,i)
    diff = dates[wt_id]-dates[wt_id-i]
    Ti = Ti[diff==i,]
  })
  Ti = do.call(rbind, Ti)
  colnames(Ti) = c("t1", "t2", "u")
  return(Ti)
}
#' Generate Variable Index Pairs
#'
#' Generates a data frame containing all unique pairs of variable names, including pairs of the same variable. This is used for
#' generating covariance matrices and parameter matrices where interactions between different variables (or the same variable) are
#' considered.
#'
#' @param names A vector of variable names for which index pairs are to be generated.
#'
#' @return A data frame with two columns ('v1' and 'v2') listing all unique pairs of variables. The data frame includes both
#' pairs of different variables and pairs of the same variable.
#'
#' @keywords internal

generate_variable_index_pairs <- function(names) {
  u1 = sapply(names, function(v1) sapply(names, function(v2) v1))
  u2 = sapply(names, function(v1) sapply(names, function(v2) v2))
  ep = data.frame(cbind(u1[!upper.tri(u1)],u2[!upper.tri(u2)]))
  names(ep) = c("v1","v2")
  ep = rbind(ep[ep$v1==ep$v2,], ep[!ep$v1==ep$v2,])
  return(ep)
}
#' Unique Elements Extraction
#'
#' Extracts unique elements from a vector or list, effectively removing duplicates.
#'
#' @param Si A vector or list from which unique elements need to be extracted.
#'
#' @return A vector containing only unique elements from the input.
#'
#' @keywords internal

unique_elements <- function(Si) {
  Si <- unlist(Si)
  Si <- Si[!duplicated(Si)]
  return(Si)
}
#' Data Preprocessing for Spatio-Temporal Models
#'
#' Performs preprocessing on spatio-temporal data to prepare for model fitting, including computing spatial distances and temporal lags.
#'
#' @param Ti A matrix containing temporal index pairs and their respective lags.
#' @param Si A matrix containing spatial index pairs.
#' @param coordinates A matrix of spatial coordinates for each location.
#'
#' @return A list containing the preprocessed data, including a combined matrix of temporal lags, spatial distances, and indices (`uh`), as well as separate vectors of temporal lags (`u`) and spatial distances (`h`).
#'
#' @keywords internal

preprocess_data <- function(Ti, Si, coordinates) {
  e <- expand.grid(1:nrow(Ti), 1:nrow(Si))
  u <- Ti[e[, 1], 3] 
  h <- ds(Si[e[, 2], 1], Si[e[, 2], 2],coordinates)
  uh <- cbind(u, h, Ti[e[, 1], 1], Ti[e[, 1], 2], Si[e[, 2], 1], Si[e[, 2], 2])
  return(list(uh = uh, u = u, h = h))
}

#' Estimation of Gaussian Field Parameters
#'
#' Estimates the parameters of a Gaussian field model for each weather type across spatial and temporal dimensions of weather data.
#'
#' @param data A 3D array containing weather data with dimensions [time, location, variable].
#' @param wt Vector of weather type classifications for each time point in the data.
#' @param names Vector of variable names in the data array.
#' @param coordinates A matrix of geographic coordinates for the locations in the data.
#' @param tmax Maximum time lag for temporal analysis.
#' @param max_it Maximum number of iterations for optimization procedures.
#' @param n1, n2 Parameters defining the granularity of spatial analysis.
#' @param dates Vector of dates corresponding to the time points in the data.
#' @param threshold_precip Threshold values for considering precipitation, used in data preprocessing.
#'
#' @return A list of estimated Gaussian field parameters for each weather type.
#'
#' @keywords internal
#' @importFrom geosphere distHaversine
estimate_gaussian_field_params <- function(data, wt, names, coordinates, tmax, max_it, n1, n2, dates, threshold_precip) {
  K = length(unique(wt))
  # Initialize the Gaussian field parameters storage
  gf_par <- vector(mode = "list", length = K)
  
  # Compute average pairwise correlations for each pair of variables across all locations
  
  ep <- generate_variable_index_pairs(names)
  # Estimate spatial covariance structures for each pair of variables 
  dst = sapply(1:nrow(coordinates), function(i){
    sapply(1:nrow(coordinates), function(j){
      geosphere::distHaversine(coordinates[i,], coordinates[j,])/1000
    })
  })
  vgm = lapply(1:nrow(ep), function(i){
    variable = unlist(ep[i,])
    dist = sort(unique(c(floor(dst))))
    dist = dist[seq(1, length(dist)/1.5, length.out = 2)]
    vgm = spacetime_cov(data = data[,,variable],wt_id = 2:dim(data)[1], locations = coordinates, ds = dst,
                        dates = dates, lagstime = 0, dist = dist, covgm = T)
    vgm$v = paste(variable[1], variable[2], sep = "-")
    vgm$v1 = variable[1]
    vgm$v2 = variable[2]
    
    return(vgm)
  })
  vgm = do.call(rbind, vgm)
  cr <- sapply(names, function(v1) {
    sapply(names, function(v2) {
      mean(sapply(1:dim(data)[2], function(j) cor(data[, j, v1], data[, j, v2], use = "complete.obs")), na.rm = TRUE)
    })
  })
  colnames(cr) <- rownames(cr) <- names
  # For each weather type, estimate Gaussian field parameters
  for (k in 1:K) {
    wt_id <- which(wt == k)
    wt_id <- wt_id[wt_id > tmax + 1]  
    
    #Estimate Gaussian field parameters
    gf_par[[k]] <- estimation_gf(data = data, wt_id = wt_id, max_it = max_it, dates = dates, 
                                 tmax = tmax, names = names, coordinates = coordinates, n1 = n1, 
                                 n2 = n2, ax = vgm[vgm$lagtime==0&vgm$dist==max(vgm$dist),], 
                                 cr = cr, threshold_precip = threshold_precip[[k]])$parm
  }
  
  return(gf_par)
}
