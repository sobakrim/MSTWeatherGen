#' Data Compression
#'
#' Function to compress data using the PTA3 method.
#'
#' @param data The dataset to be compressed.
#' @return The compressed data.
#'
#' @importFrom PTAk PTA3
#' @keywords internal
data_compression = function(data) {
  # Data compression function
  # Arguments:
  #   data: The dataset to be compressed
  # Returns:
  #   The compressed data 
  
  # Perform the PTA3 method on the data
  pta = PTAk::PTA3(data, nbPT = 3, nbPT2 = 1, verbose = F)
  
  # Filter components based on the presence of a "*" at the beginning of their names
  out = !substr(pta[[3]]$vsnam, 1, 1) == "*"
  
  # Calculate the contribution to total variance for the filtered components
  gct = (pta[[3]]$pct * pta[[3]]$ssX / pta[[3]]$ssX[1])[out]
  
  # Select the components based on their variance contribution
  v = pta[[1]]$v[out,]
  keep = order(gct, decreasing = TRUE)
  keep = keep[1:4] # Keep the top 4 components
  
  # Perform PCA on the selected components and remove duplicates
  pca = t(v[c(keep),])
  pca = pca[, !duplicated(pca, MARGIN = 2)]
  
  return(pca)
}
#' Weather Types Identification and Visualization
#'
#' Function to identify and visualize different weather types (WTs) based on input data.
#'
#' @param data An array with dimensions [time, locations, variables] containing meteorological data.
#' @param variables Names of variables used for constructing weather types, e.g., wind, temperature.
#' @param dates Vector of dates corresponding to the first dimension of the data.
#' @param n_wt Optional; predefined number of weather types to identify.
#' @param coordinates Data frame (or matrix) containing coordinates of the locations.
#' @param max_number_wt Optional; maximum number of weather types for automatic identification.
#' @param threshold Threshold used for defining wet days in precipitation analysis (default=0).
#' @param return_plots Logical indicating whether to return plots.
#' @param names_units Units of the variables for labeling plots.
#' @param dir Directory where to save the generated plots.
#' @return A list containing the classification of each time point into weather types,
#'   and optionally, plots illustrating the weather types.
#'
#' @keywords internal
#' @import viridis
#' @import mclust 
#' @import ggplot2
weather_types = function(data, variables,dates, n_wt = NULL,coordinates,
                         max_number_wt = NULL,threshold = 0, return_plots = T, 
                         names_units, dir){
  # Function to identify and visualize different weather types (WTs) based on input data
  # Arguments:
  #   data: An array with dimensions [time, locations, variables] containing meteorological data
  #   variables: Names of variables used for constructing weather types, e.g., wind, temperature
  #   dates: Vector of dates corresponding to the first dimension of the data
  #   n_wt: Optional; predefined number of weather types to identify
  #   coordinates: Data frame (or matrix) containing coordinates of the locations
  #   max_number_wt: Optional; maximum number of weather types for automatic identification
  #   threshold: Threshold used for defining wet days in precipitation analysis (default=0)
  #   return_plots: Logical indicating whether to return plots
  #   names_units: Units of the variables for labeling plots
  #   dir: Directory where to save the generated plots
  # Returns:
  #   A list containing the classification of each time point into weather types,
  #   and optionally, plots illustrating the weather types 
  
  # Check for necessary input for plot saving
  if(return_plots & missing(dir)) stop("provide dir to save plots")
  
  # Determine the range of weather types to consider
  if(is.null(n_wt)){
    if(missing(max_number_wt)){
      stop("Provide either n_wt or max_number_wt")
    }else{
      n_w = 1:max_number_wt
    }
  }else{
    n_w = n_wt
    plot_bic = NULL
  }
  
  # Compress data using specified method
  pca = data_compression(data)
  
  # Perform clustering to classify days into weather types
  m = mclust::Mclust(pca, G = n_w, modelNames = "VVV", verbose = F)
  cluster = m$classification
  K = length(unique(cluster))
  
  # Generate and save plots if requested
  if(return_plots){
    if(is.null(n_wt)){
      df = data.frame(bic = c(m$BIC), K = 1:max_number_wt)
      plot_bic = ggplot(df, aes(x = K, y = bic)) + geom_line(color = "red") +
        geom_point(color = "red") + theme_light() + xlab("Number of WTs") +
        ylab("BIC")
      ggsave(plot = plot_bic,paste0(dir, "plot_bic.pdf"))
      
    }
    
    plots_variables = vector(mode = "list", length  = length(variables))
    names(plots_variables) = variables
    for(i in 1:length(variables)){
      gm = colMeans(data[,,i])
      df = lapply(1:K, function(k){
        df = data.frame(lon = coordinates$longitude, lat = coordinates$latitude, x = colMeans(data[cluster==k,,i])-gm, kk = paste0("WT ",k))
      })
      df = do.call(rbind, df)
      plots_variables[[i]] = ggplot(df, aes(lon, lat )) +  borders("world", colour="black",fill= "grey",xlim=range(df$lon), ylim = range(df$lat)) +coord_cartesian(xlim=range(df$lon), ylim = range(df$lat))+
        geom_point(aes(color = x),size=3, shape=15)+scale_color_gradientn(paste("mean",names_units[i],sep = " "), colours = rev(viridis(10, option="magma")))+ theme_light() +
        theme(legend.position="top",plot.title = element_text(size = 3, hjust = 0.5), panel.spacing = unit(1, "lines"),aspect.ratio=0.9)+ xlab("Longitude (degree)") + ylab("Latitude (degree)")+
        ggtitle("") + theme_light() +facet_wrap(~kk)
      ggsave(plot =plots_variables[[i]],paste0(dir,variables[i],"_wts.pdf"), width = 28, height = 18, units = "cm")
      
    }
    
    ## Wet days frequency
    
    precip = data[,,"Precipitation"]
    ind = which((precip>threshold), arr.ind = T)
    ind1 = which(!(precip>threshold), arr.ind = T)
    precip[ind] = 1
    precip[ind1] = 0
    
    df = lapply(1:K, function(i){
      df = data.frame(lon = coordinates$longitude, lat = coordinates$latitude, x = 100*colSums(precip[cluster==i,])/length(which(cluster==i)), k = paste0("WT ", i))
    })
    df = do.call(rbind, df)
    plot_wet_days = ggplot(df, aes(lon, lat )) +  borders("world", colour="black",fill= "grey",xlim=range(df$lon), ylim = range(df$lat)) +coord_cartesian(xlim=range(df$lon), ylim = range(df$lat))+
      geom_point(aes(color = x),size=3, shape=15)+scale_color_gradientn("Frequency of wet days (%)", colours = rev(viridis(10, option="magma")))+ theme_light() +
      theme(legend.position="top",plot.title = element_text(size = 3, hjust = 0.5), panel.spacing = unit(1, "lines"),aspect.ratio=0.9)+ xlab("Longitude (degree)") + ylab("Latitude (degree)")+ggtitle("")+facet_wrap(~k)
    ggsave(plot =plot_wet_days,paste0(dir,"wet_days_wts.pdf"), width = 28, height = 18, units = "cm")
    
    ## Dry days
    
    precip = data[,,"Precipitation"]
    ind = which((precip<=threshold), arr.ind = T)
    ind1 = which(!(precip<=threshold), arr.ind = T)
    precip[ind] = 1
    precip[ind1] = 0
    
    df = lapply(1:K, function(i){
      df = data.frame(lon = coordinates$longitude, lat = coordinates$latitude, x = 100*colSums(precip[cluster==i,])/length(which(cluster==i)), k = paste0("WT ", i))
    })
    df = do.call(rbind, df)
    plot_dry_days = ggplot(df, aes(lon, lat )) +  borders("world", colour="black",fill= "grey",xlim=range(df$lon), ylim = range(df$lat)) +coord_cartesian(xlim=range(df$lon), ylim = range(df$lat))+
      geom_point(aes(color = x),size=3, shape=15)+scale_color_gradientn("Frequency of dry days (%)", colours = rev(viridis(10, option="magma")))+ theme_light() +
      theme(legend.position="top",plot.title = element_text(size = 3, hjust = 0.5), panel.spacing = unit(1, "lines"),aspect.ratio=0.9)+ xlab("Longitude (degree)") + ylab("Latitude (degree)")+ggtitle("")+facet_wrap(~k)
    ggsave(plot =plot_dry_days,paste0(dir,"dry_days_wts.pdf"), width = 28, height = 18, units = "cm")
    
    return(list(cluster = cluster, plot_bic = plot_bic, plots_variables = plots_variables, 
                plot_wet_days = plot_wet_days, plot_dry_days = plot_dry_days))
  }else{
    return(list(cluster=cluster, model = m))
  }
}
#' Transition Matrix Estimation
#'
#' Function to estimate transition matrices based on weather types.
#'
#' @param wt A vector indicating the weather type.
#' @param dates A vector of dates corresponding to each weather type.
#' @param K The total number of unique weather types.
#' @return A square matrix of size K x K where each element [i, j] represents the normalized probability
#'   of transitioning from weather type i to weather type j.
#' @keywords internal
#' @importFrom lubridate year
estimtransition = function(wt, dates, K){
  # Function to estimate transition matrices based on weather types 
  # Arguments:
  #   wt: A vector indicating the weather type 
  #   dates: A vector of dates corresponding to each weather type 
  #   K: The total number of unique weather types 
  # Returns:
  #   A square matrix of size K x K where each element [i, j] represents the normalized probability 
  #   of transitioning from weather type i to weather type j.
  
  # Initialize variables
  nwt = length(unique(wt)) # Number of unique weather types
  n = length(wt) # Total number of observations
  years = lubridate::year(dates) # Extract years from dates
  y = unique(years) # Unique years
  M = matrix(0, K, K) # Initialize transition matrix with zeros
  
  # Loop over each year
  for (i in y) {
    wty = wt[years == i] # Extract weather types for the current year
    # Loop over each day in the year except the last day
    for (j in 1:(length(wty) - 1)){
      # Increment the transition count from the current day's weather type to the next day's weather type
      M[wty[j], wty[j + 1]] = M[wty[j], wty[j + 1]] + 1
    }  
  }
  
  # Normalize the matrix rows to get transition probabilities
  M = M / apply(M, 1, sum)
  
  return(M)
}
#' Estimate Transition Matrices
#'
#' Function to estimate transition matrices between clusters over time.
#'
#' @param cluster A vector of cluster assignments for each time point.
#' @param dates A vector of dates corresponding to each time point.
#' @param nb Base number of neighbors to consider for each year, default is 30.
#' @param K The number of clusters, which determines the size of the transition matrices.
#' @return A list of transition matrices, each corresponding to a time point. These matrices contain
#'   the estimated probabilities of transitioning from one cluster to another based on the nearest neighbors' approach.
#' @keywords internal
#' @importFrom FNN get.knn
#' @importFrom lubridate year
estimate_transitions = function(cluster, dates, nb = 30, K){
  # Function to estimate transition matrices between clusters over time
  # Arguments:
  #   cluster: A vector of cluster assignments for each time point
  #   dates: A vector of dates corresponding to each time point
  #   nb: Base number of neighbors to consider for each year, default is 30. 
  #   K: The number of clusters, which determines the size of the transition matrices
  # Returns:
  #   A list of transition matrices, each corresponding to a time point. These matrices contain 
  #   the estimated probabilities of transitioning from one cluster to another based on the nearest neighbors' approach.
  
  
  # Adjust the number of neighbors 
  nb = nb * length(unique(lubridate::year(dates)))
  
  # Find the nearest neighbors for each date based on the day and month
  neighbors = FNN::get.knn(as.numeric(factor(substr(dates, 6, 11))), nb)$nn.index
  
  # Estimate transition matrices for each point in time
  tm = lapply(1:length(cluster), function(i){
    # Sort neighbors by date
    neighbors[i,] = neighbors[i, order(neighbors[i,])]
    
    # Reduce the neighbor set until there are no unique transitions left
    while (length(which(cluster[neighbors[i,]] == cluster[neighbors[i, ncol(neighbors)]])) == 1) {
      neighbors = neighbors[, -ncol(neighbors)]
    }
    
    # Estimate transition matrix for the current set of neighbors
    M = estimtransition(cluster[neighbors[i,]], dates[neighbors[i,]], K = K)
    
    # Set column and row names of the matrix
    colnames(M) = as.character(1:length(unique(cluster)))
    rownames(M) = colnames(M)
    
    # Initialize a matrix to hold transition probabilities
    ind = unique(cluster)
    tm = matrix(0, nrow = length(unique(cluster)), ncol = length(unique(cluster)))
    diag(tm) = 1  # Set diagonal to 1 (self-transition)
    
    # Fill in the transition probabilities
    for (i in 1:length(ind)) {
      for (j in 1:length(ind)) {
        tm[ind[i], ind[j]] = M[as.character(ind[i]), as.character(ind[j])]
      }
    }
    return(tm)
  })
  
  return(tm)
}
#' Simulate Weather Types
#'
#' Function to simulate weather types over specified dates using transition matrices.
#'
#' @param first_state The initial weather type from which to start the simulation.
#' @param dates_sim A vector of dates for which the weather types are to be simulated.
#' @param dates A vector of dates corresponding to the transition matrices in 'tm'.
#' @param tm A list of transition matrices, each representing the transition probabilities between weather types for a specific date.
#' @param K The total number of unique weather types or states.
#' @return A vector of simulated weather types or states for each date in 'dates_sim'.
#' @keywords internal

simulate_weathertypes = function(first_state, dates_sim, dates, tm, K) {
  # Function to simulate weather types over specified dates using transition matrices
  # Arguments:
  #   first_state: The initial weather type from which to start the simulation
  #   dates_sim: A vector of dates for which the weather types are to be simulated
  #   dates: A vector of dates corresponding to the transition matrices in 'tm'
  #   tm: A list of transition matrices, each representing the transition probabilities between weather types for a specific date
  #   K: The total number of unique weather types or states
  # Returns:
  #   A vector of simulated weather types or states for each date in 'dates_sim'
  
  # Initialize variables
  n = length(dates_sim) # Number of dates to simulate
  states = rep(first_state, n) # Initialize states vector with the first state
  # Convert dates to Date class for comparison
  dates_sim = as.Date(dates_sim)
  dates = as.Date(dates)
  
  # Simulate weather types for each day
  for (i in 2:n) {
    # Find the index of the current date in the original 'dates' vector
    id = which(dates == dates_sim[i-1])
    
    # Handle the case where the date is not found in 'dates'
    if (length(id) == 0) {
      # Find the closest date in 'dates' to the missing 'dates_sim' date
      closest = which.min(abs(dates - dates_sim[i-1]))
      warning(sprintf("Date '%s' not found in 'dates'. Using transition matrix from the closest available date: '%s'.", dates_sim[i-1], dates[closest]))
      id = closest # Use the index of the closest date
    }
    
    # Sample the next state based on the transition probabilities from the previous state
    # Adjusted to use the closest available date's transition matrix if the exact date is missing
    states[i] = sample(1:K, 1, prob = tm[[id]][states[i-1],])
  }
  
  # Return the vector of simulated states
  return(states)
}
 
