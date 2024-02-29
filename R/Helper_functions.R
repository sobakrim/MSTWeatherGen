# This file is part of MSTWeatherGen package.
#
# MSTWeatherGen is free software: you can redistribute it and/or modify it under
# the terms of the GNU General Public License as published by the Free Software
# Foundation, either version 3 of the License, or (at your option) any later version.
#
# MSTWeatherGen is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
# without even the implied warranty of MERCHANTABILITY or
# FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along with MSTWeatherGen.
# If not, see <https://www.gnu.org/licenses/>.

#' Data Compression Function
#'
#' This function performs data compression using either PCA or PTA3 methods.
#'
#' @param data The dataset to be compressed.
#' @param method A character string indicating the compression method to use.
#'               Defaults to c("pca", "pta3").
#' @return The compressed data. If "pta3" is used, returns the PCA of selected components
#'         after filtering and variance contribution analysis. If PCA is used, returns the
#'         combined and scaled PCA results from each slice of the data.
#' @examples
#' \dontrun{
#' # Assuming `my_data` is your dataset
#' compressed_data_pca = data_compression(my_data, method = "pca")
#' compressed_data_pta3 = data_compression(my_data, method = "pta3")
#' }
#' @importFrom PTAk PTA3
#' @export
data_compression = function(data, method = c("pca", "pta3")) {
  # Data compression function
  # Arguments:
  #   data: The dataset to be compressed
  #   method: A character string indicating the compression method to use
  # Returns:
  #   The compressed data

  if (method == "pta3") {
    #require(PTAk)
    # Perform the PTA3 method on the data
    pta = PTAK::PTA3(data, nbPT = 3, nbPT2 = 1)

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
  } else {
    # If method is not "pta3", perform PCA on each slice of the data and keep the first 5 components
    pca = lapply(1:dim(data)[3], function(v) {
      pc = prcomp(data[,,v], center = TRUE, scale. = TRUE)
      return(predict(pc, data[,,v])[, 1:5])
    })

    # Combine the PCA results from each slice and scale them
    pca = do.call(cbind, pca)
    pca = scale(pca)
  }


  return(pca)
}

#' Identify and Visualize Different Weather Types
#'
#' This function identifies and visualizes different weather types (WTs) based on meteorological data.
#' It supports data compression using PCA or PTA3 methods and classifies days into weather types using clustering.
#' The function can optionally return plots illustrating the weather types, including the distribution
#' of variables across weather types and the frequency of wet and dry days.
#'
#' @param data An array with dimensions time*location*variable containing meteorological data.
#' @param variables Names of variables used for constructing weather types, e.g., wind, temperature.
#' @param dates Vector of dates corresponding to the first dimension of the data.
#' @param n_wt Optional; predefined number of weather types to identify.
#' @param lonlat Data frame (or matrix) containing longitude and latitude of the locations.
#' @param max_number_wt Optional; maximum number of weather types for automatic identification.
#' @param data_compression_method Method to compress data, default options are "pca" and "pta3".
#' @param threshold Threshold used for defining wet days in precipitation analysis (default=0).
#' @param return_plots Logical indicating whether to return plots.
#' @param names_units Units of the variables for labeling plots.
#' @param dir Directory where to save the generated plots. This parameter is mandatory if return_plots is TRUE.
#' @import lubridate
#' @import PTAk
#' @import abind
#' @importFrom mclust Mclust
#' @import ggplot2
#' @import viridis
#' @return Depending on the value of return_plots, returns a list containing the classification of each time
#'         point into weather types and, optionally, plots illustrating the weather types.
#' @export

weather_types = function(data, variables,dates, n_wt = NULL,lonlat,
                         max_number_wt = NULL, data_compression_method = c("pca", "pta3"),
                         threshold = 0, return_plots = T, names_units,
                         dir){
  # Function to identify and visualize different weather types (WTs) based on input data
  # Arguments:
  #   data: An array with dimensions time*location*variable containing meteorological data
  #   variables: Names of variables used for constructing weather types, e.g., wind, temperature
  #   dates: Vector of dates corresponding to the first dimension of the data
  #   n_wt: Optional; predefined number of weather types to identify
  #   lonlat: Data frame (or matrix) containing longitude and latitude of the locations
  #   max_number_wt: Optional; maximum number of weather types for automatic identification
  #   data_compression_method: Method to compress data, default options are "pca" and "pta3"
  #   threshold: Threshold used for defining wet days in precipitation analysis (default=0)
  #   return_plots: Logical indicating whether to return plots
  #   names_units: Units of the variables for labeling plots
  #   dir: Directory where to save the generated plots
  # Returns:
  #   A list containing the classification of each time point into weather types,
  #   and optionally, plots illustrating the weather types

  # Load required packages
  # require(lubridate)
  # require(PTAk)
  # require(abind)
  # require(mclust)
  # require(ggplot2)
  # require(viridis)

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
  pca = data_compression(data, method = data_compression_method)

  # Perform clustering to classify days into weather types
  m = mclust::Mclust(pca, G = n_w, modelNames = "VVV")
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
        df = data.frame(lon = lonlat$longitude, lat = lonlat$latitude, x = colMeans(data[cluster==k,,i])-gm, kk = paste0("WT ",k))
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
      df = data.frame(lon = lonlat$longitude, lat = lonlat$latitude, x = 100*colSums(precip[cluster==i,])/length(which(cluster==i)), k = paste0("WT ", i))
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
      df = data.frame(lon = lonlat$longitude, lat = lonlat$latitude, x = 100*colSums(precip[cluster==i,])/length(which(cluster==i)), k = paste0("WT ", i))
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
#' Estimate Transition Matrices for Weather Types
#'
#' This function estimates the transition matrices based on weather types. It calculates the normalized
#' probabilities of transitioning from one weather type to another over a given set of dates. This is
#' useful for understanding the dynamics and transitions between different weather conditions over time.
#'
#' @param wt A vector indicating the weather type for each observation.
#' @param dates A vector of dates corresponding to each weather type. The dates should be in a format
#'              that allows extraction of years, as transitions are computed on an annual basis.
#' @param K The total number of unique weather types. This should be specified to ensure the transition
#'          matrix is correctly sized, even if some weather types do not occur in the data.
#' @return A square matrix of size K x K where each element [i, j] represents the normalized probability
#'         of transitioning from weather type i to weather type j. Each row in the matrix sums to 1,
#'         representing the probability distribution of transitioning to each possible weather type
#'         from a given starting type.
#' @examples
#' \dontrun{
#' # Assuming `wt_vector` is your vector of weather types and `date_vector` corresponds to dates
#' # for each weather type, and you have identified 5 unique weather types:
#' transition_matrix = estimtransition(wt_vector, date_vector, 5)
#' }
#' @export

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
  years = year(dates) # Extract years from dates
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

  # Load required packages
  require(FNN)
  require(lubridate)

  # Adjust the number of neighbors
  nb = nb * length(unique(year(dates)))

  # Find the nearest neighbors for each date based on the day and month
  neighbors = get.knn(as.numeric(factor(substr(dates, 6, 11))), nb)$nn.index

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

############# Transformation function #############

orderNorm <- function(x,left, n_logit_fit = min(length(x), 100000), ..., warn = TRUE) {
  # Function for ordered normalization Ordered Quantile normalizing transformation
  # Function inspired from "bestNormalize" package and adapted for truncated normal distributions.


  # Arguments:
  #   x: Numeric vector to be normalized.
  #   left: Numeric value specifying the left truncation.
  #   n_logit_fit: The number of points to use in the logistic regression for extrapolation,
  #                default is the lesser of the length of x or 100,000.
  #   ...: Additional arguments passed to glm function.
  #   warn: Logical flag indicating whether to issue a warning when ties are detected in the data.

  # Returns:
  #   A list with class 'orderNorm' containing the normalized data (x.t), the original data (x),
  #   and a fitted model object (fit) for future extrapolation.
  stopifnot(is.numeric(x))
  ties_status <- 0
  nunique <- length(unique(x))
  na_idx <- is.na(x)

  if (nunique < length(x)) {
    if(warn) warning('Ties in data, Normal distribution not guaranteed
                     Ties are treated by averaging\n')
    ties_status <- 1
  }
  n = length(x)
  qc = 0.5
  q.x <- (rank(x)-qc) / ((1/(1-pnorm(left)))*n) + pnorm(left)
  x.t <- crch::qcnorm(q.x,left = left)

  # fit model for future extrapolation
  # create "reduced" x with n_logit_fit equally spaced observations
  keep_idx <- round(seq(1, length(x), length.out = min(length(x), n_logit_fit)))
  x_red <- x[order(x[!na_idx])[keep_idx]]
  n_red = length(x_red)
  q_red <- (rank(x_red, na.last = 'keep', ties.method = "min")) / n_red

  # fit model for future extrapolation
  fit <- suppressWarnings(
    glm(q_red ~ x_red , family = 'binomial',
        weights = rep(n_red, n_red))
  )
  fit <- list(coef = fit$coefficients, x_red = x_red)
  val <- list(
    x.t = x.t,
    x = x,
    fit = fit
  )

  class(val) <- 'orderNorm'
  val
}
predict_binomial = function(fit, newdata){
  if(missing(newdata)) newdata = fit$x_red
  pred = fit$coef[1]+fit$coef[2]*newdata
  return(1/(1+exp(-pred)))
}
orderNorm_all <- function(data, j, lonlat, left) {
  # Function to perform normalization on a selected variable 'j' from 'data', considering spatial information
  # from 'lonlat' and a threshold 'left' for the transformation. This is particularly useful for variables with
  # many zero values, ensuring robust normalization across spatially correlated variables.

  # Arguments:
  #   data: A matrix or data frame containing the variables to be normalized. Rows represent observations,
  #         and columns represent different variables.
  #   j: The index of the variable within 'data' to be normalized.
  #   lonlat: A matrix or data frame of spatial coordinates (longitude and latitude) corresponding to each
  #           variable in 'data'. It is used to calculate spatial distances and determine proximity.
  #   left: A parameter for the orderNorm function, specifying the transformation threshold.

  # Returns:
  #   The result of applying the orderNorm function to the selected variable 'j', considering spatial proximity
  #   and handling variables with a high proportion of zero values.

  require(FNN) # Load the FNN package for k-nearest neighbors calculations.

  D = as.matrix(dist(lonlat))
  kn = order(D[j,])
  j = kn[1]
  x = data[,j]
  j = 2
  if(length(which(x==0))/length(x)<0.8){
    return(orderNorm(x[!x==0],left = left))
  }else{
    while ((length(which(x!=0))<50)) {
      x = c(x,data[,kn[j]])
      j = j + 1
    }
    return(orderNorm(x[!x==0],left = left))
  }
}

#' predict order norm Transf
#' @param object description
#' @param newdata description
#' @param inverse description
#' @param warn desc
#' @param ... more arguments
#' 
predict.orderNormTransf <- function(object,
                                    newdata = NULL,
                                    inverse = FALSE,
                                    warn = TRUE,
                                    ...) {
  stopifnot(is.null(newdata) || is.numeric(newdata))

  # Perform transformation
  if(!inverse) {
    if(is.null(newdata)) newdata <- object$x
    na_idx <- is.na(newdata)

    newdata[!na_idx] <- orderNormTransf(object, newdata[!na_idx], warn)
    return(newdata)
  }

  # Perform reverse transformation
  if (is.null(newdata)) newdata <- object$x.t

  na_idx <- is.na(newdata)
  newdata[!na_idx] <- inv_orderNorm_Transf(object, newdata[!na_idx], warn)

  return(newdata)
}
inv_orderNorm_Transf <- function(orderNorm_obj, new_points_x_t, left, warn = FALSE) {
  # Reverses the normalization or transformation applied by the orderNorm function.
  # Arguments:
  #   orderNorm_obj: An object containing details of the original normalization or transformation,
  #                  including the transformed and original data points and the fitting model.
  #   new_points_x_t: Transformed data points for which the original values are to be estimated.
  #   left: Parameter used in the original transformation, affecting the normalization method.
  #   warn: A logical flag indicating whether to issue warnings for transformations outside the observed domain.

  # Extract transformed and original data points from the orderNorm object.
  x_t <- orderNorm_obj$x.t
  old_points <- orderNorm_obj$x

  if(min(old_points)>0){
    vals <- suppressWarnings(
      approx(x_t, old_points, xout = new_points_x_t, rule = 2:1)
    )
  }else{
    vals <- suppressWarnings(
      approx(x_t, old_points, xout = new_points_x_t, rule = 1)
    )
  }
  # If predictions have been made outside observed domain
  if (any(is.na(vals$y))) {
    if(warn) warning('Transformations requested outside observed domain; logit approx. on ranks applied')

    fit <- orderNorm_obj$fit
    p <- crch::qcnorm(predict_binomial(fit), left = left)
    if(min(old_points)>0){
      l_idx = FALSE
    }else{
      l_idx <- vals$x < min(x_t, na.rm = TRUE)
    }
    h_idx <- vals$x > max(x_t, na.rm = TRUE)
    if (any(l_idx)) {
      # Solve algebraically from original transformation
      pp = crch::pcnorm(vals$x[l_idx]+ min(p, na.rm = TRUE) - min(x_t, na.rm = TRUE) , left = left)
      pp[pp==1] = 0.99999
      logits <- log(pp / (1 - pp))
      vals$y[l_idx] <-
        unname((logits - fit$coef[1] ) / fit$coef[2])
    }

    if (any(h_idx)) {
      pp = crch::pcnorm(vals$x[h_idx]+ max(p, na.rm = TRUE) - max(x_t, na.rm = TRUE) , left = left)
      pp[pp==1] = 0.99999
      logits <- log(pp / (1 - pp))
      vals$y[h_idx] <-
        unname((logits - fit$coef[1] ) / fit$coef[2])
    }
  }

  return(vals$y)
}
orderNormTransf <- function(orderNorm_obj, new_points, warn, left) {
  # Transforms new data points based on a previously fitted normalization or transformation model.
  # This function is used to apply the same transformation to new data that was applied to the original data.

  # Arguments:
  #   orderNorm_obj: An object containing the original and transformed data points, as well as the fit model.
  #   new_points: The new data points to be transformed.
  #   warn: A logical flag indicating whether to issue warnings for transformations outside the observed domain.
  #   left: A parameter affecting the normalization method used in the original transformation.

  # Extract transformed and original data points from the orderNorm object.
  x_t <- orderNorm_obj$x.t
  old_points <- orderNorm_obj$x
  vals <- suppressWarnings(
    approx(old_points, x_t, xout = new_points, rule = 1)
  )

  # If predictions have been made outside observed domain
  if (any(is.na(vals$y))) {
    if (warn) warning('Transformations requested outside observed domain; logit approx. on ranks applied')
    fit <- orderNorm_obj$fit
    p <- crch::qcnorm(predict_binomial(fit), left = left)
    l_idx <- vals$x < min(old_points, na.rm = TRUE)
    h_idx <- vals$x > max(old_points, na.rm = TRUE)

    # Check
    if (any(l_idx)) {
      xx <- data.frame(x_red = vals$x[l_idx])
      vals$y[l_idx] <- crch::qcnorm(predict_binomial(fit, newdata = xx), left = left) -
        (min(p, na.rm = TRUE) - min(x_t, na.rm = TRUE))

    }
    if (any(h_idx)) {
      xx <- data.frame(x_red = vals$x[h_idx])
      vals$y[h_idx] <- crch::qcnorm(predict_binomial(fit, newdata = xx), left = left) -
        (max(p, na.rm = TRUE) - max(x_t, na.rm = TRUE))
    }
  }

  return(vals$y)
}

############# Covariance function #############



Matern <- function(h, r, v) {
  # Calculates the Matern covariance function for a given vector of distances.

  # Arguments:
  #   h: Numeric vector of distances between points.
  #   r: Scalar range parameter of the Matern function, affecting spatial correlation decay.
  #   v: Scalar smoothness parameter of the Matern function, controlling the smoothness of the resulting field.

  # Returns:
  #   Numeric vector representing the covariance values calculated using the Matern function for the distances in h.

  rt <- (2 ^ (1 - v)) / gamma(v) * (r * abs(h)) ^ v * besselK(r * abs(h), nu = v)
  rt[h == 0] <- 1  # Ensures that the covariance at distance 0 is 1.
  return(rt)
}

Gneiting <- function(h, u, par, dij) {
  # Multivariate space-time Gneiting's covariance function
  # From the paper: https://doi.org/10.1016/j.spasta.2022.100706

  # Arguments:
  #   h: Numeric vector of spatial distances.
  #   u: Numeric vector of temporal distances.
  #   par: Numeric vector of parameters used in the covariance function.
  #   dij: correlation parameter between variable i and j.

  # Returns:
  #   Covariance value(s) calculated using Gneiting's spatio-temporal covariance model.

  if(!is.numeric(par)) par <- as.numeric(par)

  # Unpack parameters from the 'par' vector for clarity.
  a <- par[1]
  b <- par[2]
  c <- par[3]
  d <- par[4]
  e <- par[5]
  ci <- par[6]
  cj <- par[7]
  rii <- par[12]
  rjj <- par[13]
  vii <- par[14]
  vjj <- par[15]
  ax <- par[16]
  aii <- par[17]
  ajj <- par[18]
  bii <- par[19]
  bjj <- par[20]

  # Calculated intermediate parametrs for the covariance calculation.
  vij = (vii + vjj) / 2
  rij = sqrt((rii^2 + rjj^2) / 2)
  aij = sqrt((aii^2 + ajj^2) / 2)
  bij = sqrt((bii^2 + bjj^2) / 2)

  eij <- dij * ((rii^vii * rjj^vjj) / rij^(2*vij)) * (gamma(vij) / (gamma(vii)^(1/2) * gamma(vjj)^(1/2))) *
    sqrt((1-ci^2)*(1-cj^2)) * (aii^(1/2) * ajj^(1/2)) / aij

  muij <- ((a * abs(u))^(2*b) + 1)^c - (ci*cj * ((d * abs(u))^(2*e) + 1)^(-c))

  A1 <- eij / muij
  A2 <- Matern(h, r = (rij^2 / muij)^(1/2), v = vij) * exp(-aij * abs(u))
  A3 <- ax * ((bii^(1/2) * bjj^(1/2)) / bij) * exp(-bij * abs(u))

  return(A1 * A2 + A3)
}

param <- function(par, names) {
  # Function to construct a data frame with covariance parameters

  # Arguments:
  #   par: Named vector of parameters
  #   names: Character vector of variable names, such as "temperature", "wind", etc.

  # Returns:
  #   A data frame with columns for each combination of variables and their associated parameters,
  #   filled with values specified in the 'par' vector.

  # Generate all possible pairs of variable names, including self-pairs, for parameter definitions
  ep <- generate_variable_index_pairs(names)
  pairs <- paste(ep[,1],ep[,2], sep = "-")
  # Initialize a data frame to be populated with parameter values
  J <- length(pairs)
  u <- data.frame(v1 = ep$v1, v2 = ep$v2, stringsAsFactors = FALSE)

  # Assign common temporal parameters to all pairs
  u$a <- par["a"]
  u$b <- par["b"]
  u$c <- par["c"]
  u$d <- par["d"]
  u$e <- par["e"]

  # Loop through each pair to populate the data frame with corresponding parameter values
  for (i in seq_len(J)) {
    # Extract and assign specific parameters for each pair based on naming convention
    u$ci[i] <- par[paste(u$v1[i],  "ci", sep = ":")]
    u$cj[i] <- par[paste(u$v2[i],  "ci", sep = ":")]

    u$aij[i] <- par[paste(paste(u$v1[i],u$v2[i],  sep = "-"),  "aij", sep = ":")]
    u$bij[i] <- par[paste(paste(u$v1[i],u$v2[i],  sep = "-"),  "bij", sep = ":")]

    u$rij[i] <- par[paste(paste(u$v1[i],u$v2[i],  sep = "-"),  "rij", sep = ":")]
    u$vij[i] <- par[paste(paste(u$v1[i],u$v2[i],  sep = "-"),  "vij", sep = ":")]

    u$rii[i] <- par[paste(paste(u$v1[i],u$v1[i],  sep = "-"),  "rij", sep = ":")]
    u$rjj[i] <- par[paste(paste(u$v2[i],u$v2[i],  sep = "-"),  "rij", sep = ":")]
    MSTSWG_simulation_season
    u$vii[i] <- par[paste(paste(u$v1[i],u$v1[i],  sep = "-"),  "vij", sep = ":")]
    u$vjj[i] <- par[paste(paste(u$v2[i],u$v2[i],  sep = "-"),  "vij", sep = ":")]

    u$ax[i] <- par[paste(paste(u$v1[i],u$v2[i],  sep = "-"),  "ax", sep = ":")]
    u$aii[i] <- par[paste(paste(u$v1[i],u$v1[i],  sep = "-"),  "aij", sep = ":")]
    u$ajj[i] <- par[paste(paste(u$v2[i],u$v2[i],  sep = "-"),  "aij", sep = ":")]

    u$bii[i] <- par[paste(paste(u$v1[i],u$v1[i],  sep = "-"),  "bij", sep = ":")]
    u$bjj[i] <- par[paste(paste(u$v2[i],u$v2[i],  sep = "-"),  "bij", sep = ":")]

    u$dij[i] <- par[paste(paste(u$v1[i],u$v2[i],  sep = "-"),  "dij", sep = ":")]
  }
  return(u)
}

compute_beta <- function(parm, names, cr) {
  # Function for calculating correlations (dij) based on the Gneiting function

  # Arguments:
  #   parm: A data frame or similar structure containing parameters for the Gneiting function,
  #   names: A vector of variable names (e.g., "temperature", "wind") to calculate correlations between.
  #   cr: A matrix containing initial correlation values between the variables.

  # Returns:
  #   A symmetric matrix (beta) where each element [i, j] represents the correlation coefficient
  #   between variables i and j, adjusted based on the Gneiting function and a correction term.

  J = length(names)  # Number of variables
  beta <- matrix(0, ncol = J, nrow = J)  # Initialize the beta matrix with zeros
  colnames(beta) <- rownames(beta) <- names  # Set the row and column names of the matrix to variable names

  # Create a map to fetch parameters quickly using a two-level list structure
  parm_map <- split(parm, list(parm$v1, parm$v2))

  # Function to retrieve parameters for a given pair of variables v1 and v2
  get_parameters <- function(v1, v2) {
    # Attempt to fetch parameters based on the naming convention, handling both v1,v2 and v2,v1 cases
    if (exists(paste0(v1, ".", v2), parm_map)) {
      par <- as.numeric(parm_map[[paste0(v1, ".", v2)]][-c(1, 2)])  # Exclude the first two elements (variable names)
    } else {
      par <- as.numeric(parm_map[[paste0(v2, ".", v1)]][-c(1, 2)])
    }
    return(par)
  }

  # Iterate over pairs of variables to compute and assign the beta values
  for (j in seq_along(names)) {
    for (k in seq(j, J)) {
      v1 <- names[j]  # First variable in the pair
      v2 <- names[k]  # Second variable in the pair
      par <- get_parameters(v1, v2)  # Retrieve parameters for the current pair

      # Calculate the correlation coefficient using the Gneiting function and correction term
      cc = Gneiting(0, 0, par, dij = 1)  # Gneiting function calculation for the pair
      ax = par[16] * (par[19]^(1/2) * par[20]^(1/2)) / sqrt((par[19]^2 + par[20]^2) / 2)  # Correction term calculation
      beta_val <- (cr[v1, v2] - ax) / (cc - ax)  # Adjusted correlation coefficient

      beta[v1, v2] <- beta[v2, v1] <- beta_val  # Symmetric assignment to ensure the matrix is symmetric
    }
  }

  return(beta)
}
compute_ax <- function(parm, names, cr) {
  # Function to calculate correction terms (ax) for spatial covariance model

  # Arguments:
  #   parm: A data frame or similar structure containing parameters for the Gneiting function
  #   names: A vector of variable names (e.g., "temperature", "wind")
  #   cr: A matrix containing initial correlation values between the variables.

  # Returns:
  #   A symmetric matrix (ax) where each element [i, j] represents the correction term
  #   between variables i and j, computed based on the beta coefficients and the Gneiting function.

  J = length(names)  # Number of variables
  ax <- matrix(0, ncol = J, nrow = J)  # Initialize the ax matrix with zeros
  colnames(ax) <- rownames(ax) <- names  # Set the row and column names of the matrix to variable names

  # Create a map for quick parameter fetching using a two-level list structure
  parm_map <- split(parm, list(parm$v1, parm$v2))

  # Compute the beta matrix using nearPD for ensuring positive definiteness
  beta <- Matrix::nearPD(compute_beta(parm, names, cr))$mat

  # Function to retrieve parameters for a given pair of variables v1 and v2
  get_parameters <- function(v1, v2) {
    # Attempt to fetch parameters based on naming convention, handling both v1,v2 and v2,v1 cases
    if (exists(paste0(v1, ".", v2), parm_map)) {
      par <- as.numeric(parm_map[[paste0(v1, ".", v2)]][-c(1, 2)])  # Exclude variable names from parameters
    } else {
      par <- as.numeric(parm_map[[paste0(v2, ".", v1)]][-c(1, 2)])
    }
    return(par)
  }

  # Iterate over pairs of variables to compute and assign the ax values
  for (j in seq_along(names)) {
    for (k in seq(j, J)) {
      v1 <- names[j]  # First variable in the pair
      v2 <- names[k]  # Second variable in the pair
      par <- get_parameters(v1, v2)  # Retrieve parameters for the current pair

      # Calculate the ax value using the beta matrix and the Gneiting function
      cc <- Gneiting(0, 0, par, dij = 1)  # Gneiting function calculation for the pair
      ax_val <- (cr[v1, v2] - (beta[v1, v2] * cc)) / (1 - beta[v1, v2])

      ax[v1, v2] <- ax[v2, v1] <- ax_val  # Symmetric assignment
    }
  }

  return(ax)
}

compute_ax <- function(parm, names) {
  # Extract a matrix of correction terms ('ax') for a set of variables based on parameters
  # provided in 'parm'.

  # Arguments:
  #   parm: A data frame or list containing the model parameters.
  #   names: A vector of variable names.

  # Returns:
  #   A matrix where each element [i, j] represents the correction term ('ax') between the ith and jth variables

  ax = sapply(names, function(v1){
    sapply(names, function(v2){
      parm$ax[parm$v1==v1&parm$v2==v2|parm$v1==v2&parm$v2==v1]
    })
  })
  return(ax)
}

loglik_pair <- function(par, parms, pair, par_all, data, names, Vi, h, u, uh, ep, cr) {
  # Function to compute the log-likelihood for a pair of variables using the Gneiting spatio-temporal
  # covariance model.

  # Arguments:
  #   par: Current parameters being optimized.
  #   parms: Indices or names of parameters in 'par' to be updated.
  #   pair: A string indicating the pair of variables (e.g., "temperature-wind") being analyzed.
  #   par_all: Complete set of parameters for the model.
  #   data: 3D array of observed data, with dimensions corresponding to  times, locations, and variables.
  #   names: Vector of variable names, indicating the variables' names (e.g., "temperature" and "wind").
  #   Vi: Matrix where each line corresponds for a possible combination of variables in "names"
  #   h: Vector of spatial distances for the pair.
  #   u: Vector of temporal distances for the pair.
  #   uh: Matrix containing pairs of spatial and temporal distances, and additional information.
  #   ep: A matrix or data frame defining pairs of variables.
  #   beta: Matrix of beta coefficients, precomputed.
  #   cr: Correlation matrix, initial or base correlations between variables.

  # Returns:
  #   The log-likelihood value for the given pair of variables based on the current model parameters.

  require(VGAM)

  J = length(names)  # Number of variables
  pairs = paste(ep[,1], ep[,2], sep = "-")  # Constructing pairs from 'ep' data frame

  par_all[parms] = par  # Update specific parameters in the complete set

  par = par_all  # Use the updated parameter set for computations
  sp = unlist(strsplit(pair, "-"))  # Split the pair string to individual variables
  v = which(Vi[,1] == sp[1] & Vi[,2] == sp[2])  # Find the index of the pair in 'Vi'

  # Update and compute model parameters
  parm <- param(par, names)
  beta <- compute_beta(parm, names, cr)  # Compute beta coefficients
  ax <- compute_ax(parm, names)  # Compute ax correction terms

  # Attempt Cholesky decompositions for 'ax' and 'beta', checking for positive definiteness
  ae <- try(chol(ax), silent = TRUE)
  be <- try(chol(beta), silent = TRUE)

  if (!is.character(be) & (!is.character(ae) | length(which(ax == 0)) == length(names)^2)) {
    # Proceed if both 'ax' and 'beta' matrices are valid for further computations

    # Map parameters to each variable pair in 'Vi'
    parmm = lapply(1:nrow(Vi), function(v) {
      as.numeric(parm[(parm$v1 == Vi[v,1] & parm$v2 == Vi[v,2]) | (parm$v1 == Vi[v,2] & parm$v2 == Vi[v,1]), ][-c(1,2)])
    })
    u = uh[,1]
    h = uh[,2]

    # Initializing log-likelihood components
    l1 = l2 = l3 = l4 = 0
    par = parmm[[v]]  # Parameters for the current pair

    # Parameter constraints check; return a large value if constraints are violated
    if (any(par[c(1,3,4,6,7,8,9,10,11)] < 0) | any(par[c(2,3,5,6,7)] > 1)) {
      return(abs(rnorm(1)) * 1e+20)
    } else {
      # Compute covariance for the pair
      cij = Gneiting(h = h, u = u, par = par, dij = beta[sp[1], sp[2]])
      delta = 1 - cij^2
      # Extract observed values for the pair from 'data'
      v1 = data[,,Vi[v,1]]; v1 = v1[cbind(uh[,3], uh[,5])]
      v2 = data[,,Vi[v,2]]; v2 = v2[cbind(uh[,4], uh[,6])]

      # Detailed computation of log-likelihood components for various cases
      # Identifying cases based on the variable type (e.g., Precipitation) and the presence of zero values

      dz = !(h==0 & u ==0 & Vi[v,1] == Vi[v,2])
      cij = cij[dz]; delta = delta[dz]; v1= v1[dz]; v2 = v2[dz]; uh = uh[dz,]
      id1 = (v1 == 0)&(!v2 == 0)&(sp[1]=="Precipitation"); id2 = (!v1 == 0)&(v2 == 0)&(sp[2]=="Precipitation")
      id4 = (!v1 == 0)&(!v2 == 0) ; id3 = (v1 == 0)&(v2 == 0)&(sp[1]=="Precipitation")&(sp[2]=="Precipitation")


      # Adjustments for infinite values in setting a practical lower limit
      uh[, 8][which(uh[, 8] == -Inf)] = -2.282295  # Adjusting for log transform lower bounds
      uh[, 7][which(uh[, 7] == -Inf)] = -2.282295

      # Case 1: Both variables have non-zero values
      if (!length(which(id4 == TRUE)) == 0) {
        l4 = sum((-1 / 2) * (log(delta[id4]) + (v1[id4]^2 - (2 * cij[id4] * v1[id4] * v2[id4]) + v2[id4]^2) / delta[id4]))
      }

      # Case 2: First variable is non-zero and the second is zero, and the second is Precipitation
      if (!length(which(id2 == TRUE)) == 0) {
        l2 = sum(log(pnorm((uh[id2, 8] - cij[id2] * v1[id2]) / sqrt(delta[id2]))))
      }

      # Case 3: First variable is zero and the second is non-zero, and the first is Precipitation
      if (!length(which(id1 == TRUE)) == 0) {
        l1 = sum(log(pnorm((uh[id1, 7] - cij[id1] * v2[id1]) / sqrt(delta[id1]))))
      }

      # Case 4: Both variables are zero, and both are Precipitation
      if (!length(which(id3 == TRUE)) == 0) {
        l3 = try(sum(log(pbinorm(uh[id3, 7], uh[id3, 8], var1 = 1, var2 = 1, cov12 = cij[id3]))), silent = TRUE)
        if (is.character(l3)) l3 = -abs(rnorm(1)) * 1e+20  # Handle errors in computing bivariate normal CDF
      }

      # The negative log-likelihood
      return(-(l1 + l2 + l3 + l4))

    }
  } else {
    # Return a large value if Cholesky decomposition fails, indicating non-positive definiteness
    return(abs(rnorm(1)) * 1e+20)
  }
}


loglik <- function(par, parms, par_all, data, names, Vi, h, u, uh, ep, cr) {
  # Calculates the total log-likelihood for spatial or spatio-temporal data based on model parameters,
  # integrating across all variable pairs to support model fitting and optimization.

  # Arguments:
  #   par: Vector of current parameter estimates to be optimized.
  #   parms: Indices or names of parameters in 'par' that are to be updated.
  #   par_all: Complete vector of all model parameters, including those not currently being optimized.
  #   data: 3D array of observed data, with dimensions for locations, times, and variables.
  #   names: Vector of variable names corresponding to the third dimension of 'data'.
  #   Vi: Matrix where each line corresponds for a possible combination of variables in "names"
  #   h: Vector of spatial distances for the pairs being analyzed.
  #   u: Vector of temporal distances for the pairs being analyzed.
  #   uh: Matrix containing combined spatial and temporal distances along with additional identifiers.
  #   ep: Data frame defining pairs of variables for analysis.
  #   beta: Precomputed beta coefficients matrix for all pairs.
  #   cr: Initial correlation matrix for the variables.

  # Returns:
  #   The total log-likelihood value for the data based on the current set of parameters.

  require(VGAM)  # Load VGAM for statistical modeling capabilities.

  J = length(names)  # Number of variables in the analysis.
  pairs = paste(ep[,1], ep[,2], sep = "-")  # Construct pairs from 'ep' for parameter naming.

  # Update parameter names and values based on current estimates.
  names(par_all) <- c(paste(pairs, "dij", sep = ":"), "a", "b", "c", "d", "e",
                      paste(names, "ci", sep = ":"), paste(pairs, "aij", sep = ":"),
                      paste(pairs, "bij", sep = ":"), paste(pairs, "rij", sep = ":"),
                      paste(pairs, "vij", sep = ":"), paste(pairs, "ax", sep = ":"))
  par_all[parms] = par  # Update specified parameters.

  par = par_all  # Use updated parameters for calculations.
  parm = param(par, names)  # Organize parameters for each pair.
  beta = compute_beta(parm, names, cr)  # Compute beta coefficients matrix.
  ax = compute_ax(parm, names)  # Compute ax correction terms matrix.

  # Attempt Cholesky decomposition to ensure positive definiteness.
  ae <- try(chol(ax), silent = TRUE)
  be <- try(chol(beta), silent = TRUE)

  if (!is.character(be) & (!is.character(ae) | length(which(ax == 0)) == length(names)^2)) {
    # Proceed if both 'ax' and 'beta' matrices are valid for further computations.

    # Map parameters to each variable pair in 'Vi'.
    parmm = lapply(1:nrow(Vi), function(v){
      as.numeric(parm[(parm$v1==Vi[v,1] & parm$v2 ==Vi[v,2])|(parm$v1==Vi[v,2] & parm$v2 ==Vi[v,1]),][,-c(1,2)])
    })
    u = uh[,1]
    h = uh[,2]

    # Parallel computation of log-likelihood for each pair using mclapply (if multicore is intended, else lapply).
    ll = mclapply(1:nrow(Vi), function(v) {
      # Initialize log-likelihood components for the current pair.
      l1 = l2 = l3 = l4 = 0
      par = parmm[[v]]  # Parameters for the current pair.

      # Validate parameter constraints; return a large penalty if violated.
      if (any(par[c(1,3,4,6,7,8,9,10,11)] < 0) | any(par[c(2,3,5,6,7)] > 1)) {
        return(-abs(rnorm(1)) * 1e+20)
      } else {
        # Calculate pairwise log-likelihood using Gneiting function and parameter adjustments.
        cij = Gneiting(h = h, u = u, par = par, dij = beta[Vi[v,1], Vi[v,2]])
        delta = 1 - cij^2
        v1 = data[,,Vi[v,1]]; v1 = v1[cbind(uh[,3], uh[,5])]
        v2 = data[,,Vi[v,2]]; v2 = v2[cbind(uh[,4], uh[,6])]
        dz = !(h == 0 & u == 0 & Vi[v,1] == Vi[v,2])
        cij = cij[dz]; delta = delta[dz]; v1 = v1[dz]; v2 = v2[dz]; uh = uh[dz,]

        # Detailed computations for log-likelihood components based on variable presence and types.
        id1 = (v1 == 0)&(!v2 == 0)&( Vi[v,1]=="Precipitation"); id2 = (!v1 == 0)&(v2 == 0)&( Vi[v,2]=="Precipitation")
        id4 = (!v1 == 0)&(!v2 == 0) ; id3 = (v1 == 0)&(v2 == 0)&( Vi[v,1]=="Precipitation")&( Vi[v,2]=="Precipitation")
        uh[,8][which(uh[,8]==-Inf)] = -2.282295
        uh[,7][which(uh[,7]==-Inf)] = -2.282295

        # l1: Case where the first variable is zero and the second is non-zero
        # where the first variable is "Precipitation"
        if (!length(which(id1 == TRUE)) == 0) {
          l1 = sum(log(pnorm((uh[id1, 7] - cij[id1] * v2[id1]) / sqrt(delta[id1]))))
        }

        # l2: Case where the first variable is non-zero and the second is zero
        # where the second variable is "Precipitation"
        if (!length(which(id2 == TRUE)) == 0) {
          l2 = sum(log(pnorm((uh[id2, 8] - cij[id2] * v1[id2]) / sqrt(delta[id2]))))
        }

        # l3: Case where both variables are zero and both are "Precipitation"
        if (!length(which(id3 == TRUE)) == 0) {
          l3 = sum(log(pbinorm(uh[id3, 7], uh[id3, 8], var1 = 1, var2 = 1, cov12 = cij[id3])))
        }

        # l4: Case where both variables have non-zero values
        if (!length(which(id4 == TRUE)) == 0) {
          l4 = sum((-1 / 2) * (log(delta[id4]) + (v1[id4]^2 - (2 * cij[id4] * v1[id4] * v2[id4]) + v2[id4]^2) / delta[id4]))
        }

        return(l1 + l2 + l3 + l4)
      }
    })

    # Sum and negate the log-likelihood contributions from all pairs.
    return(-sum(unlist(ll)))
  } else {
    # Return a large penalty if Cholesky decomposition fails, indicating issues with matrix definiteness.
    return(abs(rnorm(1)) * 1e+20)
  }
}
loglik_spatial <- function(par, data, h, uh, v) {
  # Function to calculate the log-likelihood for spatial data using the Matérn covariance function.
  # This is crucial for estimating geostatistical parameters in spatial models.
  # Arguments:
  #   par: Vector of parameters for the Matérn covariance function, specifically the range (par[1]) and
  #        smoothness (par[2]). Both parameters must be positive.
  #   data: A 3D array of observed data, with dimensions representing different spatial locations, time points,
  #         and variables.
  #   h: A vector of spatial distances between locations, used in the covariance function.
  #   uh: A matrix with indices to map the data array to spatial-temporal pairs for which the log-likelihood is calculated.
  #   v: The index of the variable in 'data' for which the log-likelihood is being calculated.
  # Returns:
  #   The log-likelihood value of the spatial data under the specified Matérn covariance model

  # Penalize negative parameters to enforce model constraints.
  if(par[1] < 0 | par[2] < 0) {
    return(abs(rnorm(1)) * 1e+20)
  } else {
    # Initialize components of the log-likelihood calculation.
    l1 = l2 = l3 = l4 = 0

    # Compute covariances using the Matérn function based on spatial distances 'h'.
    cij = Matern(h, r = par[1], v = par[2])
    delta = 1 - cij^2

    # Extract paired observations for variable 'v' based on spatial-temporal indices in 'uh'.
    v1 = data[,,v]; v1 = v1[cbind(uh[,3], uh[,5])]
    v2 = data[,,v]; v2 = v2[cbind(uh[,4], uh[,6])]

    # Exclude stationary points to focus on spatial variation.
    dz = !(h == 0)
    cij = cij[dz]; delta = delta[dz]; v1 = v1[dz]; v2 = v2[dz]

    # Identify scenarios based on zero and non-zero observations and compute respective components.
    id1 = (v1 == 0) & (!v2 == 0)
    id2 = (!v1 == 0) & (v2 == 0)
    id4 = (!v1 == 0) & (!v2 == 0)
    id3 = (v1 == 0) & (v2 == 0)

    # Aggregate log-likelihood components considering the identified scenarios.
    if(!length(which(id4==T))==0){
      l4 = sum((-1/2)*(log(delta[id4])+(v1[id4]^2-(2*cij[id4]*v1[id4]*v2[id4])+v2[id4]^2)/delta[id4]))
    }else if(!length(which(id2==T))==0){
      l2 = sum(log(pnorm((-cij[id2]*v1[id2])/sqrt(delta[id2]))))
    }else if(!length(which(id1==T))==0){
      l1 = sum(log(pnorm((-cij[id1]*v2[id1])/sqrt(delta[id1]))))
    }else if(!length(which(id3==T))==0){
      l3 = sum(pbinorm(0,0,var1 = 1, var2 = 1, cov12 = cij[id3]))
    }

    # Return the aggregated negative log-likelihood, adjusting for errors or infinite values.
    ll = try(-(l1 + l2 + l3 + l4), silent = TRUE)
    if(is.character(ll) || is.infinite(ll)) ll = abs(rnorm(1)) * 1e+20
    return(ll)
  }
}

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
  par = mclapply(names, function(v) {
    optim(
      par = c(1, 0.1),  # Initial parameter guesses
      fn = loglik_spatial,     # Objective function to minimize (negative log-likelihood)
      data = data,
      v = v,            # Current variable being optimized
      h = h,
      uh = uh,
      control = list(maxit = max_it)  # Optimization control settings
    )$par  # Extract the optimized parameters
  })

  return(par)
}
spacetime_cov <- function(data, wt_id, locations, ds = NULL, dates, lagstime, dist, covgm = TRUE) {
  # Computes spatial and temporal covariances for spatio-temporal data.

  # Arguments:
  #   data: A 3D array with dimensions time*location*variable
  #   wt_id: Weather type indices
  #   locations: Spatial locations of the data points.
  #   ds: Precomputed distance matrix. If NULL, distances are computed from 'locations'.
  #   dates: Time points corresponding to the observations.
  #   lagstime: Vector of time lags for which covariances are computed.
  #   dist: Vector of spatial distances for which covariances are computed.
  #   covgm: Flag indicating whether a cross-covariances should be calculated (default TRUE).

  # Returns:
  #   A data frame of computed covariances across specified spatial distances and time lags.

  # Validate input data structure
  if (covgm && length(dim(data)) < 3) {
    stop("data must be a 3D array when 'covgm' flag is set.")
  }

  # Compute distance matrix if not provided
  if (is.null(ds)) {
    ds <- round(as.matrix(dist(locations)), 3)
  }

  # Compute covariance for zero distance to establish a baseline
  id = cbind(wt_id, wt_id)
  idx <- which(ds == 0, arr.ind = TRUE)
  e <- expand.grid(1:nrow(id), 1:nrow(idx))
  ide <- id[e[,1],]
  idxe <- idx[e[,2],]

  if (covgm) {
    x1 <- data[,,2]
    x2 <- data[,,1]
  } else {
    x1 <- c(data[cbind(ide[, 1], idxe[, 1])])
    x2 <- c(data[cbind(ide[, 2], idxe[, 2])])
  }

  # Baseline covariances for normalization
  c1 = cov(x1[cbind(ide[, 1], idxe[, 1])], x1[cbind(ide[, 2], idxe[, 2])])
  c2 = cov(x2[cbind(ide[, 1], idxe[, 1])], x2[cbind(ide[, 2], idxe[, 2])])

  # Loop over lag times to compute covariances at different spatial distances
  vgm <- lapply(lagstime, function(u) {
    id = cbind(wt_id - u, wt_id)
    diff = dates[wt_id] - dates[wt_id - u]
    id = id[diff == u,]

    cv <- sapply(1:length(dist), function(i) {
      d = dist[i]
      idx <- which(ds < d + 0.003 & ds > d - 0.003, arr.ind = TRUE)
      e <- expand.grid(1:nrow(id), 1:nrow(idx))
      ide <- id[e[,1],]
      idxe <- idx[e[,2],]

      if (covgm) {
        x1 <- data[,,2]
        x2 <- data[,,1]
      } else {
        x1 <- c(data[cbind(ide[, 1], idxe[, 1])])
        x2 <- c(data[cbind(ide[, 2], idxe[, 2])])
      }

      # Compute normalized covariance
      return(cov(x1[cbind(ide[, 1], idxe[, 1])], x2[cbind(ide[, 2], idxe[, 2])]) / sqrt(c1 * c2))
    })

    # Data frame with lag time, distances, and corresponding covariance values
    return(data.frame(lagtime = u, dist = dist, cov = cv))
  })

  # Combine and return results
  return(do.call(rbind, vgm))
}

generate_variable_index_pairs <- function(names) {
  # This function creates a matrix with all combinations of variable names
  u1 = sapply(names, function(v1) sapply(names, function(v2) v1))
  u2 = sapply(names, function(v1) sapply(names, function(v2) v2))
  ep = data.frame(cbind(u1[!upper.tri(u1)],u2[!upper.tri(u2)]))
  names(ep) = c("v1","v2")
  ep = rbind(ep[ep$v1==ep$v2,], ep[!ep$v1==ep$v2,])

  return(ep)
}
initialize_par_all_if_missing <- function(par_all, names, pairs, par_s) {
  # Initialize the `par_all` vector if it is missing, with default values or using `par_s`
  if (is.null(par_all)) {
    names_par_all <- c(paste(pairs, "dij", sep = ":"), "a", "b", "c", "d", "e", paste(names, "ci", sep = ":"),
                       paste(pairs, "aij", sep = ":"), paste(pairs, "bij", sep = ":"),
                       paste(pairs, "rij", sep = ":"), paste(pairs, "vij", sep = ":"),
                       paste(pairs, "ax", sep = ":"))

    par_all <- setNames(rep(1, length(names_par_all)), names_par_all)

    par_all[paste(pairs, "dij", sep = ":")] = 1
    par_all[paste(pairs, "bij", sep = ":")] = 1
    par_all[paste(names, "ci", sep = ":")] = 0
    par_all[paste(pairs[1:length(names)], "rij", sep = ":")] = par_s[1,]
    par_all[paste(pairs[1:length(names)], "vij", sep = ":")] = par_s[2,]
    par_all[paste(pairs, "ax", sep = ":")] = 0
    parms = c("a","b","c","d","e",paste(pairs, "aij", sep = ":"))
    par_all[parms] = c(0.1, 0.1, 0.1, 0.1,0.1,rep(1, length(pairs)))
  }
  return(par_all)
}
update_ax_parameters <- function(par_all, names, ax) {
  # Update the `ax` parameters in `par_all` based on the covariance information in `ax`
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
  colnames(a) = rownames(a) = names
  a = Matrix::nearPD(a)$mat
  for (v1 in names) {
    for(v2 in names){
      par_all[paste(paste(v1,v2,sep = "-"), "ax", sep = ":")] = a[v1,v2]
    }
  }
  return(par_all)
}
optimize_pairs_spatial <- function(par_all, data, names, Vi, uh, cr, max_it, ep) {
  pairs <- paste(ep[,1],ep[,2], sep = "-")

  # Optimize model parameters for each pair of variables using the log-likelihood function
  for (i in seq(nrow(ep))) {
    pair <- pairs[i]
    sp <- unlist(str_split(pair, "-"))
    if (sp[1] == sp[2]) {
      parms <- c(paste(pair, "rij", sep = ":"),
                 paste(pair, "vij", sep = ":"))
      par_all[parms] <- optim(par_all[parms], fn = loglik_pair, data = data, pair = pair, parms = parms,
                              par_all = par_all, ep = ep, names = names,
                              Vi = Vi, uh = uh, cr = cr,
                              control = list(maxit = max_it, trace = 2))$par
    } else {
      parms <- c(paste(pair, "ax", sep = ":"))
      par_all[parms] = par_all[parms]
    }

  }
  return(par_all)
}
optimize_pairs_spatiotemporal <- function(par_all, data, names, Vi, uh, cr, max_it, ep) {
  pairs <- paste(ep[,1],ep[,2], sep = "-")
  # Optimize model parameters for each pair of variables using the log-likelihood function
  for (i in seq(nrow(ep))) {
    pair <- pairs[i]
    sp <- unlist(str_split(pair, "-"))
    if (sp[1] == sp[2]) {
      parms = c(paste(sp[1], "ci", sep = ":"),paste(sp[2], "ci", sep = ":"),
                paste(pair, "bij", sep = ":"), paste(pair, "aij", sep = ":"),
                paste(pair, "rij", sep = ":"),paste(pair, "vij", sep = ":"))
      par_all[parms] <- optim(par_all[parms], fn = loglik_pair, data = data, pair = pair, parms = parms,
                              par_all = par_all, ep = ep, names = names,
                              Vi = Vi, uh = uh, cr = cr,
                              control = list(maxit = max_it, trace = 2))$par
    }
  }
  return(par_all)
}
optimize_temporal_parameters <- function(par_all, data, names, Vi, uh, cr, max_it, ep) {
  # Final optimization step for the subset of parameters across all variable pairs
  parms <- c("a", "b", "c", "d", "e")
  optimized_par <- optim(par_all[parms], fn = loglik, data = data, parms = parms,
                         par_all = par_all, ep = ep, names = names,
                         Vi = Vi, uh = uh, cr = cr,
                         control = list(maxit = max_it, trace = 2))$par
  par_all[parms] <- optimized_par
  return(par_all)
}


selectUniformPointsIndices <- function(lonlat, N) {
  library(stats)

  colnames(lonlat) <- c("lon", "lat")
  # K-means clustering
  clusters <- kmeans(lonlat, centers = N, nstart = 25)

  # For each cluster center, find the index of the closest original data point
  indices <- sapply(1:N, function(cluster_num) {
    subset_indices <- which(clusters$cluster == cluster_num)
    distances <- sqrt((lonlat$lon[subset_indices] - clusters$centers[cluster_num, "lon"])^2 +
                        (lonlat$lat[subset_indices] - clusters$centers[cluster_num, "lat"])^2)
    return(subset_indices[which.min(distances)])
  })

  return(indices)
}



selectPoints <- function(lonlat, betaIndex, v) {
  if (v > nrow(lonlat) || betaIndex > nrow(lonlat)) {
    stop("Invalid v or betaIndex")
  }

  calculateProbability <- function(alphaIndex) {
    1 / (sum((lonlat[alphaIndex,] - lonlat[betaIndex,])^2) + 1)
  }

  probabilities <- sapply(1:nrow(lonlat), calculateProbability)

  selectedIndices <- unique(c(betaIndex, sample(1:nrow(lonlat), size = v, prob = probabilities)))
  return(selectedIndices)
}

generate_spatial_index_pairs <- function(lonlat,n1, n2) {
  require(stringr)
  D = as.matrix(dist(lonlat))
  Ns = nrow(lonlat)
  rs <- selectUniformPointsIndices(lonlat, n1)
  Si <- sapply(rs, function(p1) {

    rs1 = selectPoints(lonlat, p1, v=n2)

    sapply(rs1, function(p2) {
      return(paste(min(c(p1, p2)), max(c(p1, p2)), sep = "-"))
    })
  })
  Si <- unique_elements(Si)
  Si <- matrix(unlist(lapply(str_split(Si, pattern = "-"), as.numeric)), byrow = T, nrow = length(Si))
  return(Si)
}

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

generate_variable_index_pairs <- function(names) {
  u1 = sapply(names, function(v1) sapply(names, function(v2) v1))
  u2 = sapply(names, function(v1) sapply(names, function(v2) v2))
  ep = data.frame(cbind(u1[!upper.tri(u1)],u2[!upper.tri(u2)]))
  names(ep) = c("v1","v2")
  ep = rbind(ep[ep$v1==ep$v2,], ep[!ep$v1==ep$v2,])
  return(ep)
}

unique_elements <- function(Si) {
  Si <- unlist(Si)
  Si <- Si[!duplicated(Si)]
  return(Si)
}
preprocess_data <- function(Ti, Si, lonlat) {
  e <- expand.grid(1:nrow(Ti), 1:nrow(Si))
  u <- Ti[e[, 1], 3]
  h <- ds(Si[e[, 2], 1], Si[e[, 2], 2],lonlat)
  uh <- cbind(u, h, Ti[e[, 1], 1], Ti[e[, 1], 2], Si[e[, 2], 1], Si[e[, 2], 2])
  return(list(uh = uh, u = u, h = h))
}
season_indices = function(dates, season, Year){

  years = unique(year(dates))
  years = c(years, max(years) + 1)
  md = lapply(years, function(Year) {
    if(Year == years[1]){
      is_leap = leap_year(Year)

      # Adjust min_day for non-leap years if necessary
      adjusted_min_day = ifelse(!is_leap && season$min_month == 2 && season$min_day == 29, 28, season$min_day)
      adjusted_max_day = ifelse(!is_leap && season$max_month == 2 && season$max_day == 29, 28, season$max_day)

      # Check if the season spans the end of the year
      if (season$min_month > season$max_month || (season$min_month == season$max_month && adjusted_min_day > adjusted_max_day)) {
        d1 = as.Date(paste(Year, min(month(dates[year(dates) %in% years[1]])), adjusted_min_day, sep = "-"))
        d2 = as.Date(paste(Year, season$max_month, adjusted_max_day, sep = "-"))
      }else{
        d1 = as.Date(paste(Year, season$min_month, adjusted_min_day, sep = "-"))
        d2 = as.Date(paste(Year, season$max_month, adjusted_max_day, sep = "-"))
      }

      return(seq(d1, d2, by = "day"))
    }else{
      if(Year == max(years)){
        Year = Year - 1
        is_leap = leap_year(Year)

        # Adjust min_day for non-leap years if necessary
        adjusted_min_day = ifelse(!is_leap && season$min_month == 2 && season$min_day == 29, 28, season$min_day)
        adjusted_max_day = ifelse(!is_leap && season$max_month == 2 && season$max_day == 29, 28, season$max_day)

        # Check if the season spans the end of the year
        if (season$min_month > season$max_month || (season$min_month == season$max_month && adjusted_min_day > adjusted_max_day)) {
          d1 = as.Date(paste(Year, season$min_month, adjusted_min_day, sep = "-"))
          d2 = max(dates)
          if(d1>d2){
            return(NULL)
          }else{
            return(seq(d1, d2, by = "day"))
          }
        }else{
          d1 = as.Date(paste(Year, season$min_month, adjusted_min_day, sep = "-"))
          d2 = as.Date(paste(Year, season$max_month, adjusted_max_day, sep = "-"))

        }

      }else{

        # Check for leap year
        is_leap = leap_year(Year)

        # Adjust min_day for non-leap years if necessary
        adjusted_min_day = ifelse(!is_leap && season$min_month == 2 && season$min_day == 29, 28, season$min_day)
        d1 = as.Date(paste(Year, season$min_month, adjusted_min_day, sep = "-"))

        # Adjust max_day for non-leap years if necessary
        if (season$min_month > season$max_month){
          is_leap = leap_year(Year)
          adjusted_max_day = ifelse(!is_leap && season$max_month == 2 && season$max_day == 29, 28, season$max_day)
          is_leap = leap_year(Year-1)
          adjusted_min_day = ifelse(!is_leap && season$min_month == 2 && season$min_day == 29, 28, season$min_day)
          d2 = as.Date(paste(Year, season$max_month, adjusted_max_day, sep = "-"))
          d1 = as.Date(paste(Year-1, season$min_month, adjusted_min_day, sep = "-"))
        }else{
          is_leap = leap_year(Year)
          adjusted_max_day = ifelse(!is_leap && season$max_month == 2 && season$max_day == 29, 28, season$max_day)
          d2 = as.Date(paste(Year, season$max_month, adjusted_max_day, sep = "-"))
        }

        return(seq(d1, d2, by = "day"))
      }
    }
  })
  md = do.call(c, md)
  md = format(md, "%Y-%m-%d CET")
  dates = format(dates, "%Y-%m-%d CET")
  if(missing(Year)){
    return(which(dates %in% md))
  }else{
    return(which(dates %in% md & year(dates)==Year))
  }
}

########### swg estimation ###############
filter_season_data <- function(data, dates, season, names) {
  # Ensure 'dates' is converted to Date class
  dates <- as.Date(dates)

  # Obtain indices for the specified season
  # Assuming 'season_indices' is a function that returns indices of dates within the specified season
  season_indices <- season_indices(dates, season)

  # Filter data and dates based on season indices
  data_filtered <- data[season_indices,,]
  dates_filtered <- dates[season_indices]

  # Further filter data to include only the specified variables (names)
  data_filtered <- data_filtered[,,names]

  # Order the filtered data by dates
  data_filtered <- data_filtered[order(dates_filtered),,, drop = FALSE]

  # Return the filtered and ordered data along with the corresponding dates
  return(list(data_filtered = data_filtered, dates_filtered = dates_filtered))
}
scale_data <- function(data, names, dates, window_size = 10) {
  # Initialize scale parameters storage
  scale_parm <- list(mu = list(), sd = list())

  # Convert dates to a consistent format for day extraction
  days <- substr(as.Date(dates), 6, 11)

  # Define a moving average function
  moving_average <- function(x, n) {
    len <- length(x)
    avg <- rep(NA, len)

    for (i in 1:len) {
      lower <- max(1, i - n %/% 2)
      upper <- min(len, i + n %/% 2)
      window <- x[lower:upper]
      avg[i] <- mean(window, na.rm = TRUE)
    }

    return(avg)
  }

  # Iterate over each variable except "Precipitation"
  for (v in names[!names %in% "Precipitation"]) {
    # Calculate daily mean and standard deviation for each variable
    mu <- sapply(unique(days), function(day) colMeans(data[days == day, , v], na.rm = TRUE))
    sd <- sapply(unique(days), function(day) apply(data[days == day, , v], 2, sd, na.rm = TRUE))

    # Smooth mu and sd using the moving average
    mu_smooth <- apply(mu, 1, function(x) moving_average(x, window_size))
    sd_smooth <- apply(sd, 1, function(x) moving_average(x, window_size))

    # Store smoothed parameters
    scale_parm$mu[[v]] <- mu_smooth
    scale_parm$sd[[v]] <- sd_smooth

    # Standardize data using smoothed parameters
    for (i in 1:nrow(data)) {
      day_index <- which(unique(days) == days[i])
      data[i, , v] <- (data[i, , v] - mu_smooth[day_index]) / sd_smooth[day_index]
    }
  }

  # Return the standardized data and scale parameters
  return(list(data = data, scale_parm = scale_parm))
}
estimate_lambda_transformations <- function(data, wt, names, lonlat) {
  ns = dim(data)[2]
  K = length(unique(wt))
  # Iterate over each weather type
  lambda_transformations <- lapply(1:K, function(k) {
    # For each variable
    variable_transforms <- lapply(names, function(v) {
      # For each spatial location
      location_transforms <- lapply(1:ns, function(j) {
        # Extract data for the current weather type, location, and variable
        x <- data[wt == k, j, v]
        # Calculate the proportion of zeros and its corresponding normal quantile
        proportion_zeros <- length(which(x == 0)) / length(x)
        q <- qnorm(proportion_zeros)

        # Apply the transformation, assuming orderNorm_all is defined elsewhere
        bx <- orderNorm_all(data = data[wt == k, , v], j = j, lonlat = lonlat, left = q)
        bx$q <- q  # Store the quantile in the result

        return(bx)
      })

      return(location_transforms)
    })
    return( variable_transforms)
  })
  threshold_precip = lapply(1:K, function(k){
    sapply(1:ns, function(j) lambda_transformations[[k]][[which(names=="Precipitation")]][[j]]$q)
  })
  return(list(lambda_transformations = lambda_transformations, threshold_precip = threshold_precip))
}
transformations <- function(data, wt, names, lonlat, lmbd) {
  ns = dim(data)[2]
  K = length(unique(wt))
  # Iterate over each weather type
  for (k in 1:K) {
    # Iterate over each variable
    for (v in names) {
      # Iterate over each spatial location
      for (j in 1:ns) {
        # Identify indices for the current weather type and non-zero values
        ind <- which(wt == k & data[, j, v] != 0)

        # Retrieve the lambda transformation parameters for the current combination
        m <- lmbd[[k]][[which(names == v)]][[j]]

        # Check if the second coefficient of the fit is NA; find an alternative if necessary
        if (is.na(m$fit$coef[2])) {
          # Find an index of a location with a valid second coefficient
          valid_index <- which(!sapply(1:ns, function(x) is.na(lmbd[[k]][[which(names == v)]][[x]]$fit$coefficients[2])))
          if (length(valid_index) > 0) {
            # Use the first valid transformation parameters found
            m <- lmbd[[k]][[which(names == v)]][[valid_index[1]]]
          } else {
            next  # Skip if no valid transformation parameters are found
          }
        }

        # Apply the transformation to the data
        # Assuming predict.orderNormTransf is defined elsewhere and applies the transformation
        if (length(ind) > 0) {  # Check if there are indices to update
          data[ind, j, v] <- predict.orderNormTransf(m, newdata = data[ind, j, v], inverse = FALSE)
        }
      }
    }
  }

  return(data)
}
estimate_gaussian_field_params <- function(data, wt, names, lonlat, tmax, max_it, n1, n2, dates, threshold_precip) {
  K = length(unique(wt))
  # Initialize the Gaussian field parameters storage
  gf_par <- vector(mode = "list", length = K)

  # Compute average pairwise correlations for each pair of variables across all locations
  cr <- sapply(names, function(v1) {
    sapply(names, function(v2) {
      mean(sapply(1:dim(data)[2], function(j) cor(data[, j, v1], data[, j, v2], use = "complete.obs")), na.rm = TRUE)
    })
  })
  colnames(cr) <- rownames(cr) <- names

  ep <- generate_variable_index_pairs(names)
  # Estimate spatial covariance structures for each pair of variables
  vgm <- lapply(1:nrow(ep), function(i) {
    variable <- unlist(ep[i, ])
    ds <- round(as.matrix(dist(lonlat)), 3)
    dist <- sort(unique(c(ds)))
    dist <- dist[seq(1, length(dist), length.out = 10)]
    vgm <- spacetime_cov(data = data[, , variable], wt_id = 2:dim(data)[1], locations = lonlat, ds = ds,
                         dates = dates, lagstime = 0, dist = dist, covgm = TRUE)
    vgm$v <- paste(variable[1], variable[2], sep = "-")
    vgm$v1 <- variable[1]
    vgm$v2 <- variable[2]

    return(vgm)
  })
  vgm <- do.call(rbind, vgm)

  # For each weather type, estimate Gaussian field parameters
  for (k in 1:K) {
    wt_id <- which(wt == k)
    wt_id <- wt_id[wt_id > tmax + 1]

    #Estimate Gaussian field parameters
    gf_par[[k]] <- estimation_gf(data = data, wt_id = wt_id, max_it = max_it, dates = dates,
                                 tmax = tmax, names = names, lonlat = lonlat, n1 = n1,
                                 n2 = n2, ax = vgm[vgm$lagtime==0&vgm$dist==max(vgm$dist),],
                                 cr = cr, threshold_precip = threshold_precip[[k]])$parm
  }

  return(gf_par)
}



######### Simulation #############

cov_matrices = function(par, lonlat, names, M) {
  # Function for generating covariance matrices for a multivariate space-time model.
  # The covariance is calculated using a specified covariance function, such as Gneiting's function.
  #
  # Arguments:
  #   par: Parameters for the covariance function, including spatial and temporal ranges, sill, and nugget among others.
  #   lonlat: A matrix or data frame of longitude and latitude coordinates for each spatial location.
  #   names: Names of the variables involved in the covariance calculation.
  #   M: The maximum time lag considered in the model.
  #
  # Returns:
  #   A list of covariance matrices for each time lag (up to M) and for each pair of variables.
  #   Each matrix represents the spatial covariance structure for a given time lag and variable pair.

  Nt = M + 1  # Number of time points considered
  Ns = nrow(lonlat)  # Number of spatial locations
  Nv = length(names)  # Number of variables

  # Generate all combinations of time points and spatial locations
  d = expand.grid(t1 = 1:Nt, t2 = 1:Nt, s1 = 1:Ns, s2 = 1:Ns)

  # Generate all combinations of variable pairs
  Vi = expand.grid(names, names)

  # Calculate time lags (u) and spatial distances (h) between all pairs of points
  u = d$t1 - d$t2
  h = ds(d$s1, d$s2, lonlat)  # Assuming ds is a function to calculate distances based on lonlat

  # Initialize a list to store covariance matrices
  cp = lapply(1:Nt, function(t1) {
    cp_v1 = lapply(names, function(v1) {
      cp_v2 = lapply(names, function(v2) {
        # Retrieve parameters for the current pair of variables and calculate covariance
        cov_params = par[(par$v1 == v1 & par$v2 == v2) | (par$v2 == v1 & par$v1 == v2), -c(1, 2)]
        dij = par$dij[(par$v1 == v1 & par$v2 == v2) | (par$v2 == v1 & par$v1 == v2)]
        cov = Gneiting(h, u, cov_params, dij)

        # Filter to the current time point and reshape the covariance values into a matrix
        up = (d$t1 == t1) & (d$t2 == 1)  # Assuming t2 is fixed for simplicity
        dd = d[up, ]
        co = cov[up]
        cv = matrix(0, ncol = Ns, nrow = Ns)
        for (i in 1:nrow(dd)) {
          cv[dd$s1[i], dd$s2[i]] = co[i]
        }
        return(cv)
      })
      return(do.call(rbind, cp_v2))
    })
    return(do.call(cbind, cp_v1))
  })
  return(cp)
}
calculate_Bk_matrices <- function(C_k_matrices) {
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
  Bk_matrices_solution <- t(v) %*% solve(A)

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
  Bk_0 <- t(chol(Bk_0_rhs))  # Perform Cholesky decomposition to obtain Bk_0

  Bk_matrices_list <- list(Bk_0 = Bk_0, bk = Bk_matrices_list)

  return(Bk_matrices_list)
}
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
  Z_initial <- simulate_Z(Bk = bk, M = AR_lag, num_steps = AR_lag+1, Z_initial = Z_initial, wt = wt)
  Z_initial <- as.matrix(do.call(cbind, Z_initial)[,(length(Z_initial)-AR_lag+1):length(Z_initial)])
  return(Z_initial)
}
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
apply_inverse_transformations <- function(sim, wt, parm, names) {
  # Reverses transformations applied to the simulated data for each weather type and variable.
  #
  # Arguments:
  #   sim: A 3D array of simulated weather data, with dimensions [time, locations, variables].
  #   wt: A vector of simulated weather types for each time step.
  #   parm: Parameters object containing lambda transformations and other model parameters.
  #   names: Names of the variables included in the simulation.
  #
  # Returns:
  #   The sim array with inverse transformations applied, bringing the data back to its original scale.


  lmbd = parm$lmbd
  nt = dim(sim)[1]
  ns = dim(sim)[2]
  for(k in unique(wt)){
    for(v in names) {
      for(j in 1:ns){
        if(v == "Precipitation"){
          q = lmbd[[k]][[which(names==v)]][[j]]$q
        }else{
          q = -Inf
        }
        m = lmbd[[k]][[which(names==v)]][[j]]
        if(is.na(m$fit$coef[2])){
          m = lmbd[[k]][[which(names==v)]][[which(!is.na(sapply(1:ns, function(j) lmbd[[k]][[which(names==v)]][[j]]$fit$coefficients[[2]])))[1]]]
        }
        ind = wt==k & sim[,j,v]  > q
        sim[wt==k & sim[,j,v]  <= q,j,v] = 0
        if(!length(which(ind))==0){
          sim[ind,j,v] = inv_orderNorm_Transf(m ,new_points_x_t=sim[ind,j,v], left =q)
        }
      }
    }
  }

  return(sim)
}
rescale_data <- function(sim, parm, names, dates) {
  # Rescales simulated weather data for variables other than "Precipitation" to their original scale.
  #
  # Arguments:
  #   sim: A 3D array of simulated weather data, with dimensions [time, locations, variables].
  #   parm: An object containing model parameters, including scale parameters (mu and sd) for each variable.
  #   names: Names of the variables included in the simulation.
  #   dates: Vector of dates corresponding to the time dimension in the sim array.
  #
  # Returns:
  #   The sim array with data rescaled back to its original scale for variables other than "Precipitation".

  nt <- dim(sim)[1]  # Number of time steps
  ns <- dim(sim)[2]  # Number of spatial locations
  # Extract day identifiers from dates to match with scale parameters
  days <- substr(dates, 6, 11)
  dayso <- substr(parm$dates, 6, 11)
  for (v in names[!names=="Precipitation"]) {
    mu = parm$scale_parm$mu[[v]]
    sd = parm$scale_parm$sd[[v]]
    for (i in 1:nt) {
      sim[i,,v] = (sim[i,,v]*sd[which(unique(dayso)==days[i]),])+mu[which(unique(dayso)==days[i]),]
    }
  }

  return(sim)
}

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

