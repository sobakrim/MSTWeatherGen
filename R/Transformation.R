#' Ordered Normalization
#'
#' Function for ordered normalization using Ordered Quantile Normalizing transformation.
#'
#' @param x Numeric vector to be normalized.
#' @param left Numeric value specifying the left truncation.
#' @param n_logit_fit The number of points to use in the logistic regression for extrapolation, default is the lesser of the length of x or 100,000.
#' @param ... Additional arguments passed to the glm function.
#' @param warn Logical flag indicating whether to issue a warning when ties are detected in the data.
#'
#' @return A list with class 'orderNorm' containing the normalized data (x.t), the original data (x), and a fitted model object (fit) for future extrapolation.
#'
#' @importFrom crch qcnorm
#' @importFrom stats glm
#'
#' @keywords internal
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
    ties_status <- 1
  }
  if(length(unique(x))==1) x = x + rnorm(length(x), sd = 1e-7)
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
    stats::glm(q_red ~ x_red , family = 'binomial', 
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
#' Predict Binomial
#'
#' Function to predict binomial probabilities using a fitted model object.
#'
#' @param fit Fitted model object containing coefficients.
#' @param newdata New data for prediction.
#'
#' @return Predicted binomial probabilities.
#'
#' @keywords internal
predict_binomial = function(fit, newdata){
  if(missing(newdata)) newdata = fit$x_red
  pred = fit$coef[1]+fit$coef[2]*newdata
  return(1/(1+exp(-pred)))
}
#' Order Normalization for All Variables
#'
#' Function to perform normalization on a selected variable 'j' from 'data', considering spatial information 
#' from 'coordinates' and a threshold 'left' for the transformation. This is particularly useful for variables with 
#' many zero values, ensuring robust normalization across spatially correlated variables.
#'
#' @param data A matrix or data frame containing the variables to be normalized. Rows represent observations, 
#'             and columns represent different variables.
#' @param j The index of the variable within 'data' to be normalized.
#' @param coordinates A matrix or data frame of spatial coordinates corresponding to each 
#'                    variable in 'data'. It is used to calculate spatial distances and determine proximity.
#' @param left A parameter for the orderNorm function, specifying the transformation threshold.
#'
#' @return The result of applying the orderNorm function to the selected variable 'j', considering spatial proximity 
#'         and handling variables with a high proportion of zero values.
#'
#' @keywords internal
orderNorm_all <- function(data, j, coordinates, left) {
  # Function to perform normalization on a selected variable 'j' from 'data', considering spatial information 
  # from 'coordinates' and a threshold 'left' for the transformation. This is particularly useful for variables with 
  # many zero values, ensuring robust normalization across spatially correlated variables.
  
  # Arguments:
  #   data: A matrix or data frame containing the variables to be normalized. Rows represent observations, 
  #         and columns represent different variables.
  #   j: The index of the variable within 'data' to be normalized.
  #   coordinates: A matrix or data frame of spatial coordinates (coordinates) corresponding to each 
  #           variable in 'data'. It is used to calculate spatial distances and determine proximity.
  #   left: A parameter for the orderNorm function, specifying the transformation threshold.
  
  # Returns:
  #   The result of applying the orderNorm function to the selected variable 'j', considering spatial proximity 
  #   and handling variables with a high proportion of zero values.
  
  D = as.matrix(dist(coordinates))
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
#' Predict Method for orderNormTransf Objects
#'
#' Function to predict the normalized values for new data using a fitted orderNormTransf object or 
#' to perform inverse transformation.
#'
#' @param object Fitted orderNormTransf object.
#' @param newdata Numeric vector of new data to be transformed. If NULL, transformation is applied to the original data.
#' @param inverse Logical indicating whether to perform inverse transformation. Default is FALSE.
#' @param warn Logical indicating whether to issue a warning when ties are detected in the data. Default is TRUE.
#' @param ... Additional arguments to be passed.
#'
#' @return Transformed or inverse-transformed data based on the specified operation.
#'
#' @keywords internal
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
    
    newdata[!na_idx] <- orderNormTransf(orderNorm_obj = object, new_points = newdata[!na_idx],
                                        warn = warn, left = object$q)
    return(newdata)
  } 
  
  # Perform reverse transformation
  if (is.null(newdata)) newdata <- object$x.t
  
  na_idx <- is.na(newdata)
  newdata[!na_idx] <- inv_orderNorm_Transf(object, newdata[!na_idx], warn)
  
  return(newdata)
}
#' Inverse Transformation for orderNorm Objects
#'
#' Function to perform inverse transformation or normalization applied by the orderNorm function.
#'
#' @param orderNorm_obj Fitted orderNorm object containing details of the original normalization or transformation.
#' @param new_points_x_t Transformed data points for which the original values are to be estimated.
#' @param left Parameter used in the original transformation, affecting the normalization method.
#' @param warn Logical indicating whether to issue warnings for transformations outside the observed domain. Default is FALSE.
#'
#' @return Estimated original data values corresponding to the transformed data points.
#'
#' @keywords internal
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
      stats::approx(x_t, old_points, xout = new_points_x_t, rule = 2:1)
    )
  }else{
    vals <- suppressWarnings(
      stats::approx(x_t, old_points, xout = new_points_x_t, rule = 1)
    )
  }
  # If predictions have been made outside observed domain
  if (any(is.na(vals$y))) {
    
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
#' Transforms New Data Points using orderNorm Object
#'
#' Function to transform new data points based on a previously fitted normalization or transformation model.
#'
#' @param orderNorm_obj Fitted orderNorm object containing details of the original normalization or transformation.
#' @param new_points New data points to be transformed.
#' @param warn Logical indicating whether to issue warnings for transformations outside the observed domain.
#' @param left Parameter affecting the normalization method used in the original transformation.
#'
#' @return Transformed data points based on the fitted model.
#'
#' @keywords internal
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
    stats::approx(old_points, x_t, xout = new_points, rule = 1)
  )
  
  # If predictions have been made outside observed domain
  if (any(is.na(vals$y))) {
    fit <- orderNorm_obj$fit
    p <- crch::qcnorm(predict_binomial(fit), left = left)
    l_idx <- vals$x < min(old_points, na.rm = TRUE)
    h_idx <- vals$x > max(old_points, na.rm = TRUE)
    
    # Check 
    if (any(l_idx)) {
      xx <- data.frame(x_red = vals$x[l_idx])
      q <- predict_binomial(fit, newdata = xx)$x_red 
      vals$y[l_idx] <- crch::qcnorm(q, left = left) - 
        (min(p, na.rm = TRUE) - min(x_t, na.rm = TRUE))
      
    }
    if (any(h_idx)) {
      xx <- data.frame(x_red = vals$x[h_idx])
      q <- predict_binomial(fit, newdata = xx)$x_red 
      vals$y[h_idx] <- crch::qcnorm(q, left = left) - 
        (max(p, na.rm = TRUE) - max(x_t, na.rm = TRUE))
    }
  }
  
  return(vals$y)
}

#' Scale Data
#'
#' Function to scale the input data based on daily mean and standard deviation, and optionally smooth them using a moving average.
#'
#' @param data Input weather data, organized as a 3D array where the dimensions are [time, locations, variables].
#' @param names Names of the variables in the data.
#' @param dates Vector of dates corresponding to the time dimension in the data.
#' @param window_size Size of the window for the moving average used to smooth the daily mean and standard deviation. Default is 30.
#'
#' @return A list containing the scaled data and scale parameters (mean and standard deviation) for each variable.
#' @importFrom lubridate year
#' @keywords internal
scale_data <- function(data, names, dates, window_size = 30) {
  
  if(length(unique(lubridate::year(dates))) > 5){
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
        data[i, , v] <- (data[i, , v] - mu_smooth[day_index,]) / sd_smooth[day_index,]
      }
    }
  }else{
    warning('No scale is used. The number of years considered is too small to estimate the seasonal means and standard deviations.')
    scale_parm = NULL
  }
  # Return the standardized data and scale parameters
  return(list(data = data, scale_parm = scale_parm))
}
#' Estimate Lambda Transformations
#'
#' Function to estimate lambda transformations for each variable, location, and weather type based on the data.
#'
#' @param data Input weather data, organized as a 3D array where the dimensions are [time, locations, variables].
#' @param wt Vector indicating the weather type for each observation.
#' @param names Names of the variables in the data.
#' @param coordinates Spatial coordinates corresponding to each location in the data.
#'
#' @return A list containing lambda transformations for each variable, location, and weather type, along with the thresholds for precipitation.
#'
#' @keywords internal
estimate_lambda_transformations <- function(data, wt, names, coordinates) {
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
        if(v=="Precipitation"){
          proportion_zeros <- length(which(x == 0)) / length(x)
          q <- qnorm(proportion_zeros)
        }else{
          q <- -Inf
        }
        
        # Apply the transformation, assuming orderNorm_all is defined elsewhere
        bx <- orderNorm_all(data = data[wt == k, , v], j = j, coordinates = coordinates, left = q)
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
#' Transform Data using Lambda Transformations
#'
#' Function to transform data using lambda transformations based on the provided lambda parameters.
#'
#' @param data Input weather data, organized as a 3D array where the dimensions are [time, locations, variables].
#' @param wt Vector indicating the weather type for each observation.
#' @param names Names of the variables in the data.
#' @param coordinates Spatial coordinates corresponding to each location in the data.
#' @param lmbd Lambda transformation parameters estimated for each variable, location, and weather type.
#'
#' @return Transformed data using lambda transformations.
#'
#' @keywords internal
transformations <- function(data, wt, names, coordinates, lmbd) {
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
        i = 1
        while(is.na(m$fit$coef[2])) {
          # Find an index of a location with a valid second coefficient
          valid_indices <- which(!sapply(1:ns, function(x) is.na(lmbd[[k]][[which(names == v)]][[x]]$fit$coef[2])))
          valid_indices <- order(apply(coordinates[valid_indices,], 1, function(point) {
            sqrt(sum((point - coordinates[j,])^2))
          }))
          if (length(valid_indices) > 0) {
            # Use the first valid transformation parameters found
            m$fit <- lmbd[[k]][[which(names == v)]][[valid_indices[i]]]$fit
          } else {
            next  # Skip if no valid transformation parameters are found
          }
          i = i+1
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


#' Apply Inverse Transformations
#'
#' Function to apply inverse transformations to the simulated data, reversing the lambda transformations applied during simulation.
#'
#' @param sim Simulated weather data, organized as a 3D array where the dimensions are [time, locations, variables].
#' @param wt Vector indicating the weather type for each time step.
#' @param parm Parameters object containing lambda transformations and other model parameters.
#' @param names Names of the variables included in the simulation.
#'
#' @return Simulated data with inverse transformations applied to bring it back to its original scale.
#'
#' @keywords internal

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
#' Rescale Simulated Data
#'
#' Function to rescale simulated weather data for variables other than "Precipitation" to their original scale.
#'
#' @param sim Simulated weather data, organized as a 3D array where the dimensions are [time, locations, variables].
#' @param parm Parameters object containing model parameters, including scale parameters (mu and sd) for each variable.
#' @param names Names of the variables included in the simulation.
#' @param dates Vector of dates corresponding to the time dimension in the sim array.
#'
#' @return Simulated data with variables other than "Precipitation" rescaled back to their original scale.
#'
#' @keywords internal
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

