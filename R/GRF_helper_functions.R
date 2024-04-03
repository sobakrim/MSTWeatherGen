generate_variable_index_pairs <- function(names) {
  # This function creates a matrix with all combinations of variable names
  u1 = sapply(names, function(v1) sapply(names, function(v2) v1))
  u2 = sapply(names, function(v1) sapply(names, function(v2) v2))
  ep = data.frame(cbind(u1[!upper.tri(u1)],u2[!upper.tri(u2)]))
  names(ep) = c("v1","v2")
  ep = rbind(ep[ep$v1==ep$v2,], ep[!ep$v1==ep$v2,])
  
  return(ep)
}
initialize_par_all_if_missing <- function(par_all, names, pairs, par_s, ax, cr) {
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
    par_all[paste(pairs[1:length(names)], "rij", sep = ":")] <- par_s[1,] 
    par_all[paste(pairs[1:length(names)], "vij", sep = ":")] <- par_s[2,] 
    par_all[paste(pairs, "ax", sep = ":")] <- 0
    parms <- c("a","b","c","d","e",paste(pairs, "aij", sep = ":"))
    par_all[parms] <- c(0.1, 0.1, 0.1, 0.1,0.1,rep(1, length(pairs)))
  }
  
  # Update ax parameters based on covariance information
  par_all <- update_ax_parameters(par_all, names, ax)
  
  parm <- param(par_all, names)
  beta <- try(compute_beta(parm, names, cr), silent = T)
  ch <- try(chol(beta), silent = T)
  if(is.character(ch)){
    par_s <- matrix(rep(1, length(names)^2), ncol = length(names), nrow= length(names))
    par_all[paste(pairs[1:length(names)], "rij", sep = ":")] <- par_s[1,] 
    par_all[paste(pairs[1:length(names)], "vij", sep = ":")] <- par_s[2,] 
    par_all <- update_ax_parameters(par_all, names, ax)
  } 
  return(par_all)
}
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
  }
  colnames(ax) = rownames(ax) = names
  for (v1 in names) {
    for(v2 in names){
      par_all[paste(paste(v1,v2,sep = "-"), "ax", sep = ":")] = ax[v1,v2]
    }
  }
  return(par_all)
}
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
#' @importFrom stringr str_split
optimize_pairs_spatial <- function(par_all, data, names, Vi, uh, cr, max_it, ep) {
  pairs <- paste(ep[,1],ep[,2], sep = "-")
  
  # Optimize model parameters for each pair of variables using the log-likelihood function
  for (i in seq(nrow(ep))) {
    pair <- pairs[i]
    sp <- unlist(stringr::str_split(pair, "-"))
    if (sp[1] == sp[2]) {
      parms <- c(paste(pair, "rij", sep = ":"),paste(pair, "ax", sep = ":"),
                 paste(pair, "vij", sep = ":"))
      par_all[parms] <- optim(par_all[parms], fn = loglik_pair, data=data, pair = pair, parms = parms,
                              par_all = par_all,ep = ep,
                              names = names, 
                              Vi = Vi, uh = uh[uh[,1]==0,], cr = cr,
                              control = list(maxit = max_it))$par
    } else {
      parms <- c(paste(pair, "ax", sep = ":"))
      par_all[parms] <- optim(par_all[parms], fn = loglik_pair, data=data, pair = pair, parms = parms,
                              par_all = par_all,ep = ep,
                              names = names, 
                              Vi = Vi, uh = uh[uh[,1]==0,], cr = cr,
                              control = list(maxit = max_it))$par   
    }
    
  }
  return(par_all)
}
#' @importFrom stringr str_split
optimize_pairs_spatiotemporal <- function(par_all, data, names, Vi, uh, cr, max_it, ep) {
  pairs <- paste(ep[,1],ep[,2], sep = "-")
  # Optimize model parameters for each pair of variables using the log-likelihood function
  for (i in seq(nrow(ep))) {
    pair <- pairs[i]
    sp <- unlist(stringr::str_split(pair, "-"))
    if (sp[1] == sp[2]) {
      parms <- c(paste(sp[1], "ci", sep = ":"),paste(sp[2], "ci", sep = ":"),paste(pair, "ax", sep = ":"),
                paste(pair, "bij", sep = ":"), paste(pair, "aij", sep = ":"),
                paste(pair, "rij", sep = ":"),paste(pair, "vij", sep = ":"))
      par_all[parms] <- optim(par_all[parms], fn = loglik_pair, data = data, pair = pair, parms = parms,
                              par_all = par_all, ep = ep, names = names, 
                              Vi = Vi, uh = uh, cr = cr,
                              control = list(maxit = max_it))$par
    }else{
      pair <- paste(ep[i,1],ep[i,1], sep = "-")
      parms <- c(paste(sp[1], "ci", sep = ":"),paste(sp[2], "ci", sep = ":"),
                paste(pair, "bij", sep = ":"), paste(pair, "aij", sep = ":"))
      pair <- paste(ep[i,2],ep[i,2], sep = "-")
      parms <- c(parms, paste(pair, "bij", sep = ":"), paste(pair, "aij", sep = ":"))
      par_all[parms] <- optim(par_all[parms], fn = loglik_pair, data = data, pair = pairs[i], parms = parms,
                              par_all = par_all, ep = ep, names = names, 
                              Vi = Vi, uh = uh, cr = cr,
                              control = list(maxit = max_it))$par
      
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
                         control = list(maxit = max_it))$par
  par_all[parms] <- optimized_par
  return(par_all)
}

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
  Si <- generate_spatial_index_pairs(coordinates, n1, n2)
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

  
  for (v in 1:3) {
    # Optimize spatial parameters for each pair
    par_all <- optimize_pairs_spatial(par_all, data, names, Vi, uh, cr, max_it, ep)
    
    # Optimize temporal parameters
    par_all <- optimize_temporal_parameters(par_all, data, names, Vi, uh, cr, max_it, ep)

    # Optimize spatiotemporal parameters for each pair
    par_all <- optimize_pairs_spatiotemporal(par_all, data, names, Vi, uh, cr, max_it, ep)
  }
  
  # Construct parameter and beta matrices
  parm <- param(par_all, names)
  beta <- compute_beta(parm, names, cr)
  beta <- sapply(1:nrow(ep), function(i) beta[ep[i,1], ep[i,2]])
  par_all[1:length(beta)] <- beta
  
  return(list(parm = param(par_all,names), par_all = par_all))
}

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
preprocess_data <- function(Ti, Si, coordinates) {
  e <- expand.grid(1:nrow(Ti), 1:nrow(Si))
  u <- Ti[e[, 1], 3] 
  h <- ds(Si[e[, 2], 1], Si[e[, 2], 2],coordinates)
  uh <- cbind(u, h, Ti[e[, 1], 1], Ti[e[, 1], 2], Si[e[, 2], 1], Si[e[, 2], 2])
  return(list(uh = uh, u = u, h = h))
}

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
