#' Matern Covariance Function
#'
#' Calculates the Matern covariance function for a given vector of distances. This function is intended for internal use within the package to compute spatial covariances with the Matern model. It is not exported for end-user interaction.
#'
#' @param h Numeric vector of distances between points.
#' @param r Scalar range parameter of the Matern function, affecting the spatial correlation decay.
#' @param v Scalar smoothness parameter of the Matern function, controlling the smoothness of the resulting field.
#'
#' @return Numeric vector representing the covariance values calculated using the Matern function for the distances in `h`.
#'
#' @keywords internal
Matern <- function(h, r, v) {
  # Calculates the Matern covariance function for a given vector of distances.
  
  # Arguments:
  #   h: Numeric vector of distances between points.
  #   r: Scalar range parameter of the Matern function, affecting spatial correlation decay.
  #   v: Scalar smoothness parameter of the Matern function, controlling the smoothness of the resulting field.
  
  # Returns:
  #   Numeric vector representing the covariance values calculated using the Matern function for the distances in h.
  
  rt <- (2 ^ (1 - v)) / gamma(v) * ((r * abs(h)) ^ v) * besselK(r * abs(h), nu = v)
  rt[h == 0] <- 1  # Ensures that the covariance at distance 0 is 1.
  return(rt)
}
#' Gneiting's Spatio-Temporal Covariance Model
#'
#' Computes the covariance based on Gneiting's spatio-temporal model. Intended for internal package use.
#'
#' @param h Numeric vector of spatial distances.
#' @param u Numeric vector of temporal distances.
#' @param par Numeric vector of parameters for the covariance function.
#' @param dij Correlation parameter between variables i and j.
#'
#' @return Numeric vector of covariance values calculated using Gneiting's model.
#'
#' @keywords internal

Gneiting <- function(h, u, par, dij) {
  # Multivariate space-time Gneiting's covariance function
  # From the paper: https://doi.org/10.1016/j.spasta.2022.100706
  
  # Arguments:loglik
  #   h: Numeric vector of spatial distances.
  #   u: Numeric vector of temporal distances.
  #   par: Numeric vector of parameters used in the covariance function.
  #   dij: correlation parameter between variable i and j.
  
  # Returns:
  #   Covariance value(s) calculated using Gneiting's spatio-temporal covariance model.
  
  if(!is.numeric(par)) par <- as.numeric(par)
  
  # Unpack parameters from the 'par' vector for clarity.
  a1 <- par[1]
  d1 <- par[2]
  g1 <- par[3]
  a2 <- par[4]
  d2 <- par[5]
  g2 <- par[6]
  b1 <- par[7]
  e1 <- par[8]
  l1 <- par[9]
  b2 <- par[10]
  e2 <- par[11]
  l2 <- par[12]
  c <- par[13]
  f <- par[14]
  m <- par[15]
  ai <- par[16]
  aj <- par[17]
  bi <- par[18]
  bj <- par[19]
  ci <- par[20]
  cj <- par[21]
  rii <- par[22]
  rjj <- par[23]
  vii <- par[24]
  vjj <- par[25]
  ax <- par[26]
  
  # Calculated intermediate parameters for the covariance calculation.
  vij <- (vii + vjj) / 2
  rij <- sqrt((rii^2 + rjj^2) / 2)
  
  eij <- dij * ((rii^vii * rjj^vjj) / rij^(2*vij)) * 
    (gamma(vij) / (gamma(vii)^(1/2) * gamma(vjj)^(1/2))) * 
    sqrt((1-ai^2)*(1-aj^2)) * sqrt((1-bi^2)*(1-bj^2))
  
  muij <- (((a1 * abs(u))^(2*b1) + 1)^(c) - (ai*aj * ((a2 * abs(u))^(2*b2) + 1)^(-c))) 
  rhoij <- 1 / (((d1 * abs(u))^(2*e1) + 1)^(f) - (bi*bj * ((d2 * abs(u))^(2*e2) + 1)^(-f)))
  
  A1 <- eij / muij
  A2 <- Matern(abs(h), r = (rij^2 / muij)^(1/2), v = vij) * rhoij
  
  A3 <- ax * 1 / (((g1 * abs(u))^(2*l1) + 1)^(m) - ci*cj * ((g2 * abs(u))^(2*l2) + 1)^(-m))
  
  return(A1 * A2 + A3)
}

#' Construct Covariance Parameters DataFrame
#'
#' Creates a data frame of covariance parameters for all possible pairs of variables. This function is 
#' designed for internal use, facilitating the organization of parameters for spatio-temporal modeling.
#'
#' @param par Named vector of parameters.
#' @param names Character vector of variable names.
#'
#' @return A data frame where each row corresponds to a pair of variables (including self-pairs) 
#' and their associated spatio-temporal covariance parameters.
#'
#' @keywords internal

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
  u$a1 <- par["a1"]
  u$d1 <- par["d1"]
  u$g1 <- par["g1"]
  u$a2 <- par["a2"]
  u$d2 <- par["d2"]
  u$g2 <- par["g2"]
  u$b1 <- par["b1"]
  u$e1 <- par["e1"]
  u$l1 <- par["l1"]
  u$b2 <- par["b2"]
  u$e2 <- par["e2"]
  u$l2 <- par["l2"]
  u$c <- par["c"]
  u$f <- par["f"]
  u$m <- par["m"]
  # Loop through each pair to populate the data frame with corresponding parameter values
  for (i in seq_len(J)) {
    # Extract and assign specific parameters for each pair based on naming convention
    u$ai[i] <- par[paste(u$v1[i],  "ai", sep = ":")]
    u$aj[i] <- par[paste(u$v2[i],  "ai", sep = ":")]
    
    u$bi[i] <- par[paste(u$v1[i],  "bi", sep = ":")]
    u$bj[i] <- par[paste(u$v2[i],  "bi", sep = ":")]
    
    u$ci[i] <- par[paste(u$v1[i],  "ci", sep = ":")]
    u$cj[i] <- par[paste(u$v2[i],  "ci", sep = ":")]
    
    u$rii[i] <- par[paste(paste(u$v1[i],u$v1[i],  sep = "-"),  "rij", sep = ":")]
    u$rjj[i] <- par[paste(paste(u$v2[i],u$v2[i],  sep = "-"),  "rij", sep = ":")]
    
    u$vii[i] <- par[paste(paste(u$v1[i],u$v1[i],  sep = "-"),  "vij", sep = ":")]
    u$vjj[i] <- par[paste(paste(u$v2[i],u$v2[i],  sep = "-"),  "vij", sep = ":")]
    
    u$ax[i] <- par[paste(paste(u$v1[i],u$v2[i],  sep = "-"),  "ax", sep = ":")]
    
    u$dij[i] <- par[paste(paste(u$v1[i],u$v2[i],  sep = "-"),  "dij", sep = ":")]
  }
  return(u)
}
#' Compute Beta Correlations
#'
#' This function calculates the beta correlation coefficients between variables based on the Gneiting function, adjusted for a correction term. It is intended for internal use within package functions to adjust initial correlation values using specified parameters.
#'
#' @param parm A data frame or list containing parameters for the Gneiting function.
#' @param names Character vector of variable names to calculate correlations between.
#' @param cr Matrix of initial correlation values between the variables.
#'
#' @return Symmetric matrix of adjusted correlation coefficients (beta) between the variables.
#'
#' @importFrom Matrix nearPD
#' @keywords internal

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
  ax = sapply(names, function(v1){
    sapply(names, function(v2){
      par <- get_parameters(v1, v2)
      if(any(is.na(par))) par <- get_parameters(v2, v1)
      return(par[26]/ (1 - par[20] * par[21]))
    })
  })
  cr = Matrix::nearPD(cr-ax)$mat
  # Iterate over pairs of variables to compute and assign the beta values
  for (j in seq_along(names)) {
    for (k in seq(j, J)) {
      v1 <- names[j]  # First variable in the pair
      v2 <- names[k]  # Second variable in the pair
      par <- get_parameters(v1, v2)  # Retrieve parameters for the current pair
      
      # Calculate the correlation coefficient using the Gneiting function and correction term
      cc = Gneiting(0, 0, par, dij = 1)  # Gneiting function calculation for the pair
      ax = par[26] / (1 - par[20] * par[21]) # Correction term calculation
      beta_val <- (cr[v1, v2]) / (cc - ax)  # Adjusted correlation coefficient
      
      beta[v1, v2] <- beta[v2, v1] <- beta_val  # Symmetric assignment to ensure the matrix is symmetric
    }
  }
  
  return(beta)  
}
#' Extract Correction Terms Matrix
#'
#' Extracts a matrix of correction terms ('ax') for each pair of variables based on the model parameters. Designed for internal use to facilitate calculations involving correction terms in spatial or spatio-temporal modeling.
#'
#' @param parm A data frame or list containing the model parameters, including 'ax' values.
#' @param names Character vector specifying the variable names for which correction terms are to be calculated.
#'
#' @return A square matrix where each element [i, j] represents the correction term ('ax') between the ith and jth variables, facilitating the adjustment of correlations or covariances between them.
#'
#' @keywords internal


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
      ax = parm$ax[parm$v1==v1&parm$v2==v2|parm$v1==v2&parm$v2==v1]
      return(ax)
    })
  })
  return(ax)
}
#' Extract Beta Coefficients Matrix
#'
#' Extracts a matrix of beta coefficients ('dij') for each pair of variables from the provided model parameters. Intended for internal use, this function supports spatial and spatio-temporal modeling by organizing pairwise beta coefficients into a structured format.
#'
#' @param parm A data frame or list containing the model parameters, which must include 'dij' values representing beta coefficients between pairs of variables.
#' @param names Character vector of variable names for which beta coefficients are to be extracted.
#'
#' @return A square matrix where each element [i, j] contains the beta coefficient ('dij') between the ith and jth variables. This matrix is crucial for modeling the interactions and dependencies between different variables in the model.
#'
#' @keywords internal

extract_beta <- function(parm, names) {
  
  ax = sapply(names, function(v1){
    sapply(names, function(v2){
      parm$dij[parm$v1==v1&parm$v2==v2|parm$v1==v2&parm$v2==v1]
    })
  })
  return(ax)
}
#' Compute Log-likelihood for Variable Pair
#'
#' Calculates the log-likelihood for a given pair of variables using the Gneiting spatio-temporal covariance model. This function is part of the internal mechanism for optimizing model parameters based on observed data.
#'
#' @param par Current parameters being optimized.
#' @param parms Indices or names of parameters in 'par' to be updated.
#' @param pair A string indicating the pair of variables (e.g., "temperature-wind") being analyzed.
#' @param par_all Complete set of parameters for the model.
#' @param data 3D array of observed data, with dimensions corresponding to times, locations, and variables.
#' @param names Vector of variable names, indicating the variables' names (e.g., "temperature" and "wind").
#' @param Vi Matrix where each line corresponds for a possible combination of variables in "names".
#' @param h Vector of spatial distances for the pair.
#' @param u Vector of temporal distances for the pair.
#' @param uh Matrix containing pairs of spatial and temporal distances, and additional information.
#' @param ep A matrix or data frame defining pairs of variables.
#' @param cr Correlation matrix, initial or base correlations between variables.
#'
#' @return The log-likelihood value for the given pair of variables based on the current model parameters.
#'
#' @importFrom VGAM pbinorm
#' @keywords internal

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
  
  J = length(names)  # Number of variables
  pairs = paste(ep[,1], ep[,2], sep = "-")  # Constructing pairs from 'ep' data frame
  
  par_all[parms] = par  # Update specific parameters in the complete set
  
  par = par_all  # Use the updated parameter set for computations
  sp = unlist(strsplit(pair, "-"))  # Split the pair string to individual variables
  v = which(Vi[,1] == sp[1] & Vi[,2] == sp[2])  # Find the index of the pair in 'Vi'
  
  # Update and compute model parameters
  parm <- param(par, names)
  ax <- Matrix::nearPD(compute_ax(parm, names))$mat  # Compute ax correction terms
  beta <- try(compute_beta(parm, names, cr), silent = T)  # Compute beta coefficients
  
  # Attempt Cholesky decompositions for 'ax' and 'beta', checking for positive definiteness
  ae <- try(chol(ax), silent = TRUE)
  be <- try(chol(beta), silent = TRUE)
  
  if (!is.character(be) & (!is.character(ae) )) {
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
    if (any(par[c(1:21,23:24)] < 0) | any(par[c(2,3,5,7:13)] > 1)) {
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

#' Total Log-Likelihood Calculation
#'
#' Calculates the total log-likelihood for spatial or spatio-temporal data across all variable pairs. Utilizes the Gneiting spatio-temporal covariance model to integrate log-likelihood contributions from each variable pair. This function is core to the optimization process within model fitting.
#'
#' @param par Vector of parameter estimates currently being optimized.
#' @param parms Indices or names of parameters within 'par' that are subject to update.
#' @param par_all Comprehensive list of all model parameters.
#' @param data 3D array of observed data across locations, times, and variables.
#' @param names Character vector of variable names.
#' @param Vi Matrix indicating all combinations of variables for analysis.
#' @param h Spatial distances vector for variable pairs.
#' @param u Temporal distances vector for variable pairs.
#' @param uh Combined matrix of spatial and temporal distances with additional identifiers.
#' @param ep Data frame defining variable pairs for analysis.
#' @param cr Initial correlation matrix across variables.
#'
#' @return Total log-likelihood value for the observed data given the current model parameters.
#'
#' @importFrom VGAM pbinorm
#' @importFrom parallel mclapply
#' @keywords internal



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
  
  
  J = length(names)  # Number of variables in the analysis.
  pairs = paste(ep[,1], ep[,2], sep = "-")  # Construct pairs from 'ep' for parameter naming.
  
  par_all[parms] = par  # Update specified parameters.
  
  parm <- param(par_all, names)
  #ax <- Matrix::nearPD(compute_ax(parm, names))$mat  # Compute ax correction terms
  parm <- param(update_ax_parameters(par_all, names, compute_ax(parm, names)), names)
  beta <- try(compute_beta(parm, names, cr), silent = T)  # Compute beta coefficients
  # Attempt Cholesky decomposition to ensure positive definiteness.
  #ae <- try(chol(ax), silent = TRUE)
  be <- try(chol(beta), silent = TRUE)
  
  if (!is.character(be)) {
    # Proceed if both 'ax' and 'beta' matrices are valid for further computations.
    
    # Map parameters to each variable pair in 'Vi'.
    parmm = lapply(1:nrow(Vi), function(v){
      as.numeric(parm[(parm$v1==Vi[v,1] & parm$v2 ==Vi[v,2])|(parm$v1==Vi[v,2] & parm$v2 ==Vi[v,1]),][,-c(1,2)])
    })
    u = uh[,1]
    h = uh[,2]
    
    # Parallel computation of log-likelihood for each pair using mclapply (if multicore is intended, else lapply).
    ll = parallel::mclapply(1:nrow(Vi), function(v) {
      # Initialize log-likelihood components for the current pair.
      l1 = l2 = l3 = l4 = 0
      par = parmm[[v]]  # Parameters for the current pair.
      # Validate parameter constraints; return a large penalty if violated.
      if (any(par[c(1:25)] < 0) | any(par[c(7:21)] > 1)) {
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
#' Log-Likelihood for Spatial Data
#'
#' Calculates the log-likelihood for spatial data based on the Matérn covariance function. This function plays a pivotal role in estimating spatial parameters for geostatistical models.
#'
#' @param par Vector containing parameters for the Matérn covariance function: range (`par[1]`) and smoothness (`par[2]`). Both parameters must be positive.
#' @param data 3D array of observed spatial data.
#' @param h Vector of spatial distances between observations.
#' @param uh Matrix specifying indices for pairing spatial observations.
#' @param v Index of the variable within `data` for which the log-likelihood is computed.
#'
#' @return Log-likelihood value for the spatial data under the Matérn covariance model.
#'
#' @importFrom VGAM pbinorm
#' @keywords internal

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
      l3 = sum(pbinorm(uh[id3, 7], uh[id3, 8],var1 = 1, var2 = 1, cov12 = cij[id3]))
    }
    
    # Return the aggregated negative log-likelihood, adjusting for errors or infinite values.
    ll = try(-(l1 + l2 + l3 + l4), silent = TRUE)
    if(is.character(ll) || is.infinite(ll)) ll = abs(rnorm(1)) * 1e+20
    return(ll)
  }
}
#' Compute Spatio-Temporal Covariances
#'
#' Calculates spatial and temporal covariances for given spatio-temporal data, facilitating the understanding of spatial and temporal variability in the context of different weather types.
#'
#' @param data 3D array representing time, location, and variable dimensions of the spatio-temporal data.
#' @param wt_id Indices of weather types for which covariances are computed.
#' @param locations Matrix of spatial locations for the data points.
#' @param ds Precomputed distance matrix or NULL to compute distances from 'locations'.
#' @param dates Vector of dates corresponding to the time dimension of the data.
#' @param lagstime Vector of time lags for covariance computation.
#' @param dist Vector of spatial distances for covariance computation.
#' @param covgm Logical flag to compute cross-covariances (default TRUE).
#'
#' @return Data frame containing computed covariances for specified spatial distances and time lags, facilitating the analysis of spatial and temporal patterns in the data.
#'
#' @keywords internal


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
      idx <- which(ds > d-1 & ds < d+1, arr.ind = TRUE)
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
#' Generate Covariance Matrices for Spatio-Temporal Model
#'
#' Creates covariance matrices for each time lag and pair of variables using Gneiting's function, based on provided model parameters and spatial locations. These matrices are essential for multivariate space-time modeling, reflecting the covariance structure across space and time.
#'
#' @param par Parameters for the Gneiting covariance function, including details for variable pairs.
#' @param coordinates Matrix or data frame containing spatial coordinates for each location.
#' @param names Vector of variable names involved in the covariance calculations.
#' @param M Maximum time lag considered in the model.
#'
#' @return A list of covariance matrices for each time lag up to M, and for each pair of variables, where each matrix represents the spatial covariance structure for a given time lag and variable pair.
#'
#' @keywords internal

cov_matrices = function(par, coordinates, names, M) {
  # Function for generating covariance matrices for a multivariate space-time model.
  # The covariance is calculated using Gneiting's function.
  #
  # Arguments:
  #   par: Parameters for the covariance function
  #   coordinates: A matrix or data frame of coordinates coordinates for each spatial location.
  #   names: Names of the variables involved in the covariance calculation.
  #   M: The maximum time lag considered in the model.
  #
  # Returns:
  #   A list of covariance matrices for each time lag (up to M) and for each pair of variables.
  #   Each matrix represents the spatial covariance structure for a given time lag and variable pair.
  
  Nt = M + 1  # Number of time points considered
  Ns = nrow(coordinates)  # Number of spatial locations
  Nv = length(names)  # Number of variables
  
  # Generate all combinations of time points and spatial locations
  d = expand.grid(t1 = 1:Nt, t2 = 1:Nt, s1 = 1:Ns, s2 = 1:Ns)
  
  # Calculate time lags (u) and spatial distances (h) between all pairs of points
  u = d$t1 - d$t2
  h = ds(d$s1, d$s2, coordinates)  # Calculate distances based on coordinates
  
  # Initialize a list to store covariance matrices
  cp = lapply(1:Nt, function(t1) {
    cp_v1 = lapply(names, function(v1) {
      cp_v2 = lapply(names, function(v2) {
        # Retrieve parameters for the current pair of variables and calculate covariance
        cov_params = par[(par$v1 == v1 & par$v2 == v2) | (par$v2 == v1 & par$v1 == v2), -c(1, 2)]
        dij = par$dij[(par$v1 == v1 & par$v2 == v2) | (par$v2 == v1 & par$v1 == v2)]
        cov = Gneiting(h, u, cov_params, dij)
        
        # Filter to the current time point and reshape the covariance values into a matrix
        up = (d$t1 == t1) & (d$t2 == 1)
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
#' Check Positive Definiteness Condition
#'
#' Verifies the positive definiteness of the covariance matrix constructed from model parameters, which is crucial for ensuring valid covariance structures in spatio-temporal modeling.
#'
#' @param parm A data frame or list containing the model parameters.
#' @param names Character vector of variable names for which the condition is checked.
#'
#' @return Logical value indicating whether the covariance matrix, constructed based on the parameters and variable names, is positive definite.
#'
#' @keywords internal


pd_condition = function(parm, names){
  eij = sapply(names, function(v1){
    sapply(names, function(v2){
      dij = parm$dij[parm$v1==v1 & parm$v2==v2 |parm$v1==v2 & parm$v2==v1]
      vii = parm$vii[parm$v1==v1 & parm$v2==v2 |parm$v1==v2 & parm$v2==v1]
      vjj = parm$vjj[parm$v1==v1 & parm$v2==v2 |parm$v1==v2 & parm$v2==v1]
      rii = parm$rii[parm$v1==v1 & parm$v2==v2 |parm$v1==v2 & parm$v2==v1]
      rjj = parm$rjj[parm$v1==v1 & parm$v2==v2 |parm$v1==v2 & parm$v2==v1]
      aii = parm$aii[parm$v1==v1 & parm$v2==v2 |parm$v1==v2 & parm$v2==v1]
      ajj = parm$ajj[parm$v1==v1 & parm$v2==v2 |parm$v1==v2 & parm$v2==v1]
      ci = parm$ci[parm$v1==v1 & parm$v2==v2 |parm$v1==v2 & parm$v2==v1]
      cj = parm$cj[parm$v1==v1 & parm$v2==v2 |parm$v1==v2 & parm$v2==v1]
      vij <- (vii + vjj) / 2
      rij <- sqrt((rii^2 + rjj^2) / 2)
      aij <- sqrt((aii^2 + ajj^2) / 2)
      
      dij <- dij *((rii^vii * rjj^vjj) / rij^(2*vij)) *
        (gamma(vij) / (gamma(vii)^(1/2) * gamma(vjj)^(1/2))) * 
        (2-ci^2)*(2-cj^2) 
      return(dij / gamma(vij))
    })
  })
  pd = try(chol(eij), silent = T)
  return(!is.character(pd))
}
#' Modify Beta Parameters in Model Parameters
#'
#' Adjusts the 'dij' parameters in the model parameter set based on computed beta coefficients, ensuring that the covariance structure reflects these adjustments.
#'
#' @param parm A data frame or list representing the current set of model parameters, including 'v1', 'v2', and 'dij' among others.
#' @param beta A matrix of beta coefficients computed to adjust the correlations or dependencies between variables.
#'
#' @return The modified set of model parameters with updated 'dij' values based on the beta coefficients.
#'
#' @keywords internal


modify_beta_parm = function(parm, beta){
  for (i in 1:nrow(parm)) {
    parm$dij[i] = beta[parm$v1[i], parm$v2[i]] 
  }
  return(parm)
}