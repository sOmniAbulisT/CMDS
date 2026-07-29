#' Generate Random Variables
#'
#' A universal wrapper function for generating random variables using five classical statistical methods: 
#' Inverse Transform, Acceptance-Rejection, Transformation, Sums, and Mixtures. 
#' 
#' @param n Integer. The number of random variables to generate. 
#' @param method Character. The sampling algorithm to use. Default is "inverse". 
#'    Supported methods: "inverse", "rejection", "transformation", "sum", "mixture". 
#' @param ... Dynamic arguments (Dots). Required parameters depend on the chosen method:
#'   \itemize{
#'     \item{\strong{inverse}:} \code{InvCDF} (Inverse function of the target CDF).
#'     \item{\strong{rejection}:} \code{f} (Target PDF), \code{g} (Proposal PDF), \code{rg} (Proposal random generator), \code{c} (Bounding constant).
#'     \item{\strong{transformation}:} \code{base_rg} (List of base random generators), \code{transform} (Mathematical transformation function).
#'     \item{\strong{sum}:} \code{base_rg} (Base random generator function), \code{k} (Number of variables to sum per sample).
#'     \item{\strong{mixture}:} \code{mix_rg} (Stage-1 parameter generator), \code{base_rg} (Stage-2 conditional generator).
#'   }
#'   
#' @return A numeric vector of length n containing the generated random variables. 
#' @export
#' 
GenerateRV <- function(n, method = "inverse", ...){
  arg <- list(...)
  
  # method of inverse function
  if(method == "inverse"){
    if(is.null(arg$InvCDF)) stop("The inverse function method requires the 'InvCDF' function parameter! ")
    
    u <- runif(n) 
    x <- arg$InvCDF(u)
    
    return(x)
    
  # method of acceptance-rejection
  } else if (method == "rejection"){
    if(is.null(arg$f) || is.null(arg$g) || is.null(arg$rg) || is.null(arg$c)){
      stop("The rejection method requires 'f', 'g', 'rg', and 'c' parameters! ")
    }
    
    k <- 0 # counter for accepted
    Iter <- 0 # iterations
    y <- numeric(n) 
    
    while(k < n){
      u <- runif(1)
      Iter <- Iter + 1
      x <- arg$rg(1) # random from g
      
      accept <- arg$f(x) / (arg$g(x) * arg$c)
      
      if(u < accept){
        # accept x
        k <- k + 1
        y[k] <- x
      }
    }
    
    return(y)
  # method of transformation 
  } else if (method == "transformation"){
    if(is.null(arg$transform) || is.null(arg$base)){
      stop("The transformation method requires 'transform' (a formula function) and 'base' (a list of generator functions)!")
    }
    
    samples <- lapply(arg$base, function(f) f(n))
    y <- do.call(arg$transform, samples)
    
    return(y)
  # method of sums
  }else if (method == "sum"){
    if(is.null(arg$gf) || is.null(arg$k)){
      stop("The sum method requires 'gf' (generator function) and 'k' (number of variables to sum)!")
    }
    
    mat <- matrix(arg$gf(n * arg$k), ncol = arg$k)
    y <- rowSums(mat)
    
    return(y) 
  # method of mixtures
  } else if (method == "mixture"){
    if(is.null(arg$rg1) || is.null(arg$rg2)){
      stop("The mixture method requires 'rg1' (stage-1 parameter generator) and 'rg2' (stage-2 conditional generator)!")
    }
    
    params <- arg$rg1(n)
    y <- arg$rg2(n, params)
    
    return(y)
  }else{
    stop("Invalid method! Please choose from: 'inverse', 'rejection', 'transformation', 'sum', or 'mixture'.")
  }
  
}

