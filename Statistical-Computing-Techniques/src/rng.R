
GenerateRV <- function(n, method = "inverse", ...){
  arg <- list(...)
  
  # method of inverse function
  if(method == "inverse"){
    if(is.null(arg$InvCDF)) stop("The inverse function method requires the 'InvCDF' function parameter! ")
    
    u <- runif(n)
    x <- arg$InvCDF(u)
    
    return(x)
  }
  
}
