#' @rdname dpoibin_exact
#' @title Poisson binomial probability mass function
#' @description Evaluates Poisson binomial probability mass function exactly
#' @param eval a vector of evaluation points
#' @param alpha a vector or matrix of probabilities 
#' @param log_scale logical specifying if probabilities are evaluated on the logarithmic scale
#' @return a vector or matrix of probabilities 
#' @export

dpoibin_exact <- function(eval, alpha, log_scale = TRUE){
  if (is.vector(alpha)){
    logpmf <- dpoibin_cpp(eval, alpha)
  } 
  
  if (is.matrix(alpha)){
    logpmf <- dpoibin_vectorized_cpp(eval, alpha)
  }
  
  if (log_scale == TRUE){
    output <- logpmf
  } else{
    output <- exp(logpmf)
  }
  
  return(output)
}