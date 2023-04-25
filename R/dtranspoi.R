#' @rdname dtranspoi_poibin_approx
#' @title Probability mass function of a translated Poisson approximation of a Poisson binomial distribution
#' @description Evaluates probability mass function of a translated Poisson approximation of a Poisson binomial distribution
#' @param eval a vector of evaluation points
#' @param alpha a vector or matrix of probabilities 
#' @param log_scale logical specifying if probabilities are evaluated on the logarithmic scale
#' @return a vector or matrix of probabilities 
#' @export

dtranspoi_poibin_approx <- function(eval, alpha, log_scale = TRUE){
  
  if (is.vector(alpha)){
    logpmf <- dtranspoi_poibin_approx_cpp(eval, alpha)
  } 
  
  if (is.matrix(alpha)){
    logpmf <- dtranspoi_poibin_approx_vectorized_cpp(eval, alpha)
  }
  
  if (log_scale == TRUE){
    output <- logpmf
  } else{
    output <- exp(logpmf)
  }
  
  return(output)
}

#' @rdname dtranspoi
#' @title Probability mass function of a translated Poisson distribution
#' @description Evaluates probability mass function of a translated Poisson distribution
#' @param eval a vector of evaluation points
#' @param mu a vector of means
#' @param sigmasq a vector of variances
#' @param end_support upper end of support
#' @param log_scale logical specifying if probabilities are evaluated on the logarithmic scale
#' @return a vector or matrix of probabilities 
#' @export

dtranspoi <- function(eval, mu, sigmasq, end_support, log_scale = TRUE){
  
  if (length(mu) == 1){
    logpmf <- dtranspoi_cpp(eval, mu, sigmasq, end_support)
  } else{
    logpmf <- dtranspoi_vectorized_cpp(eval, mu, sigmasq, end_support)
  }
  
  if (log_scale == TRUE){
    output <- logpmf
  } else{
    output <- exp(logpmf)
  }
  
  return(output)
}