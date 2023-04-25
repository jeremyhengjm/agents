#' @rdname dsumbin
#' @title Probability mass function of the sum of two independent binomial random variables
#' @description Evaluates probability mass function of the sum of two independent binomial random variables
#' @param eval a vector of evaluation points
#' @param n1 number of trials of first binomial random variable
#' @param p1 probability of success of first binomial random variable
#' @param n2 number of trials of second binomial random variable
#' @param p2 probability of success of second binomial random variable
#' @param log_scale logical specifying if probabilities are evaluated on the logarithmic scale
#' @return a vector or matrix of probabilities 
#' @export

dsumbin <- function(eval, n1, p1, n2, p2, log_scale = TRUE){
  
  logpmf <- dsumbin_vectorized_cpp(eval, n1, p1, n2, p2)
  
  if (log_scale == TRUE){
    output <- logpmf
  } else{
    output <- exp(logpmf)
  }
  
  return(output)
}