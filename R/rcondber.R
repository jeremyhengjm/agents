#' @rdname rcondber_exact
#' @title Random sample from conditional Bernoulli distribution 
#' @description Generate a sample from conditional Bernoulli distribution
#' @param conditioned_sum a vector of conditioned sums
#' @param alpha a vector or matrix of probabilities 
#' @return a logical matrix 
#' @export

rcondber_exact <- function(conditioned_sum, alpha){
  if (is.vector(alpha)){
    N <- length(alpha)
    P <- length(conditioned_sum)
    output <- rcondber_vectorized_cpp(conditioned_sum, matrix(alpha, nrow = N, ncol = P))
  }
  
  if (is.matrix(alpha)){
    output <- rcondber_vectorized_cpp(conditioned_sum, alpha)
  }
  
  return(output)
}

#' @rdname rcondber_mcmc
#' @title Random sample from conditional Bernoulli distribution using Markov chain Monte Carlo
#' @description Generate a sample from conditional Bernoulli distribution using Markov chain Monte Carlo
#' @param conditioned_sum a vector of conditioned sums
#' @param alpha a vector or matrix of probabilities 
#' @param mcmc_iterations number of Markov chain Monte Carlo iterations
#' @return a logical matrix 
#' @export

rcondber_mcmc <- function(conditioned_sum, alpha, mcmc_iterations){
  if (is.vector(alpha)){
    N <- length(alpha)
    P <- length(conditioned_sum)
    output <- rcondber_mcmc_vectorized_cpp(conditioned_sum, matrix(alpha, nrow = N, ncol = P), mcmc_iterations)
  }
  
  if (is.matrix(alpha)){
    output <- rcondber_mcmc_vectorized_cpp(conditioned_sum, alpha, mcmc_iterations)
  }
  
  return(output)
}
