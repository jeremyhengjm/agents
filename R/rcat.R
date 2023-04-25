#' @rdname rcat
#' @title Random sample from categorical distribution 
#' @description Generate a sample from categorical distribution
#' @param support a vector specifying the support
#' @param prob a vector or matrix of probabilities 
#' @param num_samples number of samples required (default is one)
#' @return a vector of samples
#' @export

rcat <- function(support, prob, num_samples = 1){
  if (is.vector(prob)){
    output <- rcat_cpp(support, prob, num_samples)
  }
  
  if (is.matrix(prob)){
    output <- rcat_vectorized_cpp(support, prob)
  }
  
  return(output)
}
