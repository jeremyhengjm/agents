#' @rdname multinomial_resampling
#' @title Multinomial resampling scheme
#' @description Implements the multinomial resampling scheme
#' @param weights a vector of normalized weights
#' @param ndraws integer specifying number of ancestor indexes required
#' @return a vector of ancestor indexes
#' @export

multinomial_resampling <- function(weights, ndraws){
  ancestors <- multinomial_resampling(weights, ndraws)
  return(ancestors)
}

#' @rdname systematic_resampling
#' @title Systematic resampling scheme
#' @description Implements the systematic resampling scheme
#' @param weights a vector of normalized weights
#' @param ndraws integer specifying number of ancestor indexes required
#' @return a vector of ancestor indexes
#' @export

systematic_resampling <- function(weights, ndraws){
  ancestors <- systematic_resampling_cpp(weights, ndraws)
  return(ancestors)
}
