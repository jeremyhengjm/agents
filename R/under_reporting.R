#' @rdname under_reporting_model
#' @title Under reporting observation model
#' @description Returns objects defining under reporting observation model
#' @param rho under reporting parameter value
#' @return a list containing: 
#' \itemize{
#' \item 'rmeasurement': function to sample from observation model
#' \item 'dmeasurement: function evaluating observation density
#' \item 'infected_support': function returning support of number of infected given observation  
#' }
#' @export

under_reporting_model <- function(rho){
  
  # function to sample from observation model
  rmeasurement <- function(I){  
    return(rbinom(n = 1, size = I, prob = rho))
  }
  
  # function evaluating observation density
  dmeasurement <- function(y, I){
    return(dbinom(x = as.numeric(y), size = I, prob = rho, log = TRUE))
  }
  
  # function returning support of number of infected given observation  
  infected_support <- function(y){
    return(y:N)
  }
  
  return(list(rmeasurement = rmeasurement, 
              dmeasurement = dmeasurement, 
              infected_support = infected_support))
}