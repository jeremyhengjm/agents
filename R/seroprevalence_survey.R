#' @rdname seroprevalence_survey_model
#' @title Seroprevalence survey model
#' @description Returns objects defining seroprevalence survey model
#' @param sample_size sample size of survey
#' @param population_size population size 
#' @return a list containing: 
#' \itemize{
#' \item 'rmeasurement': function to sample from observation model
#' \item 'dmeasurement: function evaluating observation density
#' \item 'infected_support': function returning support of number of infected given observation  
#' }
#' @export

seroprevalence_survey_model <- function(sample_size, population_size){
  
  # population size
  N <- population_size
  
  # observation model for seroprevalence surveys
  rmeasurement <- function(I){
    # function to sample from observation model
    return(rbinom(n = 1, size = sample_size, prob = I / N))
  }
  
  # function evaluating observation density
  dmeasurement <- function(y, I){
    return(dbinom(x = as.numeric(y), size = sample_size, prob = I / N, log = TRUE))
  }
  
  # function returning support of number of infected given observation
  infected_support <- function(y){
    if (y == 0){
      I_support <- 0:(N-1)
    } else if (y == sample_size){
      I_support <- 1:N
    } else{
      I_support <- 1:(N-1)
    }
    return(I_support)
  }
  
  return(list(rmeasurement = rmeasurement, 
              dmeasurement = dmeasurement, 
              infected_support = infected_support))
}