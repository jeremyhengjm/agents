#' @rdname sis_simulate
#' @title Simulate SIS model 
#' @description Runs a simulation of SIS agent-based model 
#' @param model_config a list containing:
#' \itemize{
#' \item 'N': size of the population
#' \item 'alpha0': initial infection probability
#' \item 'lambda': vector of infection rates of length N 
#' \item 'gamma': vector of recovery rates of length N 
#' \item 'network_type': character "full" or "net" specifying if network is fully connected or has a given network structure
#' \item 'adjacency': matrix indicating interaction of agents of dimension N x N (if model_config$network_type == "net)
#' \item 'degree' vector of node degree of length N (if model_config$network_type == "net)
#' \item 'rmeasurement': function to sample from observation model
#' }
#' @param observations a list containing:
#' \itemize{
#' \item 'T': number of observations is T+1
#' }
#' @return a list containing: 
#' \itemize{
#' \item 'X': matrix of particle states of dimension T+1 x N 
#' \item 'I': vector of infected agents of length T+1
#' \item 'y': vector of observations of length T+1
#' }
#' @export

sis_simulate <- function(model_config, observations){
  # get model and observation settings
  N <- model_config$N
  T <- observations$T 
  
  # get model parameters
  alpha0 <- model_config$alpha0 
  
  # pre-allocate
  X <- matrix(0, nrow = T+1, ncol = N)
  I <- rep(0, T+1)
  y <- rep(0, T+1)
  
  # initialize agent states
  X[1, ] <- (runif(N) < alpha0)
  I[1] <- sum(X[1, ])
  
  # simulate observation
  y[1] <- model_config$rmeasurement(I[1])
  
  for (t in 1:T){
    # evolve agent states according to model transition
    alpha <- as.numeric(sis_alpha_prob(model_config, matrix(X[t, ], nrow = N, ncol = 1), I[t]))
    X[t+1, ] <- (runif(N) < alpha)
    I[t+1] <- sum(X[t+1, ])
    
    # simulate observation
    y[t+1] <- model_config$rmeasurement(I[t+1])
  }
  
  return(list(X = X, I = I, y = y))
  
}

