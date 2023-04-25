#' @rdname sis_missing_bpf
#' @title Bootstrap particle filter for SIS model with missing data
#' @description Runs a bootstrap particle filter for SIS agent-based model with missing data
#' @param model_config a list containing:
#' \itemize{
#' \item 'N': size of the population
#' \item 'alpha0': vector of initial infection probability of length N
#' \item 'lambda': vector of infection rates of length N 
#' \item 'gamma': vector of recovery rates of length N 
#' \item 'network_type': character "full" or "net" specifying if network is fully connected or has a given network structure
#' \item 'adjacency': matrix indicating interaction of agents of dimension N x N (if model_config$network_type == "net)
#' \item 'degree' vector of node degree of length N (if model_config$network_type == "net)
#' \item 'dmeasurement': function evaluating observation density
#' \item 'nsteps': number of time steps
#' }
#' @param observations a list containing:
#' \itemize{
#' \item 'T': number of observations is T+1
#' \item 'y': vector of observations of length T+1
#' \item 'times': vector of length nsteps+1 indicating the observation times
#' }
#' @param filter_config a list containing:
#' \itemize{
#' \item 'P': number of particles
#' \item 'store': logical indicating if particle states and their ancestries are stored
#' \item 'resampling': character specifying type of resampling scheme
#' }
#' @return a list containing: 
#' \itemize{
#' \item 'log_likelihood': vector of log-likelihood estimates of length nsteps+1  
#' \item 'ess: vector of effective sample sizes of length nsteps+1
#' \item 'runtime': algorithmic run time 
#' \item 'particles': array of particle states of dimension nsteps+1 x N x P (if filter_config$store == TRUE)
#' \item 'ancestries': matrix of ancestor indexes of dimension nsteps x P (if filter_config$store == TRUE)
#' }
#' @export

sis_missing_bpf <- function(model_config, observations, filter_config){
  # start timer
  start_time <- as.numeric(Sys.time())
  
  # get model and algorithmic settings
  N <- model_config$N
  dmeasurement <- model_config$dmeasurement
  nsteps <- model_config$nsteps
  T <- observations$T
  P <- filter_config$P
  resampling <- filter_config$resampling
  
  # get model parameters
  alpha0 <- model_config$alpha0 
  
  # get observations
  y <- observations$y
  obs_times <- observations$times
  
  # pre-allocate
  if (filter_config$store){
    particles <- array(0, dim = c(nsteps+1, N, P))
    ancestries <- matrix(0, nrow = nsteps, ncol = P)
  }
  log_likelihood <- rep(0, nsteps+1)
  ess <- rep(P, nsteps+1)
  index_obs <- 0
  log_current_likelihood <- 0
  
  # initialization
  X <- (matrix(runif(P * N), nrow = N, ncol = P) < matrix(alpha0, nrow = N, ncol = P))
  if (filter_config$store){
    particles[1, , ] <- X
  }
  I <- colSums(X)
  
  # check if there is an observation at initial time
  if (obs_times[1]){
    # counter to index observation
    index_obs <- index_obs + 1 
    
    # compute weights and ESS 
    logweights <- dmeasurement(y[1], I)
    maxlogweights <- max(logweights)
    weights <- exp(logweights - maxlogweights)
    normweights <- weights / sum(weights)
    ess[1] <- 1 / sum(normweights^2)
    
    # compute log-likelihood estimate
    log_incremental_likelihood <- log(mean(weights)) + maxlogweights  
    log_current_likelihood <- log_incremental_likelihood
    log_likelihood[1] <- log_current_likelihood
    
    # resampling
    A <- resampling(normweights, P) 
    X <- X[, A]
    I <- I[A]
  }
  
  for (t in 1:nsteps){
    # evolve each particle according to model transition
    alpha <- sis_alpha_prob(model_config, X, I)
    X <- (matrix(runif(P * N), nrow = N, ncol = P) < alpha)
    I <- colSums(X)
    
    # save particle states and their ancestors
    if (filter_config$store){
      particles[t+1, , ] <- X
      ancestries[t, ] <- A
    }
    
    # check if there is an observation at this time
    if (obs_times[t+1]){
      # counter to index observation
      index_obs <- index_obs + 1 
      
      # compute weights and ESS
      logweights <- dmeasurement(y[index_obs], I)
      maxlogweights <- max(logweights)
      weights <- exp(logweights - maxlogweights)
      normweights <- weights / sum(weights)
      ess[t+1] <- 1 / sum(normweights^2)
      
      # compute log-likelihood estimate
      log_incremental_likelihood <- log(mean(weights)) + maxlogweights  
      log_current_likelihood <- log_current_likelihood + log_incremental_likelihood
      log_likelihood[t+1] <- log_current_likelihood
      
      # resampling
      if (t < nsteps){
        A <- resampling(normweights, P) 
        X <- X[, A]
        I <- I[A]
      }
    } else{
      A <- 1:P
      log_likelihood[t+1] <- log_current_likelihood
    }
  }
  
  # end timer
  end_time <- as.numeric(Sys.time())
  run_time <- end_time - start_time
  
  # output results
  results <- list(log_likelihood = log_likelihood, 
                  ess = ess,
                  run_time = run_time)
  if (filter_config$store){
    results$particles <- particles
    results$ancestries <- ancestries
  } 
  return(results)
}
