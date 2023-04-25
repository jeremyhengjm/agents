#' @rdname sis_bpf
#' @title Bootstrap particle filter for SIS model 
#' @description Runs a bootstrap particle filter for SIS agent-based model 
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
#' }
#' @param observations a list containing:
#' \itemize{
#' \item 'T': number of observations is T+1
#' \item 'y': vector of observations of length T+1
#' }
#' @param filter_config a list containing:
#' \itemize{
#' \item 'P': number of particles
#' \item 'store': logical indicating if particle states and their ancestries are stored
#' \item 'resampling': character specifying type of resampling scheme
#' }
#' @return a list containing: 
#' \itemize{
#' \item 'log_likelihood': vector of log-likelihood estimates of length T+1  
#' \item 'ess: vector of effective sample sizes of length T+1
#' \item 'runtime': algorithmic run time 
#' \item 'particles': array of particle states of dimension T+1 x N x P (if filter_config$store == TRUE)
#' \item 'ancestries': matrix of ancestor indexes of dimension T x P (if filter_config$store == TRUE)
#' }
#' @export

sis_bpf <- function(model_config, observations, filter_config){
  # start timer
  start_time <- as.numeric(Sys.time())
  
  # get model and algorithmic settings
  N <- model_config$N
  dmeasurement <- model_config$dmeasurement
  T <- observations$T
  P <- filter_config$P
  resampling <- filter_config$resampling
  
  # get model parameters
  alpha0 <- model_config$alpha0 
  
  # get observations
  y <- observations$y
  
  # pre-allocate
  if (filter_config$store){
    particles <- array(0, dim = c(T+1, N, P))
    ancestries <- matrix(0, nrow = T, ncol = P)
  }
  log_likelihood <- rep(0, T+1)
  ess <- rep(0, T+1)
  
  # initialization
  X <- (matrix(runif(P * N), nrow = N, ncol = P) < matrix(alpha0, nrow = N, ncol = P))
  if (filter_config$store){
    particles[1, , ] <- X
  }
  I <- colSums(X)
  
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
  
  for (t in 1:T){
    # evolve each particle according to model transition
    alpha <- sis_alpha_prob(model_config, X, I)
    X <- (matrix(runif(P * N), nrow = N, ncol = P) < alpha)
    I <- colSums(X)
    
    # save particle states and their ancestors
    if (filter_config$store){
      particles[t+1, , ] <- X
      ancestries[t, ] <- A
    }
    
    # compute weights and ESS
    logweights <- dmeasurement(y[t+1], I)
    maxlogweights <- max(logweights)
    weights <- exp(logweights - maxlogweights)
    normweights <- weights / sum(weights)
    ess[t+1] <- 1 / sum(normweights^2)
    
    # compute log-likelihood estimate
    log_incremental_likelihood <- log(mean(weights)) + maxlogweights  
    log_current_likelihood <- log_current_likelihood + log_incremental_likelihood
    log_likelihood[t+1] <- log_current_likelihood
    
    # resampling
    if (t < T){
      A <- resampling(normweights, P) 
      X <- X[, A]
      I <- I[A]
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
