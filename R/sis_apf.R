#' @rdname sis_apf
#' @title Auxiliary particle filter for SIS model 
#' @description Runs an auxiliary particle filter for SIS agent-based model 
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
#' \item 'infected_support': function returning support of number of infected given observation 
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
#' \item 'dpoibin': function to compute or approximate Poisson binomial probability mass function 
#' \item 'rcondber': function to generate samples from conditional Bernoulli distribution 
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


sis_apf <- function(model_config, observations, filter_config){
  # start timer
  start_time <- as.numeric(Sys.time())
  
  # get model and algorithmic settings
  N <- model_config$N
  dmeasurement <- model_config$dmeasurement
  infected_support <- model_config$infected_support
  T <- observations$T
  P <- filter_config$P
  resampling <- filter_config$resampling
  dpoibin <- filter_config$dpoibin
  rcondber <- filter_config$rcondber
  mcmc_iterations <- filter_config$mcmc_iterations
  
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
  
  # compute distribution of number of infected given first observation
  I_support <- infected_support(y[1]) # M = length(I_support)
  log_poibin <- matrix(dpoibin(I_support, alpha0), ncol = 1) # M x 1
  log_dmeasurement <- dmeasurement(y[1], I_support) # vector of size M 
  infected_distribution <- sis_compute_infected_prob(log_poibin, log_dmeasurement) 
  normprob <- as.vector(infected_distribution$normprob)
  
  # sample number of infected and agent states given first observation 
  I <- rcat(I_support, normprob, P)
  X <- rcondber(I, alpha0) # N x P
  if (filter_config$store){
    particles[1, , ] <- X
  }
  
  # compute ESS and log-likelihood estimate (no resampling needed)
  ess[1] <- P
  log_incremental_likelihood <- infected_distribution$log_normconstant
  log_current_likelihood <- log_incremental_likelihood
  log_likelihood[1] <- log_current_likelihood
  
  for (t in 1:T){
    # compute conditional distribution of number of infected given current observation
    alpha <- sis_alpha_prob(model_config, X, I) # N x P 
    I_support <- infected_support(y[t+1]) # M = length(I_support)
    log_poibin <- dpoibin(I_support, alpha) # M x P
    log_dmeasurement <- dmeasurement(y[t+1], I_support) # vector of size M 
    infected_transition <- sis_compute_infected_prob(log_poibin, log_dmeasurement) 
    log_normconstant <- infected_transition$log_normconstant
    I_zeros <- (I == 0)
    if (any(I_zeros)){
      index_zeros <- which(I_zeros)
      if (I_support[1] == 0){
        log_normconstant[index_zeros] <- log_dmeasurement[1]
      } else{
        log_normconstant[index_zeros] <- -Inf
      }
    }
    normprob <- infected_transition$normprob # M x P
    
    # compute weights and ESS
    logweights <- log_normconstant # vector of size P
    maxlogweights <- max(logweights)
    weights <- exp(logweights - maxlogweights)
    normweights <- weights / sum(weights)
    ess[t+1] <- 1 / sum(normweights^2)
    
    # compute log-likelihood estimate
    log_incremental_likelihood <- log(mean(weights)) + maxlogweights  
    log_current_likelihood <- log_current_likelihood + log_incremental_likelihood
    log_likelihood[t+1] <- log_current_likelihood
    
    # resampling
    A <- resampling(normweights, P)
    normprob <- normprob[, A]
    alpha <- alpha[, A]
    
    # sample number of infected and agent states given current observation
    I <- rcat(I_support, normprob)
    X <- rcondber(I, alpha) # N x P
    
    # save particle states and their ancestors
    if (filter_config$store){
      particles[t+1, , ] <- X
      ancestries[t, ] <- A
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
