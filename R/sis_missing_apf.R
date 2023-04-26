#' @rdname sis_missing_apf
#' @title Auxiliary particle filter for SIS model with missing data
#' @description Runs an auxiliary particle filter for SIS agent-based model with missing data
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
#' \item 'log_bif': approximate backward information filter stored as matrix of dimension N+1 x T+1
#' \item 'dpoibin': function to compute or approximate Poisson binomial probability mass function 
#' \item 'rcondber': function to generate samples from conditional Bernoulli distribution 
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

sis_missing_apf <- function(model_config, observations, filter_config){
  # start timer
  start_time <- as.numeric(Sys.time())
  
  # get model and algorithmic settings
  N <- model_config$N
  nsteps <- model_config$nsteps
  dmeasurement <- model_config$dmeasurement
  infected_support <- model_config$infected_support
  T <- observations$T
  P <- filter_config$P
  resampling <- filter_config$resampling
  dpoibin <- filter_config$dpoibin
  rcondber <- filter_config$rcondber
  
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
  
  # check if there is an observation at initial time
  if (obs_times[1]){
    # counter to index observation
    index_obs <- index_obs + 1 
    
    # compute controlled distribution of number of infected 
    I_support <- infected_support(y[1]) # M = length(I_support)
    log_poibin <- matrix(dpoibin(I_support, alpha0), ncol = 1) # M x 1
    log_bif <- filter_config$log_bif[I_support+1, 1] # vector of size M 
    infected_distribution <- sis_compute_infected_prob(log_poibin, log_bif) 
    normprob <- as.vector(infected_distribution$normprob)
    log_exp_bif <- infected_distribution$log_normconstant
    
    # sample controlled number of infected and agent states 
    I <- rcat(I_support, normprob, P)
    X <- rcondber(I, alpha0) # N x P
    if (filter_config$store){
      particles[1, , ] <- X
    }
    
    # compute first controlled transition of number of infected 
    alpha <- sis_alpha_prob(model_config, X, I) # N x P 
    if (obs_times[2]){
      I_support <- infected_support(y[2]) # M = length(I_support)
    } else{
      I_support <- 0:N  
    }
    log_poibin <- dpoibin(I_support, alpha) # M x P
    log_bif <- filter_config$log_bif[I_support+1, 2]
    infected_transition <- sis_compute_infected_prob(log_poibin, log_bif) 
    normprob <- infected_transition$normprob # M x P
    log_condexp_bif <- infected_transition$log_normconstant # vector of size P
    I_zeros <- (I == 0)
    if (any(I_zeros)){
      index_zeros <- which(I_zeros)
      if (I_support[1] == 0){
        log_condexp_bif[index_zeros] <- log_bif[1]
      } else{
        log_condexp_bif[index_zeros] <- -Inf
      }
    }
    
    # compute controlled weights and ESS 
    log_dmeasurement <- dmeasurement(y[1], I) # vector of size P
    log_bif_values <- filter_config$log_bif[I+1, 1] # vector of size P
    logweights <- log_exp_bif + log_dmeasurement + log_condexp_bif - log_bif_values
    maxlogweights <- max(logweights)
    weights <- exp(logweights - maxlogweights)
    normweights <- weights / sum(weights)
    ess[1] <- 1 / sum(normweights^2)
    
    # compute log-likelihood estimate
    log_incremental_likelihood <- log(mean(weights)) + maxlogweights  
    log_current_likelihood <- log_incremental_likelihood
    log_likelihood[1] <- log_current_likelihood
  }
  
  for (t in 1:nsteps){
    # resampling
    A <- resampling(normweights, P)
    normprob <- normprob[, A]
    alpha <- alpha[, A]
    
    # sample controlled number of infected and agent states 
    I <- rcat(I_support, normprob)
    X <- rcondber(I, alpha) # N x P
    
    # save particle states and their ancestors
    if (filter_config$store){
      particles[t+1, , ] <- X
      ancestries[t, ] <- A
    }
    
    # check if there is an observation at this time
    if (obs_times[t+1]){ 
      # counter to index observation
      index_obs <- index_obs + 1 
      log_dmeasurement <- dmeasurement(y[index_obs], I) # vector of size P
    } else{
      log_dmeasurement <- rep(0, P)
    }
    
    # compute next controlled transition of number of infected 
    if (t < nsteps){
      alpha <- sis_alpha_prob(model_config, X, I) # N x P
      if (obs_times[t+2]){
        I_support <- infected_support(y[index_obs+1]) # M = length(I_support)
      } else{
        I_support <- 0:N
      }
      log_poibin <- dpoibin(I_support, alpha) # M x P
      log_bif <- filter_config$log_bif[I_support+1, t+2]
      infected_transition <- sis_compute_infected_prob(log_poibin, log_bif)
      normprob <- infected_transition$normprob # M x P
      log_condexp_bif <- infected_transition$log_normconstant # vector of size P
      I_zeros <- (I == 0)
      if (any(I_zeros)){
        index_zeros <- which(I_zeros)
        if (I_support[1] == 0){
          log_condexp_bif[index_zeros] <- log_bif[1]
        } else{
          log_condexp_bif[index_zeros] <- -Inf
        }
      }
    } else{
      log_condexp_bif <- rep(0, P)
    }
    
    # compute controlled weights and ESS 
    log_bif_values <- filter_config$log_bif[I+1, t+1] # vector of size P
    logweights <- log_dmeasurement + log_condexp_bif - log_bif_values
    maxlogweights <- max(logweights)
    weights <- exp(logweights - maxlogweights)
    normweights <- weights / sum(weights)
    ess[t+1] <- 1 / sum(normweights^2)
    
    # compute log-likelihood estimate
    log_incremental_likelihood <- log(mean(weights)) + maxlogweights  
    log_current_likelihood <- log_current_likelihood + log_incremental_likelihood
    log_likelihood[t+1] <- log_current_likelihood
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
