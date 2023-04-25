#' @rdname sis_missing_bif
#' @title Backward information filter approximation for SIS model with missing data
#' @description Approximate backward information filter for SIS agent-based model with missing data
#' @param model_config a list containing:
#' \itemize{
#' \item 'N': size of the population
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
#' @param conditional_expectation a character specifying how conditional expectation is computed or approximated:
#' \itemize{
#' \item 'sumbin': exact computation of conditional expectations
#' \item 'transpoi': translated Poisson approximation of conditional expectations
#' }
#' @return approximate backward information filter stored as matrix of dimension N+1 x T+1
#' @export

sis_missing_bif <- function(model_config, observations, conditional_expectation = 'sumbin'){
  # start timer
  start_time <- as.numeric(Sys.time())
  
  # get model and algorithmic settings
  N <- model_config$N
  nsteps <- model_config$nsteps
  dmeasurement <- model_config$dmeasurement
  infected_support <- model_config$infected_support
  T <- observations$T
  
  # get model parameters
  lambda_mean <- mean(model_config$lambda)
  gamma_mean <- mean(model_config$gamma)
  
  # get observations
  y <- observations$y
  obs_times <- observations$times
  
  # pre-allocate
  log_bif <- matrix(-Inf, nrow = N+1, ncol = nsteps+1)
  
  # backward information filter at the terminal time
  index_obs <- T+1
  I_current_support <- infected_support(y[index_obs])
  log_bif[I_current_support+1, nsteps+1] <- dmeasurement(y[index_obs], I_current_support) 
  
  # pre-compute transition probability of number of infected under approximate model
  I_current <- 0:N
  I_next <- 0:N
  if (conditional_expectation == 'sumbin'){
    log_transition <- dsumbin(I_next, N - I_current, lambda_mean * I_current / N, I_current, 1 - gamma_mean) # N+1 x N+1 
  }
  if (conditional_expectation == 'transpoi'){
    precompute_factor <- lambda_mean * I_current / N
    transition_mean <- precompute_factor * (N - I_current) + (1 - gamma_mean) * I_current
    transition_var <- transition_mean - precompute_factor^2 * (N - I_current) - (1 - gamma_mean)^2 * I_current
    log_transition <- dtranspoi(I_next, transition_mean, transition_var, N)
  }
  log_transition[2:(N+1),1] <- -Inf
  
  # approximate backward recursion
  for (t in nsteps:1){
    I_next_support <- I_current_support
    if (obs_times[t]){
      index_obs <- index_obs - 1 
      I_current_support <- infected_support(y[index_obs])
      log_dmeasurement <- dmeasurement(y[index_obs], I_current_support) 
    } else{
      I_current_support <- 0:N
      log_dmeasurement <- rep(0, N+1)
    }
    log_transition_support <- log_transition[I_next_support+1, I_current_support+1]
    log_policy <- log_bif[I_next_support+1, t+1]
    log_condexp_bif <- sis_approximate_infected_prob(log_transition_support, log_policy)
    if (I_current_support[1] == 0){
      if (I_next_support[1] == 0){
        log_condexp_bif[1] <- log_policy[1]
      } else{
        log_condexp_bif[1] <- -Inf
      }
    } 
    log_bif[I_current_support+1, t] <- log_dmeasurement + log_condexp_bif
  }
  
  # end timer
  end_time <- as.numeric(Sys.time())
  run_time <- end_time - start_time
  
  # output results
  results <- list(log_bif = log_bif, run_time = run_time)
  return(results)
}
