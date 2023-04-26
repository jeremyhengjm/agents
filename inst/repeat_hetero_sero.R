# Experiments to illustrate impact of heterogeneity
rm(list = ls())
library(agents)

# load simulated data
load('inst/data_hetero_sero.RData')
T <- observations$T

# SIS model
N <- 100
alpha0 <- rep(1 / N, N)
link <- function(u) 1 / (1 + exp(-u))
network_type <- 'full'
beta_slopes <- c(0.125, 0.25, 0.5, 1.0)

# observation model for seroprevalence surveys
sample_size <- ceiling(0.2 * N)
obs_model <- seroprevalence_survey_model(sample_size, N)

# number of repetitions
num_repeats <- 100

# number of particles
num_particles <- c(2^7, 2^8, 2^9, 2^10, 2^11)

# store results
results <- data.frame()

for (slope in beta_slopes){
  # vary heterogeneity with slopes
  beta_lambda <- c(1, slope) 
  beta_gamma <- c(-1, slope)
  lambda <- link(beta_lambda[1] + beta_lambda[2] * covariates)
  gamma <- link(beta_gamma[1] + beta_gamma[2] * covariates)
  
  # define model configuration
  model_config <- list(N = N, 
                       alpha0 = alpha0, 
                       lambda = lambda,
                       gamma = gamma,
                       network_type = network_type, 
                       rmeasurement = obs_model$rmeasurement,
                       dmeasurement = obs_model$dmeasurement, 
                       infected_support = obs_model$infected_support)
  
  # approximate backward information filter
  bif <- sis_bif(model_config, observations) # cost O(TN^3)
  cat('BIF run-time:', bif$run_time, 'secs', '\n')
  
  # vary number of particles
  for (P in num_particles){
    # configuration of particle filters
    bpf_config <- list(P = P, store = FALSE, resampling = systematic_resampling)
    apf_config <- list(P = P, store = FALSE, resampling = systematic_resampling, 
                       dpoibin = dpoibin_exact, rcondber = rcondber_exact)
    cpf_config <- list(P = P, store = FALSE, resampling = systematic_resampling, 
                       log_bif = bif$log_bif, dpoibin = dpoibin_exact, rcondber = rcondber_exact)
    
    # repeat each particle filter
    for (r in 1:num_repeats){
      # run auxiliary particle filter
      apf <- sis_apf(model_config, observations, apf_config)
      cat('APF:', apf$log_likelihood[T+1], '\n')
      results <- rbind(results, data.frame(slope = slope, 
                                           particles = P, 
                                           repetition = r, 
                                           filter = 'APF',
                                           loglikelihood = apf$log_likelihood[T+1], 
                                           runtime = apf$run_time))
      
      # run bootstrap particle filter
      bpf <- sis_bpf(model_config, observations, bpf_config)
      cat('BPF:', bpf$log_likelihood[T+1], '\n')
      results <- rbind(results, data.frame(slope = slope, 
                                           particles = P, 
                                           repetition = r, 
                                           filter = 'BPF',
                                           loglikelihood = bpf$log_likelihood[T+1], 
                                           runtime = bpf$run_time))
      
      # run controlled particle filter
      cpf <- sis_cpf(model_config, observations, cpf_config)
      cat('CPF:', cpf$log_likelihood[T+1], '\n')
      results <- rbind(results, data.frame(slope = slope, 
                                           particles = P, 
                                           repetition = r, 
                                           filter = 'CPF',
                                           loglikelihood = cpf$log_likelihood[T+1], 
                                           runtime = cpf$run_time))
      
      # print progress
      cat('Slope:', slope, 'No. of particles:', P, 'Repetition:', r, '\n')
    }
  }
  save('results', file = 'inst/results_hetero_sero.RData')
}





