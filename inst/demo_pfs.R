rm(list = ls())
library(agents)
library(ggplot2)
library(ggthemes)
setmytheme()

# SIS model
N <- 200
alpha0 <- rep(0.1, N)
lambda <- 0.25 * runif(N)
# lambda <- rep(runif(1), N)
gamma <- 0.25 * runif(N)
# gamma <- rep(runif(1), N)
network_type <- "full"

# observation model for under-reporting
rho <- 0.8
obs_model <- under_reporting_model(rho)

# # observation model for seroprevalence surveys
# sample_size <- ceiling(0.2 * N)
# obs_model <- seroprevalence_survey_model(sample_size, N)

# model configuration
model_config <- list(N = N, 
                     alpha0 = alpha0, 
                     lambda = lambda,
                     gamma = gamma,
                     network_type = network_type, 
                     rmeasurement = obs_model$rmeasurement,
                     dmeasurement = obs_model$dmeasurement, 
                     infected_support = obs_model$infected_support)

# observations desired
T <- 20
observations <- list(T = T)

# simulate SIS model
sis_sample <- sis_simulate(model_config, observations)
observations$y <- sis_sample$y

# plot simulated sample
sis_df <- data.frame(time = 0:T, infected = sis_sample$I, type = rep("latent", T+1))
sis_df <- rbind(sis_df, data.frame(time = 0:T, infected = sis_sample$y, type = rep("observed", T+1)))
ggplot() + geom_line(data = sis_df, aes(x = time, y = infected, color = type)) + 
  labs(x = "time", y = "no. of infected") + scale_color_colorblind()

# run bootstrap particle filter
bpf_config <- list(P = 2^10, store = FALSE, resampling = systematic_resampling)
bpf <- sis_bpf(model_config, observations, bpf_config)
cat('BPF run-time:', bpf$run_time, 'secs')

# run auxiliary particle filter
apf_config <- list(P = 2^7, store = FALSE, resampling = systematic_resampling, 
                   dpoibin = dpoibin_exact, rcondber = rcondber_exact)
apf <- sis_apf(model_config, observations, apf_config)
cat('APF run-time:', apf$run_time, 'secs')

apf_config <- list(P = 2^7, store = FALSE, resampling = systematic_resampling, 
                   dpoibin = dtranspoi_poibin_approx, rcondber = rcondber_exact)
apf <- sis_apf(model_config, observations, apf_config)
cat('Approximate APF run-time:', apf$run_time, 'secs')

mcmc_iterations <- ceiling(N * log(N))
apf_config <- list(P = 2^7, store = FALSE, resampling = systematic_resampling, 
                   dpoibin = dtranspoi_poibin_approx,
                   rcondber = function(I, a) rcondber_mcmc(I, a, mcmc_iterations))
apf <- sis_apf(model_config, observations, apf_config)
cat('Approximate APF run-time:', apf$run_time, 'secs')

# approximate backward information filter
bif <- sis_bif(model_config, observations) # cost O(TN^3)
cat('BIF run-time:', bif$run_time, 'secs')

bif <- sis_bif(model_config, observations, conditional_expectation = "transpoi") # cost O(TN^2)
cat('Approximate BIF run-time:', bif$run_time, 'secs')

# run controlled particle filter
cpf_config <- list(P = 2^7, store = FALSE, resampling = systematic_resampling, 
                   log_bif = bif$log_bif, dpoibin = dpoibin_exact, rcondber = rcondber_exact)
cpf <- sis_cpf(model_config, observations, cpf_config)
cat('CPF run-time:', cpf$run_time, 'secs')

cpf_config <- list(P = 2^7, store = FALSE, resampling = systematic_resampling, 
                   log_bif = bif$log_bif, dpoibin = dtranspoi_poibin_approx, rcondber = rcondber_exact)
cpf <- sis_cpf(model_config, observations, cpf_config)
cat('Approximate CPF run-time:', cpf$run_time, 'secs')

mcmc_iterations <- ceiling(N * log(N))
cpf_config <- list(P = 2^7, store = FALSE, resampling = systematic_resampling, 
                   log_bif = bif$log_bif, dpoibin = dtranspoi_poibin_approx, 
                   rcondber = function(I, a) rcondber_mcmc(I, a, mcmc_iterations))
cpf <- sis_cpf(model_config, observations, cpf_config)
cat('Approximate CPF run-time:', cpf$run_time, 'secs')

# collect particle filtering results
results_df <- data.frame(time = 0:T, 
                         ess = bpf$ess * 100 / bpf_config$P,  
                         log_likelihood = bpf$log_likelihood, 
                         algorithm = rep("BPF", T+1))
results_df <- rbind(results_df, data.frame(time = 0:T, 
                                           ess = apf$ess * 100 / apf_config$P, 
                                           log_likelihood = apf$log_likelihood, 
                                           algorithm = rep("APF", T+1)))
results_df <- rbind(results_df, data.frame(time = 0:T, 
                                           ess = cpf$ess * 100 / cpf_config$P, 
                                           log_likelihood = cpf$log_likelihood, 
                                           algorithm = rep("CPF", T+1)))

# plot ESS over time
ggplot() + geom_line(data = results_df, aes(x = time, y = ess, color = algorithm)) + 
  labs(x = "time", y = "ESS%") + scale_color_colorblind()

# plot log-likelihood estimates over time
ggplot() + geom_line(data = results_df, aes(x = time, y = log_likelihood, color = algorithm)) + 
  labs(x = "time", y = "log-likelihood") + scale_color_colorblind()

