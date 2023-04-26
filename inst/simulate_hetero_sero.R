# Experiments to illustrate impact of heterogeneity
rm(list = ls())
library(agents)
library(ggplot2)
library(ggthemes)

# fix random set
set.seed(21)

# SIS model
N <- 100
alpha0 <- rep(1 / N, N)
link <- function(u) 1 / (1 + exp(-u))
covariates <- rnorm(N)
network_type <- "full"
beta_slope <- 1
beta_lambda <- c(1, beta_slope) 
beta_gamma <- c(-1, beta_slope)
lambda <- link(beta_lambda[1] + beta_lambda[2] * covariates)
gamma <- link(beta_gamma[1] + beta_gamma[2] * covariates)

# observation model for seroprevalence surveys
sample_size <- ceiling(0.2 * N)
obs_model <- seroprevalence_survey_model(sample_size, N)

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

# plot distribution of infection rates
beta_slopes <- c(0.125, 0.25, 0.5, 1.0)
lambda_df <- data.frame()
gamma_df <- data.frame()
for (slope in beta_slopes){
  lambda_df <- rbind(lambda_df, data.frame(slope = factor(rep(slope, N)), 
                                           rate = link(beta_lambda[1] + slope * covariates)))
  gamma_df <- rbind(gamma_df, data.frame(slope = factor(rep(slope, N)), 
                                         rate = link(beta_gamma[1] + slope * covariates)))
}

ggplot(lambda_df, aes(x = rate, fill = slope)) + 
  geom_histogram(position = "identity", alpha = 0.4, aes(y = ..density..), bins = 15) + 
  labs(x = "infection rate", y = "density")

ggplot(gamma_df, aes(x = rate, fill = slope)) + 
  geom_histogram(position = "identity", alpha = 0.4, aes(y = ..density..), bins = 15) + 
  labs(x = "recovery rate", y = "density")

# save observations
save("observations", "covariates", file ="inst/data_hetero_sero.RData")
