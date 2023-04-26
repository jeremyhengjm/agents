# Experiments to illustrate impact of heterogeneity
rm(list = ls())
library(ggplot2)
library(ggthemes)
library(plyr)
setmytheme()

# load results
load('inst/results_hetero_under.RData')

# summarise results
summary <- ddply(results, c('slope', 'particles', 'filter'), summarise, 
                 elbo = mean(loglikelihood), variance = var(loglikelihood), cost = mean(runtime))
summary['asymptotic_var'] <- summary$particles * summary$variance

# plot asymptotic variance
ggplot(summary) + geom_line(aes(x = particles, y = asymptotic_var, colour = filter)) + 
  facet_grid(cols = vars(slope)) + scale_x_log10(breaks = c(2^7, 2^8, 2^9, 2^10, 2^11)) + 
  scale_y_log10() + scale_color_colorblind() + labs(x = 'no. of particles', y = 'asymptotic variance')

ggplot(summary[summary$particles == 2^11, ]) + geom_line(aes(x = slope, y = asymptotic_var, colour = filter)) + 
  scale_y_log10() + scale_color_colorblind() + labs(x = 'slope', y = 'asymptotic variance')

# plot evidence lower bound
ggplot(summary) + geom_line(aes(x = particles, y = elbo, colour = filter)) + 
  facet_grid(cols = vars(slope)) + scale_x_log10(breaks = c(2^7, 2^8, 2^9, 2^10, 2^11)) + 
  scale_color_colorblind() + labs(x = 'no. of particles', y = 'evidence lower bound')

ggplot(summary[summary$slope == 0.5, ]) + geom_line(aes(x = particles, y = elbo, colour = filter)) +
  scale_x_log10(breaks = c(2^7, 2^8, 2^9, 2^10, 2^11)) + scale_color_colorblind() + labs(x = 'particles', y = 'evidence lower bound')

# plot expected runtime
ggplot(summary) + geom_line(aes(x = particles, y = cost, colour = filter)) + 
  facet_grid(cols = vars(slope)) + scale_x_log10(breaks = c(2^7, 2^8, 2^9, 2^10, 2^11)) + 
  scale_y_log10() + scale_color_colorblind() + labs(x = 'no. of particles', y = 'expected runtime')

