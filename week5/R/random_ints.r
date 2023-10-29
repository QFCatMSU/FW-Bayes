library(cmdstanr)
library(tidyverse)
library(bayesplot)
library(ggqfc)

weight <- c(25, 14, 68, 79, 64, 139, 49, 119, 111) # obs
population <- factor(c(1, 1, 1, 2, 2, 2, 3, 3, 3)) # group
length <- c(1, 14, 22, 2, 9, 20, 2, 13, 22) # covariate
my_df <- data.frame(weight, population, length)

my_df %>%
  ggplot(aes(length, weight, color = population)) +
  geom_point(pch = 16, size = 3.5) +
  theme_qfc()

# Frequentist 
library(lme4)
fit_partial <- lmer(my_df$weight ~ 1 + my_df$length + (1|my_df$population)) 
summary(fit_partial)

# now in Bayesian:
mod <- cmdstan_model("week5/src/random_ints.stan")
mod

stan_data <- list(
  n = length(weight),
  weights = weight,
  population = population,
  lengths = length 
)

inits <- function() {
  list(
    mu_alpha = 50,
    sd_alpha = 30,
    alpha_j = rep(0, 3),
    b1 = 10,
    sd_obs = 10
  )
}

# inits for Bruna model - note plotting code may not work 
#inits <- function() {
#  list(
#    mu_alpha = 50,
#    sd_alpha = 30,
#    eps_j = rep(0, 3),
#    b1 = 10,
#    sd_obs = 10
#  )
#}

fit <- mod$sample(
  data = stan_data,
  init = inits,
  seed = 1, 
  chains = 4, 
  iter_warmup = 1000, 
  iter_sampling = 1000, 
  parallel_chains = 4, 
  refresh = 50, 
  adapt_delta = 0.99999, 
  step_size = 1e-3
)

fit$summary()
coef(fit_partial)

fit$cmdstan_diagnose()

# extract the posterior and plot the chains:
posterior <- fit$draws(format = "df") 
color_scheme_set("blue")

np <- nuts_params(fit) # get the sampler parameters

p <- mcmc_trace(posterior, pars = c("mu_alpha", "b1", "sd_obs"), np = np) +
  theme_qfc()
p

# chains for stuff associated with random effects 
p <- mcmc_trace(posterior,  regex_pars = c("alpha_j", "b1", "sd_obs", "sd_alpha"),
                np = np) +
  theme_qfc()
p

# pairs plots  
mcmc_pairs(posterior, regex_pars = c("alpha"), np = np)

# plot posterior distribution vs. Frequentist estimate 
p <- mcmc_hist(posterior, par = "mu_alpha")
p + geom_vline(xintercept = fixef(fit_partial)[1])

