library(tidyverse)
library(cmdstanr)
library(ggqfc)
library(bayesplot)
#-------------------------------------------------------------------------
# can think of beta distribution as a distribution of probabilities
# a = successes
# b = failures
#-------------------------------------------------------------------------
a <- 10 # successes
b <- 10 # failures
curve(dbeta(x, a, b))

#-------------------------------------------------------------------------
# Let's do some math to set informative priors for a beta distribution
# First, recognize that a and b are sort of hard to set priors for, right?

# Do some math:
# mu = a / (a + b) -> prior mean
# eta = a + b -> prior sample size

# So what is a,b in temrs of mu, eta?


#-------------------------------------------------------------------------
# Now that we've done the math...let's construct a prior that has 95%
# of its probability mass between 0.05 and 0.45 because previous
# studies have shown chick survival is approximately in this range


# also, we want a prior that is somewhat diffuse for eta (sample size)
# but places most mass on low numbers, so we'll use the
# exponential distribution:
e <- 1 / 20
curve(dexp(x, rate = e), from = 0, to = 100, col = "blue")

# simulate some prior draws
n <- 1e5

prior_draws <- data.frame(
  mu = rbeta(n, a, b),
  eta = rexp(n, e)
) %>%
  mutate(
    alpha = eta * mu,
    beta = eta * (1 - mu)
  )

prior_draws %>%
  tidyr::gather(parameter, value) %>%
  group_by(parameter) %>%
  summarize(
    lower95 = quantile(value, prob = 0.025),
    median = quantile(value, prob = 0.5),
    upper95 = quantile(value, prob = 0.975)
  )

#---------------------------------------------------
data <- readRDS("week5/data/survival_data.rds")

data %>%
  ggplot(aes(x = group, y = n_alive / n_released)) +
  geom_point() +
  ylab("Apparent survival (alive/released)") +
  theme_qfc()

stan_data <- list(
  N = nrow(data), # data
  y = data$n_alive,
  n = data$n_released,
  a = a, # priors
  b = b,
  e = e
)

mod <- cmdstan_model("week5/src/random_binomial.stan")

inits <- function() {
  list(
    theta_mu = 0.2,
    eta = 20,
    theta_g = rep(0.3, length(unique(data$group)))
  )
}

fit <- mod$sample(
  data = stan_data,
  init = inits,
  seed = 53,
  chains = 4,
  iter_warmup = 5000,
  iter_sampling = 5000,
  parallel_chains = 4,
  refresh = 500
)

fit$cmdstan_diagnose()
fit$summary()

posterior <- fit$draws(format = "df")
np <- nuts_params(fit)

# all parameters, whole chain:
p <- mcmc_trace(posterior,
  pars = "eta",
  regex_pars = c("mu", "alpha", "beta"), np = np
) +
  theme_qfc()
p

# plot the marginal distributions of alpha and beta
posterior %>%
  mcmc_hist(pars = c("alpha", "beta")) + 
  ggtitle("marginal distributions of a,b")

# plot the joint distribution of alpha and beta
posterior %>%
  mcmc_scatter(pars = c("alpha", "beta"), np = np) +
  theme_qfc() + ggtitle("joint distribution of a,b")


posterior %>%
  mcmc_intervals(pars = "eta", regex_pars = c("mu", "alpha", "beta"))

posterior %>%
  mcmc_intervals(regex_pars = c("theta"))

p <-
  posterior %>%
  mcmc_areas(regex_pars = c("mu", "theta"), 
             prob = 0.8) + 
  xlim(0,1.0) + theme_qfc()

p <- p + ggtitle(expression(theta[mu] ~ vs ~ group ~ specific ~ theta[g]))
p

fit$summary() # estimates 
data$n_alive / data$n_released # empirical MLE

# Three things to do in groups
# 1. do a prior sensitivity test 
# 2. What is the probability that survival in group 2 is less than average? 
# 2a. What is the probability that survival in group 2 is less than group 6?
# 3. What is the Pr(survival > 0.5) if we went to a new site, stocked new critters? 
