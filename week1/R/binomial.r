library(tidyverse)
library(cmdstanr)
library(bayesplot)
# devtools::install_github("QFCatMSU/gg-qfc")
library(ggqfc)

# in the sneak turtle example we released 30 turts,
# went back and detected 2 turtles

n <- 30 # trials

y <- c(
  0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0,
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
)

#------------------------------------------------------------------------------
# compile the .stan model
mod <- cmdstan_model("week1/src/binomial.stan")
mod # look at the model
mod$exe_file() # where did Stan create an executable?

# names in tagged list correspond to the data block in the Stan program
data <- list(n = n, y = y)

# create a function that gives you starting parameter values
inits <- list(p = 0.1) # one way
inits$p
inits$p # always the same

# a smarter way that helps you initialize different chains at different
# starting values for a parameter:
inits <- function() {
  list(
    p = jitter(0.5, amount = 0.05)
  )
}
inits()$p
inits()$p # "jitters" the starting value(s)

fit <- mod$sample(
  data = data,
  init = inits,
  seed = 13, # ensure simulations are reproducible
  chains = 4, # multiple chains
  iter_warmup = 1000, # how long to warm up the chains
  iter_sampling = 1000, # how many samples after warmp
  parallel_chains = 4, # run them in parallel?
  refresh = 500 # print update every 500 iters
)

#------------------------------------------------------------------------------
# check MCMC algorithm diagnostics:
fit$diagnostic_summary()
fit$cmdstan_diagnose()

# plot the chains
posterior <- fit$draws(format = "df") # extract draws x variables data frame
str(posterior)
color_scheme_set("mix-blue-red")
mcmc_trace(posterior, pars = "p") # whole chain
mcmc_trace(posterior, pars = "p", window = c(140, 154)) # plot a chunk

#------------------------------------------------------------------------------
# once numerical diagnostics look reasonable, plot/examine the posterior

# pluck out the summary for parameter p
fit$summary(variables = c("p"), "median", "sd")

posterior %>%
  ggplot(aes(x = p)) +
  geom_histogram(bins = 20) +
  ggtitle("posterior distribution for p") +
  geom_vline(aes(
    xintercept = fit$summary(variables = "p", "median")$median,
    color = "posterior median"
  ), lty = 2, lwd = 2) +
  geom_vline(aes(xintercept = sum(y) / n, color = "analytical MLE"),
    lty = 1, lwd = 2
  ) +
  scale_color_manual(
    name = "",
    values = c(
      `posterior median` = "steelblue",
      `analytical MLE` = "darkorange2"
    )
  ) +
  theme_qfc() +
  theme(legend.position = c(0.85, 0.9))

#------------------------------------------------------------------------------
# posterior predictive distributions (important)
# Our algorithms and posteriors look reasonable at this point, but
# does our estimated model match our observed data?

set.seed(1) # make sure we can replicate the analysis
yrep <- rbinom(nrow(posterior), n, prob = posterior$p)
df <- data.frame(yrep)

df %>%
  ggplot(aes(x = yrep)) +
  geom_histogram() +
  xlab("Simulated `observations` based on our estimated model") +
  ggtitle("Posterior predictive distribution") +
  geom_vline(aes(xintercept = sum(y), color = "observed y"),
    lty = 1, lwd = 2
  ) +
  scale_color_manual(
    name = "",
    values = c(
      `observed y` = "steelblue"
    )
  ) +
  theme_qfc() +
  theme(legend.position = c(0.85, 0.9))

# does our model generate fake data similar to our observed data?

#------------------------------------------------------------------------------
# other fun things we can do with a posterior

# a manager wants to know what is the probability that p lies between 0.1 and
# 0.2 --simply calculate the number of accepted p values between these values,
# divide by number of draws (samples) of the posterior:

sum(posterior$p > 0.1 & posterior$p < 0.2) / nrow(posterior)
