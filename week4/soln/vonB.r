#--------------------------------------------------------------------
# simulate and estimate a von Bertalanffy growth equation:
# l_i = linf*(1 - exp(-vbk*(age_i - t0))) + eps_i
# eps_i ~ N(0, sig2)
# i = observation
# linf = average asymptotic size
# vbk = Brody growth coefficient (rate of ascension from t0 to Linf)
# t0 = hypothetical age at which fish length is zero (Ricker 1975)

#--------------------------------------------------------------------
library(tidyverse)
library(cmdstanr)
library(ggqfc)
library(bayesplot)

# simulate fake data
linf <- 500
vbk <- 0.2
t0 <- -1.5
sig <- 75
age_i <- rep(0:35, 5)
n_fish <- length(age_i)

set.seed(123)
l_obs_i <- rnorm(n_fish, mean = linf * (1 - exp(-vbk * (age_i - t0))), sd = sig)

data <- data.frame(
  length = l_obs_i,
  age = age_i
)

# data[181,] <- c(2312, 5) # heh heh
# data <- data[sample(nrow(data)),]
# saveRDS(data, "week4/data/age_length.rds")
data %>%
  ggplot(aes(x = age, y = length)) +
  geom_point() +
  theme_qfc()

stan_data <- list(
  n_obs = length(l_obs_i),
  lengths = l_obs_i / 10,
  ages = age_i
)

inits <- function() {
  list(
    vbk = jitter(0.2, 0.1),
    linf = jitter(55, 1),
    t0 = jitter(0, 0.1),
    sd_len = jitter(sig / 10, 1)
  )
}

mod <- cmdstan_model("week4/soln_files/vb.stan")

fit <- mod$sample(
  data = stan_data,
  init = inits,
  set.seed(71),
  chains = 4,
  parallel_chains = 4,
  iter_warmup = 1000,
  iter_sampling = 1000
)

fit$cmdstan_diagnose()

posterior <- fit$draws(format = "df") # extract draws x variables data frame

np <- nuts_params(fit) # get the sampler parameters

mcmc_pairs(posterior, pars = c("vbk", "linf", "t0", "sd_len"), np = np)
mcmc_trace(posterior, pars = c("vbk", "linf", "t0", "sd_len"), np = np) 
