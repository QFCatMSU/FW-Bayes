library(tidyverse)
library(bayesplot)
library(cmdstanr)
data <- read.csv("week4/soln_files/moose.csv")

# C = (a*N) / (1 + a*h*N)
# h is handling time 

data %>%
  ggplot(aes(x = moose_per_wolf, y = kills_per_wolf_per_month)) + 
  xlim(0, 150) + ylim(0, 2) + geom_point()

stan_data <- list(
  n_obs = nrow(data),
  C = data$kills_per_wolf_per_month,
  N = data$moose_per_wolf
)

inits <- function() {
  list(
    a = 0.3, 
    h = 2,
    sigma = 5
  )
}

mod <- cmdstan_model("week4/soln/func_response.stan")

fit <- mod$sample(
  data = stan_data,
  init = inits,
  set.seed(82),
  chains = 4,
  parallel_chains = 4,
  iter_warmup = 1000,
  iter_sampling = 1000, 
  step_size = 0.01
)
fit$cmdstan_diagnose()
posterior <- fit$draws(format = "df") # extract draws x variables data frame
mcmc_pairs(posterior, pars = c("a", "h", "sigma"))

