#------------------------------------------------------------------------------
library(tidyverse)
library(cmdstanr)
library(bayesplot)
library(ggqfc)

data <- readRDS("week2/data/linreg.rds")

# compile the .stan model
mod <- cmdstan_model("week2/src/linreg.stan")
mod # look at the model

# names in tagged list correspond to the data block in the Stan program
stan_data <- list(n = nrow(data), y = data$y, x = data$x)

inits <- function() {
  list(
    b0 = jitter(0, amount = 0.05),
    b1 = jitter(0),
    sd = jitter(1)
  )
}

fit <- mod$sample(
  data = stan_data,
  init = inits,
  seed = 13, # ensure simulations are reproducible
  chains = 4, # multiple chains
  iter_warmup = 1000, # how long to warm up the chains
  iter_sampling = 1000, # how many samples after warmp
  parallel_chains = 4, # run them in parallel?
  refresh = 500 # print update every 500 iters
)

# check diagnostics -------------------------------------------

fit$diagnostic_summary()
fit$cmdstan_diagnose()

# extract the posterior and plot the chains:
posterior <- fit$draws(format = "df") # extract draws x variables data frame
str(posterior)

np <- nuts_params(fit) # get the sampler parameters

color_scheme_set("blue")

# all parameters, whole chain:
p <- mcmc_trace(posterior, np = np) +
  theme_qfc()
p

# plot a chunk of chain
p <- mcmc_trace(posterior, pars = "sd", window = c(600, 721), np = np) + 
theme_qfc()
p

# highlight chain 2:
p <- mcmc_trace_highlight(posterior, pars = "sd", highlight = 2) +
  theme_qfc()
p

# pairs plots  - I almost always use this
mcmc_pairs(posterior, pars = c("b0", "b1", "sd"), np = np)

# two parameters only:
p <- mcmc_scatter(posterior, pars = c("b0", "sd"), np = np) +
  theme_qfc()
p

# add a 83% ellipse to it
p + stat_ellipse(level = 0.83, color = "darkorange2", size = 1) +
  theme_qfc()

p + stat_density_2d(color = "darkorange2", size = .5) +
  theme_qfc()

# view it a different way
mcmc_hex(posterior, pars = c("b0", "b1"))

# can also look at rank- I don't use a ton but another visualization
mcmc_rank_overlay(posterior, pars = "sd")

# visualizing divergent transitions
# none here, but this plot shows each iteration as a line connecting
# parameter values, and divergent iterations will show up as red lines
# sometimes this helps you find combinations of parameters that are leading
# to divergent transitions

mcmc_parcoord(posterior,
  pars = c("sd", "b0", "b1"),
  transform = function(x) {
    (x - mean(x)) / sd(x)
  },
  np = np
)

## Have the chains converged to a common distribution?
rhats <- rhat(fit)
mcmc_rhat(rhats) # should all be less than 1.05 as rule of thumb

# number of effective samples
eff <- neff_ratio(fit)
mcmc_neff(eff) # rule of thumb is worry about ratios < 0.1

# autocorrelation
mcmc_acf(posterior)

y_rep <- posterior[grepl("y_rep", names(posterior))]
ppc_hist(y = data$y, yrep = as.matrix(y_rep[1:35, ]))
ppc_dens_overlay(y = data$y, yrep = as.matrix(y_rep[1:35, ]))

#---------------------
# ppcs, another way

y_reps <- y_rep[sample(nrow(y_rep), 9), ]

ind <- sample(9, 1)
y_reps[ind, ] <- as.list(data$y)

yrep_df <- y_reps %>%
  as.data.frame() %>%
  pivot_longer(everything()) %>% # use the long format for plotting
  mutate(name = rep(1:9, each = ncol(y_reps)))

yrep_df %>%
  ggplot() +
  geom_histogram(aes(x = value),
    fill = "steelblue",
    color = "black", binwidth = 1
  ) +
  facet_wrap(~name, nrow = 3) +
  labs(x = "", y = "") +
  scale_y_continuous(breaks = NULL) +
  theme(strip.background = element_blank()) +
  ggtitle("Can you spot the real data?")
