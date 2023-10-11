nobs <- 10000
x <- rnorm(nobs, 0, 3)
v <- rnorm(nobs, 0, exp(x))



plot(x, v)

# reinstall the gg-qfc package (updates)
install.packages("devtools")
devtools::install_github("QFCatMSU/gg-qfc")
library(ggqfc)

pp_roll() # check that it worked 

library("bayesplot")
library("ggplot2")
library("cmdstanr")   
library("tidyverse")

schools_dat <- data.frame(
  J = 8,                                     # number of schools
  y = c(28,  8, -3,  7, -1,  1, 18, 12),     # standardized test scores
  sigma = c(15, 10, 16, 11,  9, 11, 10, 18)  # standard deviations 
)

schools_dat %>%
    ggplot(aes(x = y, y = 1:J)) + 
    geom_point() + 
    pp_roll()

mod_cp <- cmdstan_model("week6/src/schools_cp.stan")

schools_dat <- list(
  J = 8,                                     # number of schools
  y = c(28,  8, -3,  7, -1,  1, 18, 12),     # standardized test scores
  sigma = c(15, 10, 16, 11,  9, 11, 10, 18)  # standard deviations 
)

fit <- mod_cp$sample(
  data = schools_dat,
  seed = 1, # ensure simulations are reproducible
  chains = 4, # multiple chains
  iter_warmup = 1000, # how long to warm up the chains
  iter_sampling = 1000, # how many samples after warmp
  parallel_chains = 4, # run them in parallel?
  refresh = 500,  # print update every 500 iters, 
  adapt_delta = 0.999
)

fit$cmdstan_diagnose()
post <- fit$draws(format = "df") # extract draws x variables data frame
np <- nuts_params(fit)

color_scheme_set("darkgray")
div_style <- scatter_style_np(div_color = "red",  div_size = 4)

mcmc_pairs(post, np = np, pars = c("mu","tau","theta[1]"),
           off_diag_args = list(size = 0.75), 
           np_style = pairs_style_np(div_color = "firebrick", div_size = 2))
?mcmc_pairs

#--------------------------------------------------
# non centered parameterization

mod_ncp <- cmdstan_model("week6/src/schools_ncp.stan")

fit <- mod_ncp$sample(
  data = schools_dat,
  seed = 1, # ensure simulations are reproducible
  chains = 4, # multiple chains
  iter_warmup = 1000, # how long to warm up the chains
  iter_sampling = 1000, # how many samples after warmp
  parallel_chains = 4, # run them in parallel?
  refresh = 500 # print update every 500 iters
)

fit$cmdstan_diagnose()
post <- fit$draws(format = "df") # extract draws x variables data frame
np <- nuts_params(fit)

color_scheme_set("darkgray")

mcmc_pairs(post, np = np, pars = c("mu","tau","theta[1]"),
           off_diag_args = list(size = 0.75))
