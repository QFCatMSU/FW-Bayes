set.seed(1)
nobs <- 1000
theta_i <- rnorm(nobs, 0, 1) # local parameter
v <- rnorm(nobs, 0, exp(theta_i)) # global parameter 
plot(theta_i ~ v, main = "visualizing the devil's funnel")

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
    geom_point() 

mod_cp <- cmdstan_model("week6/src/schools_cp.stan")

schools_dat <- list(
  J = 8,                                     # number of schools
  y = c(28,  8, -3,  7, -1,  1, 18, 12),     # standardized test scores
  sigma = c(15, 10, 16, 11,  9, 11, 10, 18)  # standard deviations 
)

fit <- mod_cp$sample(
  data = schools_dat,
  seed = 1, 
  chains = 4,
  iter_warmup = 1000, 
  iter_sampling = 1000, 
  parallel_chains = 4, 
  refresh = 500,   
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
  seed = 1, 
  chains = 4, 
  iter_warmup = 1000, 
  iter_sampling = 1000, 
  parallel_chains = 4, 
  refresh = 500 
)

fit$cmdstan_diagnose()
post <- fit$draws(format = "df") # extract draws x variables data frame
np <- nuts_params(fit)

color_scheme_set("darkgray")

mcmc_pairs(post, np = np, pars = c("mu","tau","theta[1]"),
           off_diag_args = list(size = 0.75))
