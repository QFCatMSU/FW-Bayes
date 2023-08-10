library(tidyverse)
library(ggqfc)
library(cmdstanr)
library(bayesplot)

y <- c(25, 14, 68, 79, 64, 139, 49, 119, 111) # obs
A <- factor(c(1, 1, 1, 2, 2, 2, 3, 3, 3)) # group
X <- c(1, 14, 22, 2, 9, 20, 2, 13, 22) # covariate
my_df <- data.frame(y, A, X)

my_df %>%
  ggplot(aes(X, y, color = A)) + 
  geom_point(pch = 16, size = 3.5) + theme_qfc()

mod <- cmdstan_model("week3/soln_files/ANCOVA.stan")
mod # look at the model

# names in tagged list correspond to the data block in the Stan program
X_ij <- model.matrix(~ A - 1 + X)
stan_data <- list(n_obs = nrow(X_ij), n_col = ncol(X_ij), 
                  y_obs = my_df$y, X_ij = as.matrix(X_ij)
                  )

# write a function to get starting values
inits <- function() {
  list(
    b_j = jitter(rep(0, ncol(X_ij)), amount = 0.5),
    sig = jitter(10, 1)
  )
}

fit <- mod$sample(
  data = stan_data,
  init = inits,
  seed = 1, # ensure simulations are reproducible
  chains = 4, # multiple chains
  iter_warmup = 1000, # how long to warm up the chains
  iter_sampling = 1000, # how many samples after warmp
  parallel_chains = 4, # run them in parallel?
  refresh = 500 # print update every 500 iters
)
fit$diagnostic_summary()
# extract the posterior and plot the chains:
post <- fit$draws(format = "df") # extract draws x variables data frame
str(post)
np <- nuts_params(fit) # get the sampler parameters
color_scheme_set("blue")

# all parameters, whole chain:
p <- mcmc_trace(post, np = np) +
  theme_qfc()
p

fit$summary()

# compare to Frequentist version:alpha
# summary(lm(y ~ A - 1 + X)) 

#----------------------------------------
# Poisson GLM:

# leading parameters: 
n_year <- 40
beta0 <- 3.5576
beta1 <-  -0.0912
beta2 <- 0.0091 
beta3 <- -0.00014
set.seed(1)

year <- 1:n_year 
lambda <- rep(NA, n_year) # true lambda

# calculate systematic component, apply link: 
for(i in 1:n_year){
    lambda[i] = exp(beta0 + beta1*year[i] + beta2*year[i]^2 + beta3*year[i]^3)
}

count <- rep(NA, n_year) # container for observed count
count <- rpois(n_year, lambda) # add random poisson error

my_df <- data.frame(year, count, truth = beta0 + exp(beta0 + beta1*year + beta2*year^2 + beta3*year^3))

p <- my_df %>%
    ggplot(aes(x = year, y = count)) + 
    geom_point() +
    theme_qfc()
p + geom_line(aes(x = year, y = truth), col="steelblue", linetype = 1, lwd = 1.5)

#-----------------------
# estimate it in Stan: 

mod <- cmdstan_model("week3/soln_files/clams.stan")
mod # look at the model

# names in tagged list correspond to the data block in the Stan program
stan_data <- list(n_year = n_year, counts = my_df$count, year = my_df$year)

# write a function to get starting values
inits <- function() {
  list(
    b0 = jitter(0, amount = 0.5),
    b1 = jitter(0, amount = 0.5),
    b2 = jitter(0, amount = 0.5),
    b3 = jitter(0, amount = 0.5)
  )
}

fit <- mod$sample(
  data = stan_data,
  init = inits,
  seed = 1, # ensure simulations are reproducible
  chains = 4, # multiple chains
  iter_warmup = 1000, # how long to warm up the chains
  iter_sampling = 1000, # how many samples after warmp
  parallel_chains = 4, # run them in parallel?
  refresh = 500 # print update every 500 iters
)
fit$cmdstan_diagnose()
fit$summary()
fm <- glm(count ~ year + I(year^2) + I(year^3), family = poisson, data = my_df)
summary(fm)

my_df$lambda_est <- fit$summary("lambdas")$mean
my_df$lower95 <- fit$summary("lambdas")$q5
my_df$upper95 <- fit$summary("lambdas")$q95

p <- my_df %>%
    ggplot(aes(x = year, y = count)) + 
    geom_point() +
    theme_qfc()
p <- p + geom_line(aes(x = year, y = truth), col="steelblue", linetype = 1, lwd = 1.5)
p <- p + geom_line(aes(x = year, y = lambda_est), col="red", linetype = 2, lwd = 1.5)
p + geom_ribbon(aes(ymin = lower95, ymax = upper95), alpha = 0.5)

# note this is uncertainty around MEAN, not posterior predictive distribution.

#----------------------------------
posterior <- fit$draws(format = "df") # extract draws x variables data frame

y_rep <- posterior[grepl("y_rep", names(posterior))]
ppc_hist(y = stan_data$counts, yrep = as.matrix(y_rep[1:35, ]))
ppc_dens_overlay(y = stan_data$counts, yrep = as.matrix(y_rep[1:35, ]))

# Bayesian goodness of fit test 
resid_obs <- posterior[grepl("resid_obs", names(posterior))]
resid_sim <- posterior[grepl("resid_sim", names(posterior))]

# sum of squared pearson residuals for observed data:
ss_obs <- apply(resid_obs, 1, function(x) sum(x^2))

# sum of squared pearson residuals for simulated data:
ss_sim <- apply(resid_sim, 1, function(x) sum(x^2))

mean(ss_sim > ss_obs) # I have seen better, I have seen worse
plot(ss_sim~ ss_obs)
abline(0,1)

# NOTE: if the sum of squared pearson residuals for the 
# observed data is always largr than that of simulated data, 
# indicates that Poisson model is overdispersed 
