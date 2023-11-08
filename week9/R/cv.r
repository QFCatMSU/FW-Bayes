#--------------------------------------------------------------------
# introduction to cross validation 
#--------------------------------------------------------------------

library(tidyverse)
library(cmdstanr)
library(ggqfc)
library(loo)

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

data %>%
  ggplot(aes(x = age, y = length)) +
  geom_point() +
  theme_qfc()

#------------------------------------------------------------------------------
# leave one out cross validation approximation
# compute approximate LOO-CV using Pareto smoothed importance sampling (PSIS), 
# a new procedure for regularizing importance weights
# see Vehtari et al. 2019 in the lecture materials/references 
#------------------------------------------------------------------------------

stan_data <- list(
  n_obs = length(l_obs_i),
  lengths = l_obs_i / 10,
  ages = age_i
)

inits <- function() {
  list(
    vbk = 0.2,
    linf = 55,
    t0 = 0,
    sd_len = sig / 10
  )
}

mod <- cmdstan_model("week9/src/cv.stan")

fit <- mod$sample(
  data = stan_data,
  init = inits,
  set.seed(71),
  chains = 4,
  parallel_chains = 4,
  iter_warmup = 1000,
  iter_sampling = 1000,
  step_size = 0.01
)

loo1 <- fit$loo()

inits <- function() {
  list(
    vbk = 0.2,
    linf = 55,
    sd_len = sig / 10
  )
}

mod <- cmdstan_model("week9/src/cv2.stan")

fit2 <- mod$sample(
  data = stan_data,
  init = inits,
  set.seed(71),
  chains = 4,
  parallel_chains = 4,
  iter_warmup = 1000,
  iter_sampling = 1000,
  step_size = 0.01
)

loo2 <- fit2$loo()
loo_compare(loo1, loo2)

#-----------------------------------------------------------
# K fold cross validation in Stan

stan_data <- list(
  n_obs = length(l_obs_i),
  lengths = l_obs_i / 10,
  ages = age_i
)

inits <- function() {
  list(
    vbk = 0.2,
    linf = 55,
    t0 = 0,
    sd_len = sig / 10
  )
}

mod <- cmdstan_model("week9/src/cv.stan")

# prepare a matrix with the number of post-warmup iterations by number of observations:
log_pd_kfold <- matrix(nrow = 4000, ncol = nrow(data))

# split the data into k chunks
K <- 10
data$fold <- sample(x = 1:K, size = nrow(data), replace = TRUE)

# see also kfold_split_random(K = K, N = nrow(data))

for (k in 1:K) {
  stan_data_train <- list(
    n_obs = nrow(data[which(data$fold != k), ]),
    lengths = data[which(data$fold != k), "length"] / 10, # this is scaling length, not cross validation related
    ages = data[which(data$fold != k), "age"]
  )
  stan_data_test <- list(
    n_obs = nrow(data[which(data$fold == k), ]),
    lengths = data[which(data$fold == k), "length"] / 10,  # this is scaling length, not cross validation related
    ages = data[which(data$fold == k), "age"]
  )
  # train
  fit_train <- mod$sample(data = stan_data_train, init = inits, set.seed(71), step_size = 0.01, refresh = 0)
  # test
  fit_test <- mod$generate_quantities(fit_train, data = stan_data_test)
  # store it 
  log_pd_kfold[, data$fold == k] <- fit_test$draws("log_lik", format = "matrix")
}

# calculate the log predictive density across that matrix 
elpd_full <- sum(log(colSums(exp(log_pd_kfold))) - log(nrow(log_pd_kfold)))
elpd_full

# which we can do this way as well:
elpd(log_pd_kfold)

# note the elpd() function is actually using the following:
# elpd <- sum(matrixStats::colLogSumExps(log_pd_kfold) - log(nrow(log_pd_kfold)))

#-----------------------
# fit a different model

inits <- function() {
  list(
    vbk = 0.2,
    linf = 55,
    sd_len = sig / 10
  )
}

mod <- cmdstan_model("week9/src/cv2.stan")

# prepare a matrix with the number of post-warmup iterations by number of observations:
log_pd_kfold <- matrix(nrow = 4000, ncol = nrow(data))

# split the data into k chunks
#NOTE--need to predict on same chunks of data!!!!

for (k in 1:K) {
  stan_data_train <- list(
    n_obs = nrow(data[which(data$fold != k), ]),
    lengths = data[which(data$fold != k), "length"] / 10, # this is scaling length, not cross validation related
    ages = data[which(data$fold != k), "age"]
  )
  stan_data_test <- list(
    n_obs = nrow(data[which(data$fold == k), ]),
    lengths = data[which(data$fold == k), "length"] / 10, # this is scaling length, not cross validation related
    ages = data[which(data$fold == k), "age"]
  )
  # train
  fit_train <- mod$sample(data = stan_data_train, init = inits, set.seed(71), step_size = 0.01, refresh = 0)
  # test
  fit_test <- mod$generate_quantities(fit_train, data = stan_data_test)
  # record 
  log_pd_kfold[, data$fold == k] <- fit_test$draws("log_lik", format = "matrix")
}

# calculate the log predictive density across that matrix 
elpd_reduced <- sum(log(colSums(exp(log_pd_kfold))) - log(nrow(log_pd_kfold)))

# examine which model has the larged expected log posterior density: 
elpd_full
elpd_reduced

