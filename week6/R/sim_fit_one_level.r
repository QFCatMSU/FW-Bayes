library("cmdstanr")   

one_level <- cmdstan_model("week6/src/sim_one_level.stan")

# simulate data
sim <- one_level$sample(
  fixed_param = T, # look here 
  iter_warmup = 0, iter_sampling = 1,
  chains = 1, seed = 1
)

# extract it 
y <- as.vector(sim$draws("y", format = "draws_matrix"))
sigma <- as.vector(sim$draws("sigma", format = "draws_matrix"))
theta <- as.vector(sim$draws("theta", format = "draws_matrix"))
mu <- as.vector(sim$draws("mu_print", format = "draws_matrix"))
tau <- as.vector(sim$draws("tau_print", format = "draws_matrix"))

#---------------------------------
# Centered estimation of this model: 

stan_data <- list(
  J = length(y), 
  y = y, 
  sigma = sigma
)

one_level_cp <- cmdstan_model("week6/src/one_level_cp.stan")

# simulate data
fit_cp <- one_level_cp$sample(
  data = stan_data,
  iter_warmup = 1000, iter_sampling = 1000,
  chains = 4, parallel_chains = 4,
  seed = 13, 
  refresh = 1000, 
  adapt_delta = 0.99
)

fit_cp$cmdstan_diagnose()

fit_cp$summary("mu")
fit_cp$summary("tau")

one_level_ncp <- cmdstan_model("week6/src/one_level_ncp.stan")

# simulate data
fit_ncp <- one_level_ncp$sample(
  data = stan_data,
  iter_warmup = 1000, iter_sampling = 1000,
  chains = 4, parallel_chains = 4,
  seed = 13, 
  refresh = 1000, 
  adapt_delta = 0.99
)

fit_ncp$summary("mu")
fit_ncp$summary("tau")
