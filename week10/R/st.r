# spatio-temporal GLMM
library(tidyverse)
library(cmdstanr)
library(bayesplot)
library(tidybayes)

set.seed(13)
n_site <- 30 # number of sampling locations for each year
n_per_site <- 10 # number samples per site
# simulate random x,y site locations:
g <- data.frame(
  easting = runif(n_site, 0, 10),
  northing = runif(n_site, 0, 10)
)

locs <- unique(g)
dist_sites <- as.matrix(dist(locs))

# set up a vector indicating the year and site index identifiers
n_year <- 8 # number of years

year_id <- rep(1:n_year, n_site * n_per_site)
site_id <- rep(1:n_site, each = n_year * n_per_site)

# model parameters to simulate
gp_theta <- 1 # Gaussian process scale parameter
gp_sigma <- 0.15 # Gaussian process variance / spatial noise parameter
rho <- 0.6 # autocorrelation in st field
beta0 <- log(1.75) # linear regression intercept

sim_data <-
  list(
    n_sites = nrow(locs),
    n_year = n_year,
    n = n_site * n_per_site * n_year,
    site = site_id,
    year = year_id,
    dist_sites = dist_sites,
    gp_theta = gp_theta,
    gp_sigma = gp_sigma,
    beta0 = beta0,
    rho = rho
  )

# compile the model
sim_mod <- cmdstan_model("week10/src/sim_st.stan")

sim_mod$print()

sim_s <- sim_mod$sample(
  data = sim_data,
  fixed_param = TRUE, iter_warmup = 0, iter_sampling = 1,
  chains = 1, seed = 1
)

# extract the simulated data
eps_st <- matrix(sim_s$draws("eps_st", format = "draws_matrix"),
  nrow = n_site, ncol = n_year
)

y <- as.matrix(sim_s$draws("eps_st", format = "draws_matrix"))
y <- as.data.frame.table(y, responseName = "eps_st")

out <- data.frame(y,
  easting = rep(locs$easting, n_year),
  northing = rep(locs$northing, n_year),
  year = rep(1:n_year, each = n_site)
)

# plot it to make sure we aren't on drugs

out %>%
  ggplot(aes_string(x = "easting", y = "northing", colour = "eps_st")) +
  geom_point(size = 2) +
  scale_color_gradient2() +
  facet_wrap(~year)

y_obs <- as.vector(sim_s$draws("y_obs", format = "draws_matrix"))
hist(y_obs)

#-------------------------------
# Estimate a fancy model - centered version

est_mod <- cmdstan_model("week10/src/est_st_cp.stan")

stan_data <-
  list(
    n_site = nrow(locs),
    n_year = n_year,
    n = n_site * n_year * n_per_site,
    year = year_id,
    site = site_id,
    y_obs = y_obs,
    dist_sites = dist_sites
  )

inits <- function() {
  list(
    beta0 = 0.1,
    gp_sigma = 0.15,
    gp_theta = 0.15,
    rho = 0.5, 
    eps_st = rep(0, n_site*n_year)
  )
}

fit <- est_mod$sample(
  data = stan_data,
  seed = 1,
  chains = 4,
  iter_warmup = 1000,
  iter_sampling = 1000,
  parallel_chains = 4,
  step_size = 0.01,
  refresh = 500,
  adapt_delta = 0.99
)

fit$cmdstan_diagnose()

fit$print(max_rows = 5)

#-------------------------------
# Estimate a fancy model

est_mod <- cmdstan_model("week10/src/est_st_ncp.stan")

stan_data <-
  list(
    n_site = nrow(locs),
    n_year = n_year,
    n = n_site * n_year * n_per_site,
    year = year_id,
    site = site_id,
    y_obs = y_obs,
    dist_sites = dist_sites
  )

inits <- function() {
  list(
    beta0 = 0.51,
    gp_sigma = 0.13,
    gp_theta = 0.82,
    rho = 0.4, 
    z_st = rep(0, n_site*n_year)
  )
}

fit <- est_mod$sample(
  data = stan_data,
  seed = 1,
  chains = 4,
  iter_warmup = 1000,
  iter_sampling = 1000,
  parallel_chains = 4,
  step_size = 0.01,
  refresh = 500,
  adapt_delta = 0.99
)

fit$cmdstan_diagnose()

fit$print(max_rows = 5)
