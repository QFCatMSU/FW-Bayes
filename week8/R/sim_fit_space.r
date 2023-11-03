#-------------------------------------------------------------------------------------
# Introduction to spatial models
library(tidyverse)
library(cmdstanr)
library(bayesplot)
library(tidybayes)

set.seed(1)
n_site <- 30 # number of sampling locations 
# simulate random x,y site locations:
g <- data.frame(
  easting = runif(n_site, 0, 10),
  northing = runif(n_site, 0, 10)
)

Locs <- unique(g)
dist_sites <- as.matrix(dist(Locs))

# set up a vector indicating the site index identifiers
n_per_site <- 10
site_id <- rep(1:n_site, each = n_per_site)

# model parameters to simulate
gp_theta <- 1.4 # Gaussian process scale parameter
gp_sigma <- 0.2 # Gaussian process variance / spatial noise parameter
beta0 <- 0 # linear regression intercept
sd_obs <- 0.1 # observation error for likelihood

sim_data <-
  list(
    n_sites = nrow(Locs),
    n = length(site_id), 
    site = site_id,
    dist_sites = dist_sites,
    gp_theta = gp_theta,
    gp_sigma = gp_sigma,
    beta0 = beta0, 
    sd_obs = sd_obs
  )

# compile the model
sim_mod <- cmdstanr::cmdstan_model("week8/src/sim_s_grf.stan")

sim_mod$print()

sim_s <- sim_mod$sample(
  data = sim_data,
  fixed_param = TRUE, 
  iter_warmup = 0, iter_sampling = 1,
  chains = 1, 
  seed = 1
)

# extract the simulated data
y <- as.vector(sim_s$draws("y_obs", format = "draws_matrix"))
out <- data.frame(y,
                  easting = Locs$easting,
                  northing = Locs$northing
)

# plot it to make sure we aren't on drugs
out %>% 
  mutate(y = y - mean(y)) %>%
  ggplot(aes_string(x = "easting", y = "northing", colour = "y")) +
    geom_point(size = 2) +
    scale_color_gradient2()

# fit a stupid model 
fit1 <- lm(out$y ~ 1)

ggplot() +
  geom_qq(aes(sample = rstandard(fit1))) +
  geom_abline(color = "red") +
  coord_fixed() + 
  ggtitle("QQ-Plot")

E <- resid(fit1)

# Examine normality
hist(E)
plot(E, ylim = c(-2, 2), ylab = "Residuals")
abline(h = 0)

# Can we find spatial patterns in the residuals?
out$E <- resid(fit1)
out$my_cex <- 3 * abs(out$E) / max(out$E) + 0.75
out$sign <- as.numeric(out$E >=0) + 1
out$my_pch <- c(1, 16)[out$sign]

# plot it 
p1 <- out %>%
  ggplot(aes(x = easting, y = northing)) +
  geom_point(size = out$my_cex, shape = out$my_pch) +
  ylab("northing (km)") + xlab("easting (km)") + 
  ggtitle("detecting dependency via \nresiduals plotted in space")
p1 
#--------------------------------------------------------------------

# compile the estimation model
mod <- cmdstanr::cmdstan_model("week8/src/est_s_grf.stan")

stan_data <-
  list(
    n_sites = nrow(Locs),
    n = length(site_id),
    site = site_id,
    y = y,
    dist_sites = dist_sites,
    prior_gp_theta = c(0, 2.5),
    prior_gp_sigma = c(0, 1),
    prior_intercept = c(0, 1),
    prior_sd_obs = c(0,1)
  )

inits <- function() {
  list(
    gp_sigma = 0.1,
    gp_theta = 1,
    beta0 = 2,
    sd_obs = 0.4, 
    eps_s = rep(0, n_site)
  )
}

fit <- mod$sample(
  data = stan_data,
  init = inits,
  seed = 1,
  chains = 4,
  parallel_chains = 4,
  iter_warmup = 5000,
  iter_sampling = 5000, 
  step_size = 0.01, 
  max_treedepth = 12,
  adapt_delta = 0.999,
  refresh = 500
)
# saveRDS(fit, "week8/data/fit.rds")
# fit <- readRDS("week8/data/fit.rds")

fit$cmdstan_diagnose()

# The main difference - error around intercept 
# frequentist 
confint(fit1) 

# Bayesian
summarise_draws(fit, ~ quantile(.x, probs = c(0.025, 0.5, 0.975))) %>%
  filter(grepl("beta0", variable))

posterior <- fit$draws(format = "df")
color_scheme_set("gray")

p <- mcmc_hist(posterior, par = "gp_theta")
p + geom_vline(xintercept = gp_theta, lwd = 2) + ggqfc::theme_qfc()

p <- mcmc_hist(posterior, par = "gp_sigma")
p + geom_vline(xintercept = gp_sigma, lwd = 2) + ggqfc::theme_qfc()

p <- mcmc_hist(posterior, par = "eff_range")
p + geom_vline(xintercept = -log(0.05)/gp_theta, lwd = 2) + ggqfc::theme_qfc()

p <- mcmc_hist(posterior, par = "beta0")
p <- p + geom_vline(xintercept = median(posterior$beta0), lwd = 2) + ggqfc::theme_qfc()
p

# plot chain of one parameter
p <- mcmc_trace(posterior, pars = "beta0", np = np) +
     theme_qfc() + theme(text = element_text(size = 20))
p

#---------------
# Non centered -- better for data limited situations 

# compile the estimation model
mod <- cmdstanr::cmdstan_model("week8/src/est_s_grf_ncp.stan")

stan_data <-
  list(
    n_sites = nrow(Locs),
    n = length(site_id),
    site = site_id,
    y = y,
    dist_sites = dist_sites,
    prior_gp_theta = c(0, 2.5),
    prior_gp_sigma = c(0, 1),
    prior_intercept = c(0, 1),
    prior_sd_obs = c(0,1)
  )

inits <- function() {
  list(
    gp_sigma = 0.1,
    gp_theta = 1,
    beta0 = 2,
    sd_obs = 0.4, 
    z_s = rep(0, n_site)
  )
}

fit_ncp <- mod$sample(
  data = stan_data,
  init = inits,
  seed = 1,
  chains = 4,
  parallel_chains = 4,
  iter_warmup = 1000, # fewer iterations 
  iter_sampling = 1000, # fewer iterations
  step_size = 0.01, 
  max_treedepth = 12,
  adapt_delta = 0.999, 
  500
)

# saveRDS(fit, "week8/data/fit.rds")
# fit <- readRDS("week8/data/fit.rds")

fit_ncp$cmdstan_diagnose()
