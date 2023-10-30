#-------------------------------------------------------------------------------------
# Introduction to spatial models
library(tidyverse)
library(cmdstanr)
library(bayesplot)
library(tidybayes)

set.seed(17)
n_site <- 200 # number of sampling locations 
# simulate random x,y site locations:
g <- data.frame(
  easting = runif(n_site, 0, 10),
  northing = runif(n_site, 0, 10)
)

Locs <- unique(g)
dist_sites <- as.matrix(dist(Locs))

# set up a vector indicating the site index identifiers
site_id <- 1:n_site

# model parameters to simulate
gp_theta <- 1.4 # Gaussian process scale parameter
gp_sigma <- 0.5 # Gaussian process variance / spatial noise parameter
beta0 <- 7 # linear regression intercept
sd_obs <- 0.4
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
    prior_gp_sigma = c(0.5),
    prior_intercept = c(0, 10),
    prior_sd_obs = 0.5
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
  seed = 215,
  chains = 4,
  parallel_chains = 4,
  iter_warmup = 1000,
  iter_sampling = 1000, 
  step_size = 0.01
)
# saveRDS(fit, "week8/data/fit.rds")
fit <- readRDS("week8/data/fit.rds")

fit$cmdstan_diagnose()

# The main difference - error around intercept 
# frequentist 
confint(fit1) 

# Bayesian
summarise_draws(fit, ~ quantile(.x, probs = c(0.025, 0.5, 0.975))) %>%
  filter(grepl("beta0", variable))
