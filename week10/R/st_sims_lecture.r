# simulating spatio-temporal random fields
#-----------------------------------------
# IID spatiotemporal random field:

library(tidyverse)
library(cmdstanr)
library(tidybayes)

set.seed(13)
n_site <- 500 # number of sampling locations for each year
# simulate random x,y site locations:
g <- data.frame(
  easting = runif(n_site, 0, 10),
  northing = runif(n_site, 0, 10)
)

locs <- unique(g)
dist_sites <- as.matrix(dist(locs))

# set up a vector indicating the year and site index identifiers
n_year <- 8 # number of years

# model parameters to simulate
gp_theta <- 1 # Gaussian process scale parameter
gp_sigma <- 0.15 # Gaussian process variance / spatial noise parameter

sim_data <-
  list(
    n_sites = nrow(locs),
    n_year = n_year,
    dist_sites = dist_sites,
    gp_theta = gp_theta,
    gp_sigma = gp_sigma
  )

# compile the model
sim_mod <- cmdstan_model("week10/src/sim_st_iid.stan")

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
p1 <- out %>%
  ggplot(aes_string(x = "easting", y = "northing", colour = "eps_st")) +
  geom_point(size = 2) +
  scale_color_gradient2() +
  facet_wrap(~year, nrow = 2)
ggsave("week10/images/st_iid.png", width = 9, height = 5)
#-----------------------------------------
# random walk spatiotemporal random field:

# compile the model
sim_mod <- cmdstan_model("week10/src/sim_st_rw.stan")
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
p2 <- out %>%
  ggplot(aes_string(x = "easting", y = "northing", colour = "eps_st")) +
  geom_point(size = 2) +
  scale_color_gradient2() +
  facet_wrap(~year, nrow = 2)
ggsave("week10/images/st_rw.png", width = 9, height = 5)

#-----------------------------------------
# AR(1) spatiotemporal random field:

# compile the model
sim_mod <- cmdstan_model("week10/src/sim_st_ar1.stan")
sim_data$rho <- 0.5

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
p3 <- out %>%
  ggplot(aes_string(x = "easting", y = "northing", colour = "eps_st")) +
  geom_point(size = 2) +
  scale_color_gradient2() +
  facet_wrap(~year, nrow = 2)
ggsave("week10/images/st_ar1.png", width = 9, height = 5)
