library(tidyverse)
library(MASS)
library(ellipse)
library(ggqfc)
library(gghighlight)
library(bayesplot)
library(cmdstanr)
# devtools::install_github("rmcelreath/rethinking")
library(rethinking)

# -----------------------
# adventures in covariance

sigma_b0 <- 1
sigma_b1 <- 0.5
rho <- (-0.4)
cov_b0_b1 <- sigma_b0 * sigma_b1 * rho

# one way to generate the covariance matrix - note this is not the way in presentation 
SIGMA1 <- matrix(c(sigma_b0^2, cov_b0_b1, cov_b0_b1, sigma_b1^2), ncol = 2)

sigmas <- c(sigma_b0, sigma_b1) # standard deviations
rho_mat <- matrix(c(1, rho, rho, 1), nrow = 2) # correlation matrix

# another way matrix multiply to get covariance matrix
SIGMA2 <- diag(sigmas) %*% rho_mat %*% diag(sigmas)
# diag(sigmas) %*% L --> diag_pre_multiply() in stan

# show that the two versions equal each other
SIGMA1 == SIGMA2

# playing with Cholesky decompositions

U <- chol(SIGMA1)    # upper triangular matrix
L <- t(chol(SIGMA1)) # lower triangular matrix
L %*% t(L)           # recover SIGMA1 as L*L^T

# thinking about pre_diag_multiply
# diag(sigmas) %*% t(chol(rho_mat))
# SIGMA <- diag(sigmas) %*% rho_mat %*% diag(sigmas)
# t(chol(SIGMA))
#----------------------------------
# Let's create a varying effects simulation

# unstandardized means and sds
n_lakes <- 30 # number of lakes
n_visits <- 5 # number of measurements/years at each lake
total_obs <- n_lakes * n_visits

sigma <- 1 # population (likelihood) sd
betas <- c(14, -0.15) # average intercept and slope
sigmas <- c(2, 1) # intercept and slope sds
rho <- -0.3 # covariance between intercepts and slopes
rho_mat <- matrix(c(1, rho, rho, 1), nrow = 2) # correlation matrix
SIGMA <- diag(sigmas) %*% rho_mat %*% diag(sigmas) # var covar

# draw correlated slopes and intercepts:
set.seed(5621)
vary_effects <- mvrnorm(n_lakes, betas, SIGMA)

b0 <- vary_effects[, 1]
b1 <- vary_effects[, 2]

# visualize it
plot(b0, b1,
  col = "steelblue", xlab = "intercepts (b0)", ylab = "slopes(b1)",
  ylim = c(-3, 3), xlim = c(10, 19), main = "Visualizing correlation between slopes and intercepts"
)
for (l in c(0.1, 0.3, 0.5, 0.8, 0.99)) {
  lines(ellipse(SIGMA, centre = betas, level = l))
}

# simulating the rest of our data
n <- n_visits * n_lakes # total number of observations
visit <- rep(1:n_visits, n_lakes)
lake_id <- rep(1:n_lakes, each = n_visits) # create a lake ID
x <- rnorm(n, 0, 1) # create fake covariate data
mu <- b0[lake_id] + b1[lake_id] * x # mean prediction
y <- rnorm(n, mu, sigma) # add likelihood error to mean prediction

data <- data.frame(y, x, lake_id, visit)

# say bad weather means you couldn't go out 20% of the time:
keep <- rbinom(nrow(data), size = 1, prob = 0.75)
data <- data[which(keep == 1), ]

p1 <-
  data %>%
  ggplot(aes(x = x, y = y, group = lake_id)) +
  geom_point(color = "firebrick", alpha = 0.50) +
  geom_smooth(method = "lm", se = FALSE, color = "black") +
  xlab("standardized juvenile density") +
  ylab("growth rate cm/yr") +
  facet_wrap(~lake_id) +
  theme_qfc()
p1

p2 <- data %>%
  ggplot(aes(x = x, y = y)) +
  geom_point(color = "firebrick", alpha = 0.5) +
  geom_smooth(method = "lm", se = FALSE, color = "black") +
  xlab("standardized juvenile density") +
  ylab("growth rate cm/yr") +
  theme_qfc()
p2

# compare with lm()
lm(data$y ~ data$x)
betas[2]

# now we have data much like you could come across working
# for a resource management agency somewhere...

#----------------------------------

# https://mc-stan.org/docs/2_18/functions-reference/lkj-correlation.html
R <- rethinking::rlkjcorr(1e4, K = 2, eta = 2) # try eta = 1, 2, 4, 10
rethinking::dens(R[, 1, 2], xlab = "correlation")

# Now, let's do some prior predictive checks
set.seed(4)
betas_pp <- rnorm(2, mean = 0, sd = 25) # b0 and b1
sigma_pp <- rexp(1, 0.01) # population error
sigmas_pp <- rexp(2, 0.01) # sigma b0 and sigma b1
Omega <- rethinking::rlkjcorr(n = 1, K = 2, eta = 2) # from McElreath's rethinking package

Sigma_pp <- diag(sigmas_pp) %*% Omega %*% diag(sigmas_pp)
beta_sim <- mvrnorm(n_lakes, betas_pp, Sigma_pp)
b0_pp <- beta_sim[, 1]
b1_pp <- beta_sim[, 2]
visit_pp <- rep(1:n_visits, n_lakes)
lake_id_pp <- rep(1:n_lakes, each = n_visits) # create a lake ID
x_pp <- rnorm(n, 0, 1) # create fake covariate data
mu_pp <- b0_pp[lake_id_pp] + b1_pp[lake_id_pp] * x_pp # mean prediction
y_pp <- rnorm(n, mu_pp, sigma_pp)

pp_data <- data.frame(y = y_pp, x = x_pp, lake_id, visit)

pp1 <- pp_data %>%
  ggplot(aes(x = x, y = y, group = lake_id)) +
  geom_point(color = "firebrick", alpha = 0.50) +
  geom_smooth(method = "lm", se = FALSE, color = "black") +
  xlab("juvenile density") +
  ylab("growth rate cm/yr") +
  facet_wrap(~lake_id) +
  ggtitle("simulated from priors") +
  theme_qfc()
pp1
p1 # compare to original data

pp2 <-
  pp_data %>%
  ggplot(aes(x = x, y = y)) +
  geom_smooth(method = "lm", se = FALSE, color = "black") +
  geom_point(size = 3, color = "firebrick", alpha = 0.5) +
  gghighlight(lake_id == 18,
    use_group_by = FALSE, max_highlight = Inf,
    use_direct_label = FALSE
  ) +
  labs(
    title = "simulated from priors",
    subtitle = "Red points are data from lake 18 \n
       Grey points are data from all other lakes"
  ) +
  xlab("juvenile density") +
  ylab("growth rate cm/yr") +
  theme_qfc()
pp2
p2 # compare to original data

# can do this many times, but general theme is that this seems reasonable
#----------------------------------
# estimate the centered parameterization in Stan

mod_cp <- cmdstan_model("week7/src/varying_effects_cp.stan")

stan_data <- list(
  n = nrow(data),
  n_lakes = length(unique(data$lake_id)),
  J = 2, # intercept + slope
  id = data$lake_id,
  X_ij = matrix(c(rep(1, nrow(data)), data$x), ncol = 2), 
  y = data$y
)

fit <- mod_cp$sample(
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
fit$print(max_rows = 60)
fit$print("betas")
fit$print("OMEGA")
fit$print("sigma")
fit$print("sigma_lake")

fit$print("betas_lake[2,2]")
vary_effects[2,] 

post <- fit$draws(format = "df") # extract draws x variables data frame
np <- nuts_params(fit)

#-------------------------------
# visualizations

color_scheme_set("darkgray")
mcmc_pairs(post,
  np = np, pars = c("betas[1]", "betas[2]", "sigma_lake[1]", "sigma_lake[2]", "sigma"),
  off_diag_args = list(size = 0.75),
  np_style = pairs_style_np(div_color = "firebrick", div_size = 2)
)

#-------------------------------
# plot the main effects
mcmc_areas(post, pars = c("betas[1]", "betas[2]", "sigma_lake[1]", "sigma_lake[2]", "sigma"), prob = 0.80) +
  scale_y_discrete(expand = c(0, 0))

#-------------------------------
# plot the random slopes and intercepts:
# first, pluck out the betas_lake, summarise mean slope and intercept each group
betas_lake <- post[, grep("betas_lake", names(post))]
betas_lake_mu <- betas_lake %>%
  pivot_longer(cols = everything()) %>%
  group_by(name) %>%
  summarise(mean = mean(value))
betas_lake_mu$rename <- NA
betas_lake_mu$rename[grep(",1]", betas_lake_mu$name)] <- "intercept"
betas_lake_mu$rename[grep(",2]", betas_lake_mu$name)] <- "slope"
betas_lake_mu <-
  betas_lake_mu %>%
  spread(rename, mean)

intercept <- betas_lake_mu$intercept[c(TRUE, FALSE)][1:n_lakes]
slope <- betas_lake_mu$slope[c(FALSE, TRUE)][1:n_lakes]
re <- data.frame(lake_id = 1:n_lakes, intercept, slope)
re %>%
  ggplot(aes(x = intercept, y = slope)) +
  geom_point(color = "firebrick", alpha = 0.50) +
  # geom_smooth(method = "lm", se = FALSE, color = "black") +
  ggtitle(expression(negative ~ correlation ~ between ~ B[1[lake]] ~ and ~ B[0[lake]]),
    subtitle = "posterior means"
  ) +
  theme_qfc()
re           # estimated means
vary_effects # truth
#-------------------------------
# visualize correlated slopes and intercepts one more way 
beta0 <- post[, which(names(post) == "betas[1]")]
beta1 <- post[, which(names(post) == "betas[2]")]
my_df <- data.frame(beta0, beta1)
names(my_df) <- c("beta0", "beta1")

seq(-1, 1, length.out = 100) %>% # create a sequence of evenly spaced numbers from -1 to 1
  map_dfr(~ data.frame(y = my_df$beta0 + .x * my_df$beta1,
                       x = .x,
                       int = factor(ntile(my_df$beta0, 3), levels = c(1, 2, 3),
                                    labels = c("[0-33%)", "[33%-66%)", "[66%-100%)")))) %>%
  group_by(x, int) %>%
  summarise(mu = mean(y),
            lower = quantile(y, probs = 0.05),
            upper = quantile(y, probs = 0.9)) %>%
  ungroup() %>%
  ggplot() +  
  geom_ribbon(aes(x = x, ymin = lower, ymax = upper), fill = "cadetblue4", alpha = 0.50) +
  geom_line(aes(x = x, y = mu), color = "black", linewidth = 2) +
  xlab ("standardized density") + ylab("growth rate (cm/yr)") + 
  labs(title = "Relation between standardized density and growth rate\n as a function of intercept",
       subtitle = "Intercept percentile",
       caption = "Shaded area = 90% uncertainty interval") +
  facet_wrap(~int) +
  theme_qfc() + 
  theme(plot.subtitle = element_text(hjust = 0.5))

#-------------------------------
# plot posterior predictive distribution vs. data
pp_dist <-
  summarise_draws(fit, ~ quantile(.x, probs = c(0.025, 0.5, 0.975))) %>%
  filter(grepl("y_rep", variable))

data$lower <- pp_dist$`2.5%`
data$med <- pp_dist$`50%`
data$upper <- pp_dist$`97.5%`

data %>%
  ggplot(aes(x = x, y = y)) +
  geom_ribbon(aes(ymin = lower, ymax = upper),
    alpha = 0.5, fill = "cadetblue4"
  ) +
  geom_line(aes(x = x, y = med),
    linetype = 1, lwd = 0.75,
    color = "black"
  ) +
  geom_point(size = 0.75)

#-------------------------------
vary_effects # true values 
re # estimated means

#-------------------------------
# Noncentered version of the model 

mod_ncp <- cmdstan_model("week7/src/varying_effects_ncp.stan")

fit_ncp <- mod_ncp$sample(
  data = stan_data,
  seed = 1,
  chains = 4,
  iter_warmup = 1000,
  iter_sampling = 1000,
  parallel_chains = 4,
  step_size = 0.01,
  refresh = 500, 
  max_treedepth = 12,
  adapt_delta = 0.99
)

fit_ncp$print("betas")
fit_ncp$print("OMEGA")
post <- fit_ncp$draws(format = "df") # extract draws x variables data frame

#-------------------------------
# prove to yourself that both cp and ncp return same estimates...
betas_lake <- post[, grep("betas_lake", names(post))]
betas_lake_mu <- betas_lake %>%
  pivot_longer(cols = everything()) %>%
  group_by(name) %>%
  summarise(mean = mean(value))
betas_lake_mu$rename <- NA
betas_lake_mu$rename[grep(",1]", betas_lake_mu$name)] <- "intercept"
betas_lake_mu$rename[grep(",2]", betas_lake_mu$name)] <- "slope"
betas_lake_mu <-
  betas_lake_mu %>%
  spread(rename, mean)

intercept <- betas_lake_mu$intercept[c(TRUE, FALSE)][1:n_lakes]
slope <- betas_lake_mu$slope[c(FALSE, TRUE)][1:n_lakes]
re_ncp <- data.frame(lake_id = 1:n_lakes, intercept, slope)

re # centered
re_ncp # noncentered
