#-------------------------------
# intro to state-space models
#-------------------------------

library(tidyverse)
library(ggqfc)
library(cmdstanr)
library(tidybayes)
library(bayesplot)
library(cowplot)

nt <- 15 # number of time steps
x0 <- 3 # intercept
sig_p <- 0.2 # process error
sig_o <- 0.4 # observation error
alpha <- 0.07 # drift

# Simulate predictors
set.seed(2)
x_t <- y_t <- rep(NA, nt)
x_t[1] <- x0
for (t in 2:nt) { # process error
  x_t[t] <- x_t[t - 1] + rnorm(1, mean = alpha, sd = sig_p)
}
for (t in 1:nt) { # observation error
  y_t[t] <- x_t[t] + rnorm(1, mean = 0, sd = sig_o)
}

sim_data <- data.frame(y_t, x_t, t = 1:nt)

sim_data %>%
  ggplot(aes(x = t, y = y_t)) +
  geom_point(alpha = 0.5) +
  geom_line(aes(x = t, y = x_t),
    col = "black", lwd = 0.5,
    linetype = 6
  ) +
  theme_qfc()

# -------------------------
# linear model

fit_lm <- lm(sim_data$y_t ~ sim_data$t)
y_pred_t <- predict(fit_lm, se = TRUE)

sim_data$med <- y_pred_t$fit
sim_data$lower <- y_pred_t$fit - 1.96 * y_pred_t$se.fit
sim_data$upper <- y_pred_t$fit + 1.96 * y_pred_t$se.fit

p1 <- sim_data %>%
  ggplot(aes(x = t, y = y_t)) +
  geom_point(alpha = 0.5) +
  geom_ribbon(aes(ymin = lower, ymax = upper),
    alpha = 0.5, fill = "cadetblue4"
  ) +
  geom_line(aes(x = t, y = med)) +
  geom_line(aes(x = t, y = x_t),
    col = "black", lwd = 0.3,
    linetype = 1
  ) +
  ylim(2, 11.3) +
  ylab(expression(y[t])) +
  geom_ribbon(aes(ymin = lower, ymax = upper),
    alpha = 0.65, fill = "white"
  ) +
  theme_qfc()
p1

# -------------------------
# loess smoother 

fit_loess = loess( sim_data$y_t ~ sim_data$t )
y_pred_t = predict(fit_loess, se=TRUE)

sim_data$med <- y_pred_t$fit
sim_data$lower <- y_pred_t$fit - 1.96 * y_pred_t$se.fit
sim_data$upper <- y_pred_t$fit + 1.96 * y_pred_t$se.fit

p1 <- sim_data %>%
  ggplot(aes(x = t, y = y_t)) +
  geom_point(alpha = 0.5) +
  geom_ribbon(aes(ymin = lower, ymax = upper),
    alpha = 0.5, fill = "cadetblue4"
  ) +
  geom_line(aes(x = t, y = med)) +
  geom_line(aes(x = t, y = x_t),
    col = "black", lwd = 0.3,
    linetype = 1
  ) +
  ylim(2, 11.3) +
  ylab(expression(y[t])) +
  geom_ribbon(aes(ymin = lower, ymax = upper),
    alpha = 0.65, fill = "white"
  ) +
  theme_qfc()
p1

# -------------------------
# generalized additive model

library(mgcv)
fit_gam = gam(sim_data$y_t ~ s(sim_data$t)) # smoother over time
y_pred_t = predict(fit_gam, se=TRUE)

sim_data$med <- y_pred_t$fit
sim_data$lower <- y_pred_t$fit - 1.96 * y_pred_t$se.fit
sim_data$upper <- y_pred_t$fit + 1.96 * y_pred_t$se.fit

p1 <- sim_data %>%
  ggplot(aes(x = t, y = y_t)) +
  geom_point(alpha = 0.5) +
  geom_ribbon(aes(ymin = lower, ymax = upper),
    alpha = 0.5, fill = "cadetblue4"
  ) +
  geom_line(aes(x = t, y = med)) +
  geom_line(aes(x = t, y = x_t),
    col = "black", lwd = 0.3,
    linetype = 1
  ) +
  ylim(2, 11.3) +
  ylab(expression(y[t])) +
  geom_ribbon(aes(ymin = lower, ymax = upper),
    alpha = 0.65, fill = "white"
  ) +
  theme_qfc()
p1

# -------------------------
# Kalman filter / dynamic linear model 

# compile the estimation model
library(cmdstanr)
mod <- cmdstanr::cmdstan_model("week9/src/rw_cp.stan")

stan_data <-
  list(
    nt = nrow(sim_data),
    y_t = sim_data$y_t,
    nt_pred = 12
  )

inits <- function() {
  list(
    sd_o = 0.3,
    sd_p = 0.3,
    x_t = rep(0, nrow(sim_data)),
    x0 = 3,
    alpha = 0.05
  )
}

fit <- mod$sample(
  data = stan_data,
  init = inits,
  seed = 9189,
  chains = 4,
  parallel_chains = 4,
  iter_warmup = 1000,
  iter_sampling = 1000,
  step_size = 0.01,
  adapt_delta = 0.99
)

fit$cmdstan_diagnose()

# projections in Stan

post <-
  summarise_draws(fit, ~ quantile(.x, probs = c(0.025, 0.5, 0.975), na.rm = T)) %>%
  filter(grepl("pred_y_t", variable))

# add extra rows to your data frame:
sim_data <- sim_data[, 1:3] # clean up 
sim_data[nrow(sim_data) + stan_data$nt_pred, ] <- NA
sim_data$t <- 1:nrow(sim_data)

sim_data$lower <- post$`2.5%`[1:(nt + stan_data$nt_pred)]
sim_data$med <- post$`50%`[1:(nt + stan_data$nt_pred)]
sim_data$upper <- post$`97.5%`[1:(nt + stan_data$nt_pred)]

p1 <- sim_data %>%
  ggplot(aes(x = t, y = y_t)) +
  geom_ribbon(aes(ymin = lower, ymax = upper),
    alpha = 0.5, fill = "cadetblue4"
  ) +
  geom_line(aes(x = t, y = x_t),
    col = "black", lwd = 0.3,
    linetype = 1
  ) +
  ylim(2, 11.3) +
  ylab(expression(y[t])) +
  theme_qfc()

post <-
  summarise_draws(fit, ~ quantile(.x, probs = c(0.025, 0.5, 0.975), na.rm = T)) %>%
  filter(grepl("pred_x_t", variable))

sim_data$lower <- post$`2.5%`[1:(nt + stan_data$nt_pred)]
sim_data$med <- post$`50%`[1:(nt + stan_data$nt_pred)]
sim_data$upper <- post$`97.5%`[1:(nt + stan_data$nt_pred)]

# add process error distribution
p1 <- p1 + geom_ribbon(aes(ymin = sim_data$lower, ymax = sim_data$upper),
  alpha = 0.65, fill = "white"
) +
  geom_point(aes(x = sim_data$t, y = sim_data$y_t), alpha = 0.5) +
  geom_line(aes(x = sim_data$t, y = sim_data$x_t),
    col = "black", lwd = 0.5,
    linetype = 6
  ) +
  ylim(2, 12.75) +
  geom_line(aes(x = t, y = med),
    linetype = 1, lwd = 1,
    color = "cadetblue4"
  )
p1

posterior <- fit$draws(format = "df")
np <- nuts_params(fit)

p <- mcmc_pairs(posterior,
  regex_pars = c("sd_o", "sd_p", "x0", "alpha"), np = np
)
p
