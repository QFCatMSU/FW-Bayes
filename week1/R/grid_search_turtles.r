library(tidyverse)
library(ggqfc)

y <- c(
  0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0,
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
) # sucesses (turtle detections )

n <- length(y) # trials

#---------------------------------------------------------------------------
grid_size <- 10

# p ~ uniform(0, 1):
p_grid = seq(from = 1e-3, to = 0.999, length.out = grid_size)
prior_p <- dunif(p_grid, 0, 1) # note, the same as rep(1, grid_size)

my_df <- data.frame(prior_p, p_grid)
p <- my_df %>%
  ggplot(aes(x = p_grid, y = prior_p)) +
  geom_line() +
  ggtitle("uniform prior distribution for p") +
  xlab("value of p") +
  ylab("probability") +
  theme_qfc()
p

# build a posterio with p values, log(p), log(likelihood), log(posterior) prob
posterior <- data.frame(p = p_grid, log_prior = log(prior_p))

# for each p_grid value, calculate likelihood
posterior$log_lik <- dbinom(sum(y), n, p_grid, TRUE)

# compute unstandardized log posterior-product of likelihood and prior
posterior$log_unstd_posterior <- posterior$log_lik + posterior$log_prior
posterior$unstd_posterior <- exp(posterior$log_unstd_posterior)

# standardize the posterior,so it sums to 1
posterior$posterior_prob <- posterior$unstd_posterior /
  sum(posterior$unstd_posterior)

MAP <- posterior$p[which.max(posterior$posterior_prob)] # best posterior value

posterior %>%
  ggplot(aes(x = p, y = posterior_prob)) +
  geom_point() +
  geom_line() +
  ggtitle(paste0("Grid size = ", grid_size, " points")) +
  xlab(expression(detection ~ probability ~ p)) +
  ylab("posterior probability") +
  geom_vline(aes(
    xintercept = MAP,
    color = "grid estimate"
  ), lty = 2, lwd = 2) +
  geom_vline(aes(xintercept = sum(y) / n, color = "analytical MLE"),
    lty = 1, lwd = 2, alpha = 0.6
  ) +
  scale_color_manual(
    name = "",
    values = c(
      `grid estimate` = "steelblue",
      `analytical MLE` = "darkorange2"
    )
  ) +
  theme_qfc() +
  theme(legend.position = c(0.85, 0.9))
