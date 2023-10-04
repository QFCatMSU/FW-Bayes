library(tidyverse)
library(bayesplot)
library(cmdstanr)
data <- read.csv("week4/soln_files/tuna_data.csv")
head(data)
data$year <- 1967:1989

data %>%
  ggplot(aes(x = year, y = CPUE)) +
  geom_point() + 
  geom_line() + 
  labs(y = "Harvest") + 
  scale_x_continuous(breaks = c(1967, 1975, 1983, 1989)) +
  theme_qfc()

data %>%
  ggplot(aes(x = year, y = Catches)) + 
  geom_point() + geom_line()

mp <- seq(1e-12, 1e-3, length.out = 100)
B <- 1 / (2 - tmp / .001)
plot(B ~ tmp)

mod <- cmdstan_model("week4/soln/schaefer.stan")

stan_data <- list(
  n_year = nrow(data),
  catches = data$Catches,
  log_cpue = log(data$CPUE)
)

# note: starting values very important here:
inits <- function() {
  list(
    logr = log(.8), logK = log(279), iq = 5, itau2 = 100
  )
}

fit <- mod$sample(
  data = stan_data,
  init = inits,
  set.seed(71),
  chains = 4,
  parallel_chains = 4,
  iter_warmup = 2000,
  iter_sampling = 2000,
  adapt_delta = 0.999,
  max_treedepth = 12
)

fit$cmdstan_diagnose()
posterior <- fit$draws(format = "df")
np <- nuts_params(fit)

# the things we care about, misery is the river of the world!
mcmc_pairs(posterior,
           pars = c("r", "K", "q", "MSY"),
           np = np, diag_fun = "dens"
)

hist(posterior$tau2)


pp_dist <-
  summarise_draws(fit, ~ quantile(.x, probs = c(0.025, 0.5, 0.975))) %>%
  filter(grepl("log_cpue_pred", variable))
pp_dist$year <- data$year

p <-
  pp_dist %>%
  ggplot(aes(x = year, y = exp(`50%`))) +
  geom_ribbon(aes(ymin = exp(`2.5%`), ymax = exp(`97.5%`)), fill = "#023967", alpha = 0.4) +
  geom_line(aes(x = year, y = exp(`50%`)), color = "white", linetype = 1, lwd = 0.5) +
  geom_point(data = data, aes(x = year, y = CPUE), color = "black") +
  geom_line(
    data = data, aes(x = year, y = CPUE), color = "black",
    linetype = 1, lwd = 0.1
  ) +
  ylab("Catch per unit effort") +
  xlab("Year") +
  ggtitle("Southern Atlantic Albacore Tuna \nfishery dynamics 1967-1989") +
  scale_x_continuous(breaks = c(1967, 1975, 1983, 1989)) +
  theme_qfc() +
  theme(plot.title = element_text(hjust = 0.5))
p

# make it fancy
library(rphylopic)
uuid <- get_uuid(name = "Albacore Tuna")
img <- get_phylopic(uuid = uuid)
p + add_phylopic(
  img = img,
  x = 1987,
  y = 73,
  ysize = 10,
  color = "black"
)
