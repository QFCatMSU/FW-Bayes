# libraries 
library(tidyverse)
library(ggqfc)

data <- read.table(file = "week8/data/ph_dat.txt", header = TRUE, dec = ".")
plot(data$Easting, data$Northing)

data$my_factor <- factor(data$Forested,
                        levels = c(1, 2),
                        labels = c("forested", "not forested"))

data %>%
    ggplot(aes(y = pH, x = Altitude, col = my_factor)) + 
    geom_point() + 
    geom_smooth(method = "lm")

data %>%
    ggplot(aes(y = pH, x = SDI, col = my_factor)) + 
    geom_point() + 
    geom_smooth(method = "lm")

data %>%
  ggplot(aes(x = Easting, y = Northing)) +
  geom_point(shape = 16) 

fit <- lm(pH ~ SDI + SDI:my_factor, data)
summary(fit)
plot(fit)
X_ij  <- model.matrix(fit)

plot(x = data$SDI,
     y = data$pH,
     xlab = "SDI",
     ylab = "pH")
abline(fit, lwd = 5)

# Model validation
# Homogeneity
# E3 <- rstandard(fit) -> better way to calculate residuals
E3 <- resid(fit)
F3 <- fitted(fit)
plot(x = F3, y = E3)
abline(h = 0, v = 0)

# Normality
hist(E3)

# Independence due to model misfit
plot(x = data$SDI, y = E3)
abline(h = 0, v = 0)

# Can we find spatial patterns in the residuals?
data$E3 <- resid(fit)
data$my_cex <- 3 * abs(data$E3) / max(data$E3) + 0.75
data$sign <- as.numeric(data$E3 >=0) + 1
data$my_pch <- c(1, 16)[data$sign]

# Convert Easting / Northing to km
data$easting_km <- data$Easting/1000
data$northing_km <- data$Northing/1000
data %>%
  ggplot(aes(x = easting_km, y = northing_km)) +
  geom_point(size = data$my_cex, shape = data$my_pch) +
  ylab("northing (km)") + xlab("easting (km)")

#-------------------------
library(cmdstanr)

data <- data[1:100,]

# Convert Easting / Northing to km
data$easting_km <- data$Easting/1000
data$northing_km <- data$Northing/1000

mod <- cmdstan_model("week8/src/space.stan")
data$SDI_std <- scale(data$SDI)
X_ij <- model.matrix(~SDI_std + SDI_std:my_factor, data)
stan_data <- list(
  n_sites = nrow(unique(data[,c('northing_km','easting_km')])),
  n_obs = nrow(data),
  site = 1:nrow(data), # each site unique
  y = data$pH, 
  X_ij = as.matrix(X_ij),
  J = ncol(X_ij),
  dist_sites = as.matrix(dist(data[,c('easting_km', 'northing_km')]))
)

inits <- function() {
  list(
    gp_theta = 0.2,
    gp_sigma = 0.3,
    eps_s = rep(0, stan_data$n_sites),
    betas = rep(0, stan_data$J),
    sd_obs = 0.11
  )
}

fit <- mod$sample(
  data = stan_data,
  init = inits,
  seed = 1,
  chains = 4,
  parallel_chains = 4,
  iter_warmup = 1000,
  iter_sampling = 1000, 
  step_size = 0.01,
  max_treedepth = 12,
  adapt_delta = 0.99
)

saveRDS(fit, "week8/data/space_fit.rds")

fit$cmdstan_diagnose()
fit$summary()
