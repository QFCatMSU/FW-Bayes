
library(tidyverse)
library(mvtnorm)
library(ggqfc)

# -----------------------
# Playing with multivariate distributions

a <-3.5 
b <- (-1) 
sigma_a <- 1
sigma_b <- 0.5
rho <- (-0.7) 

Mu <-c(a,b)

cov_ab <-sigma_a*sigma_b*rho

# one way to generate the covariance matrix:
SIGMA1 <-matrix(c(sigma_a^2,cov_ab,cov_ab,sigma_b^2),ncol=2)

sigmas <-c(sigma_a,sigma_b) # standard deviations
rho <-matrix(c(1,rho,rho,1),nrow=2) # correlation matrix

# another way matrix multiply to get covariance matrix
SIGMA2 <-diag(sigmas) %*% rho %*% diag(sigmas)

# show that the two versions equal each other
SIGMA1==SIGMA2



nobs <- 1000
# draw some standard normal N(0,1) distributions 
set.seed(1)
x <- rnorm(nobs, mean = 0, sd = 1)
y <- rnorm(nobs, mean = 0, sd = 1)

# covariance between x and y
cov_xy <- mean(x*y) - mean(x)*mean(y)
cov(x,y)
cov_xy

cov_ab <-sigma_a*sigma_b*rho
Sigma <-matrix(c(sigma_a^2,cov_ab,cov_ab,sigma_b^2),ncol=2)

# correlation is just rescaled covariance 
# sd(x)*sd(y) is the largest possible covariance
# sqrt(var(x))*sqrt(var(y))


# covariance between x and itself
cov(x,x)
mean(x^2) - mean(x)^2




# all variables centered on zero
plot(z1, z2, xlim = c(-10, 10), ylim = c(-10, 10))

# shift z1 to the right, X = (z1 + 5, z2)
plot(z1 + 5, z2, xlim = c(-10, 10), ylim = c(-10, 10))

# shift z1 to the left, Y = (z1 - 5, z2)
plot(z1 - 5, z2, xlim = c(-10, 10), ylim = c(-10, 10))

# scale z1 by a factor of 2, Z1 = (2*z1, z2)
plot(2 * z1, z2, xlim = c(-10, 10), ylim = c(-10, 10))

# scale z1 and z2 by a factor of 2, Z2 = (2*z1, 2*z2)
plot(2 * z1, 2 * z2, xlim = c(-10, 10), ylim = c(-10, 10))

# scale z1 and z2 and shift z2 - 2, Z3 = (2*z1, 2*z2 - 2)
plot(2 * z1, 2 * z2 - 2, xlim = c(-10, 10), ylim = c(-10, 10))

# even more shifts, W = (z1 + 4, z2 - z1 - 2)
plot(z1 + 4, z2 - z1 - 2, xlim = c(-10, 10), ylim = c(-10, 10))

# Note that we can write this in linear algebra as: 
Z <- rbind(z1, z2) # z1, z2 are just N(0,1)
A <- rbind(c(1, 0), c(-1, 1))
b <- c(4, -2)
X <- A %*% Z + b # X is now distributed as z1 + 4, z2 - z1 -2
plot(X[1, ], X[2, ], xlim = c(-10, 10), ylim = c(-10, 10))

# get the covariance matrix of W
SIGMA <- A%*%t(A)

# do the same thing as the last version with W in line 32, 41
library(mvtnorm)
print(b)
print(SIGMA)
W <- rmvnorm(n = nobs, mean = b, sigma = SIGMA) 
plot(W, xlim = c(-10, 10), ylim = c(-10, 10))

# convert variance covariance matrix SIGMA to correlation matrix 
D <- diag(1/sqrt(diag(SIGMA)))
D %*% SIGMA %*% D
cov2cor(SIGMA) # check your math 

#----------------
sd_b0 <- 1
sd_b1 <- 1

# std deviations matrix 
sig_mat <- rbind(c(sd_b0, 0), c(0, sd_b1))
sig_mat

# Corrlation matrix 
rho <- 0.6
R <- rbind(c(1, rho), c(rho, 1))

# construct variance-covariance (SIGMA) | sd and correlation matrix 
SIGMA <- sig_mat^2 %*% R
cov2cor(SIGMA)
R




cov_mat <- matrix(c(sd_b0^2, sd_b0 * sd_b1 * rho, sd_b0 * sd_b1 * rho, sd_b1^2), nrow = 2)

R
b <- c(0, 0)
W <- rmvnorm(n = nobs, mean = b, sigma = SIGMA) 


plot(W, xlim = c(-10, 10), ylim = c(-10, 10))

# 

# See also https://willhipson.netlify.app/post/stan-random-slopes/varying_effects_stan/


set.seed(1) # for reproducibility

# unstandardized means and sds
real_mean_y <- 4.25
real_mean_x <- 3.25
real_sd_y <- 0.45
real_sd_x <- 0.54

N <- 30 # number of people
n_days <- 7 # number of days
total_obs <- N * n_days

sigma <- 1 # population sd
beta <- c(0, 0.15) # average intercept and slope
sigma_p <- c(1, 1) # intercept and slope sds
rho <- -0.36 # covariance between intercepts and slopes

cov_mat <- matrix(c(sigma_p[1]^2, sigma_p[1] * sigma_p[2] * rho, sigma_p[1] * sigma_p[2] * rho, sigma_p[2]^2), nrow = 2)
beta_p <- rmvnorm(N, mean = beta, sigma = cov_mat) # participant intercepts and slopes

x <- matrix(c(rep(1, N * n_days), rnorm(N * n_days, 0, 1)), ncol = 2) # model matrix
pid <- rep(1:N, each = n_days) # participant id

sim_dat <- map_dfr(.x = 1:(N * n_days), ~ data.frame(
  mu = x[.x, 1] * beta_p[pid[.x], 1] + x[.x, 2] * beta_p[pid[.x], 2],
  pid = pid[.x],
  x = x[.x, 2]
))

sim_dat$y <- rnorm(210, sim_dat$mu, sigma) # creating observed y from mu and sigma

dat <- sim_dat %>%
  select(-mu) %>% # removing mu
  mutate(
    y = real_mean_y + (y * real_sd_y), # unstandardize
    x = real_mean_x + (x * real_sd_x)
  )
