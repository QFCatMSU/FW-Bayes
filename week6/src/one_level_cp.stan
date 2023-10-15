data {
  int<lower=0> J;
  array[J] real y;
  array[J] real sigma;
}
parameters {
  real mu;
  real<lower=0> tau;
  array[J] real theta;
}
model {
  mu ~ normal(0, 5);
  tau ~ cauchy(0, 2.5);
  theta ~ normal(mu, tau);
  y ~ normal(theta, sigma);
}