data {
  int<lower=0> J;
  array[J] real y;
  array[J] real sigma;
}
parameters {
  real mu;
  real<lower=0> tau;
  array[J] real var_theta;
}
transformed parameters {
    array[J] real theta;
  for (j in 1:J) theta[j] = tau * var_theta[j] + mu;
}
model {
  mu ~ normal(0, 5);
  tau ~ cauchy(0, 2.5);
  var_theta ~ normal(0, 1);
  y ~ normal(theta, sigma);
}