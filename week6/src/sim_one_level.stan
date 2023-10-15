transformed data {
  real mu;
  real<lower=0> tau;
  real alpha;
  int N;
  mu = 8;
  tau = 3;
  alpha = 10;
  N = 800;
}
generated quantities {
  real mu_print;
  real tau_print;
  vector[N] theta;
  vector[N] sigma;
  vector[N] y;
  mu_print = mu;
  tau_print = tau;
  for (i in 1:N) {
    theta[i] = normal_rng(mu, tau);
    sigma[i] = alpha;
    y[i] = normal_rng(theta[i], sigma[i]);
  }
}