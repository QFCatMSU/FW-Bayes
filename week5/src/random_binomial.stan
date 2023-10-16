data {
  int<lower=0> N;          // number of trials
  array[N] int<lower=0> n; // number released
  array[N] int<lower=0> y; // number alive
  real<lower=0> a;         // values for informative priors
  real<lower=0> b;
  real<lower=0> e;
}
parameters {
  real<lower=0, upper=1> theta_mu;         // mean
  real<lower=0> eta;                       // sample size
  array[N] real<lower=0, upper=1> theta_g; // survival probability for each group    
}
transformed parameters {
  real<lower=0> alpha;
  real<lower=0> beta;
  alpha = eta * theta_mu;
  beta  = eta * (1 - theta_mu);
}
model {
  theta_mu ~ beta(a, b);        // hyper-prior 
  eta ~ exponential(e);         // hyper-prior 
  // implicit joint distributions:
  theta_g ~ beta(alpha, beta);  // hyper distribution for survival pr
  y ~ binomial(n, theta_g);     // likelihood 
}

