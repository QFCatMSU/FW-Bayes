data{ 
  int<lower=0> n; // number of observations 
  vector[n] y;    // vector of responses 
  vector[n] x;    // covariate x
}
parameters {
  real b0; 
  real b1;      
  real<lower = 0> sd;
}
model {
  // priors
  b0 ~ normal(0, 10);  
  b1 ~ normal(0, 10);
  sd ~ normal(0, 10);
  
  // likelihood - many ways to do it
  // one way (vectorized, dropping constant, additive terms)
  y ~ normal(b0 + b1*x, sd);  
}
generated quantities {
 // replications for the posterior predictive distributions
 array[n] real y_rep = normal_rng(b0 + b1*x, sd);
}
