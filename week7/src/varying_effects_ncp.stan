data {
  int<lower=1> n;                      // total number of observations
  int<lower=1> n_lakes;                // number of lakes
  int<lower=1> J;                      // number of predictors + intercept
  array[n] int id;                     // lake id vector
  matrix[n, J] X_ij;                   // matrix of predictors
  vector[n] y;                         // observed y data
} 
parameters {
  matrix[J, n_lakes] z_l;              // matrix of intercepts and slopes for each lake
  vector<lower=0>[J] sigma_lake;       // standard deviation for intercept and slope among lakes
  vector[J] betas;                     // "global" intercept and slope 
  cholesky_factor_corr[J] L;           // Cholesky matrix
  real<lower=0> sigma;                 // likelihood/population standard deviation
}
transformed parameters {
  matrix[J, n_lakes] z_lake; 
  matrix[n_lakes, J] betas_lake;
  z_lake = diag_pre_multiply(sigma_lake, L) * z_l; 
  for (i in 1:n_lakes) {
    betas_lake[i, 1] = betas[1] + z_lake[1, i]; // lake specific intercept
    betas_lake[i, 2] = betas[2] + z_lake[2, i]; // lake specific slope 
  }
}
model {
  vector[n] mu;
  L ~ lkj_corr_cholesky(2);        // Cholesky correlation prior
  to_vector(z_l) ~ normal(0, 1);   // non-centering terms
  betas ~ normal(0, 25);           // hyperprior
  sigma_lake ~ exponential(0.01);  // hyperprior
  sigma ~ exponential(0.01);       // likelihood standard deviation prior
  // Calculate predictions | priors and hyperpriors
  for (i in 1:n) {
    mu[i] = betas_lake[id[i], 1] + betas_lake[id[i], 2] * X_ij[i, 2];
  }
  // Likelihood
  y ~ normal(mu, sigma);
}
generated quantities {
  matrix[2, 2] OMEGA;              // Correlation Matrix 
  OMEGA = multiply_lower_tri_self_transpose(L); // L*L^T 
}
