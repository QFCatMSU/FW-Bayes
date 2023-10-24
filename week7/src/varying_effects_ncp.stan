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
  cholesky_factor_corr[J] L_corr;      // Cholesky factor of correlation matrix
  real<lower=0> sigma;                 // likelihood/population standard deviation
}
transformed parameters {
  matrix[J, n_lakes] z_lake;                            // lake specific devaitions from global intercepts and slopes
  matrix[n_lakes, J] betas_lake;                        // lake specific intercepts and slopes
  z_lake = diag_pre_multiply(sigma_lake, L_corr) * z_l; // diag_pre_multiply returns diag(sigma_lake) * L_corr --> diagonal product of sigma_lake and L_corr
  for (i in 1:n_lakes) {
    betas_lake[i, 1] = betas[1] + z_lake[1, i];         // recover lake specific intercept
    betas_lake[i, 2] = betas[2] + z_lake[2, i];         // recover lake specific slope 
  }
  
  /*
    Unpacking Line 19 and Cholesky factorization:
    
    diag_pre_multiply(sigma_lake, L_corr) multiplies sigma_lake (as a diagonal matrix) 
    by the Cholesky factor of the random effect correlation matrix, which returns
    the Cholesky factor of the *covariance* matrix.
    
    This result is then multiplied by the independent Normal(0,1) random effects, 
    z_l, to recover the correlated random effects z_lake.

    (ノಠ益ಠ)ノ彡 ┻━┻ (angry linear algebra table flip)
  */
}
model {
  vector[n] mu;
  L_corr ~ lkj_corr_cholesky(2);   // Cholesky correlation prior -->  L ~ lkj_corr_cholesky(2.0) implies L * L' = OMEGA ~ lkj_corr(2) 
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
  matrix[2, 2] OMEGA;                                // Correlation Matrix 
  OMEGA = multiply_lower_tri_self_transpose(L_corr); // L_corr*L_corr^T 
}
