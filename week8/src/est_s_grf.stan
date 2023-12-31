data {
  int<lower=1> n_sites;                // number of sites 
  int<lower=1> n;                      // number of data points
  array[n] int site;                   // site indicator
  vector[n] y;                         // observed data 
  matrix[n_sites, n_sites] dist_sites; // distance matrix 
  vector[2] prior_gp_theta;            // priors 
  vector[2] prior_gp_sigma;
  vector[2] prior_intercept;
  vector[2] prior_sd_obs; 
}
parameters {
  real<lower=0> gp_theta;
  real<lower=0> gp_sigma;
  vector[n_sites] eps_s;  
  real beta0;
  real<lower=0> sd_obs; 
}
transformed parameters {
  matrix<lower=0>[n_sites, n_sites] SIGMA;
  vector[n] y_pred;
  real<lower=0> gp_sigma_sq;
  real delta = 1e-9; 
  
  // transformations
  gp_sigma_sq = pow(gp_sigma, 2.0);  

  // cov matrix between sites
  SIGMA = gp_sigma_sq * exp(-dist_sites / gp_theta);
  
  // Numerical trick to ensure covariance matrix is positive definite:
  SIGMA  = add_diag(SIGMA, delta); // add small offset to every diagonal element of SIGMA  

  // calculate predicted value of each observation
  for (i in 1:n) {
    y_pred[i] = beta0 + eps_s[site[i]]; 
  }
}
model {
  // priors:
  gp_theta ~ normal(prior_gp_theta[1], prior_gp_theta[2]);
  gp_sigma ~ normal(prior_gp_sigma[1], prior_gp_sigma[2]);
  beta0 ~ normal(prior_intercept[1], prior_intercept[2]);
  sd_obs ~ normal(prior_sd_obs[1], prior_sd_obs[2]); 
  
  // spatial re's
  eps_s ~ multi_normal(rep_vector(0, n_sites), SIGMA);
  // likelihood 
  y ~ normal(y_pred, sd_obs);
}
generated quantities {
   real eff_range = 3/gp_theta; // distance at which correlation drops ~ 0.05
}
