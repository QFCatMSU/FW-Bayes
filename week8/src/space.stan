data {
  int n_sites;                         // number of sites 
  int n_obs;                           // number of observations
  int J;                               // number of columns of design matrix 
  array[n_obs] int site;               // site indicator
  vector[n_obs] y;                     // observations
  matrix[n_obs, J] X_ij;               // design matrix of predictors 
  matrix[n_sites, n_sites] dist_sites; // distance among sites 
}
parameters {
  real<lower=0> gp_theta;
  real<lower=0> gp_sigma;
  vector[n_obs] eps_s;
  vector[J] betas; 
  real<lower=0> sd_obs; 
}
transformed parameters {
  matrix<lower=0>[n_sites, n_sites] SIGMA; // covariance matrix 
  vector[n_obs] fixed_effects;
  vector[n_obs] y_pred;
  vector[n_obs] resid; 
  real<lower=0> gp_sigma_sq;

  // transformation gp_sigma
  gp_sigma_sq = pow(gp_sigma, 2.0);

  // set up the covariance between sites
  SIGMA = gp_sigma_sq * exp(-dist_sites / gp_theta);
  
  // calculate predicted value of each observation
  fixed_effects = X_ij * betas;  

  for (i in 1:n_obs) {
    y_pred[i] = fixed_effects[i] + eps_s[site[i]]; 
  }
  resid = y - y_pred; // observed - expected
}
model {
  // priors:
  gp_theta ~ normal(0, 500);
  gp_sigma ~ normal(0, 1);
  betas ~ normal(0, 20);

  // spatial re's
  eps_s ~ multi_normal(rep_vector(0, n_sites), SIGMA);
  
  // likelihood 
  y ~ normal(y_pred, sd_obs);
}

