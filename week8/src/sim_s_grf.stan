data {
  int<lower=1> n_sites;
  int<lower=1> n;
  array[n] int site;
  matrix[n_sites, n_sites] dist_sites;
  real<lower=0> gp_theta;
  real<lower=0> gp_sigma;
  real beta0;
  real sd_obs; 
}
generated quantities {
  vector[n_sites] mu_zeros;
  matrix[n_sites, n_sites] SIGMA;
  vector[n] y_hat;
  vector[n] y_obs;
  real<lower=0> gp_sigma_sq;
  vector[n_sites] eps_s;  

  // transformations
  gp_sigma_sq = pow(gp_sigma, 2.0);

  // cov matrix between sites
  SIGMA = gp_sigma_sq * exp(-dist_sites / gp_theta);

  for (k in 1:n_sites) {
    mu_zeros[k] = 0;
  }

  //spatial re's first slice
  eps_s = multi_normal_rng(mu_zeros, SIGMA);

  // calculate predicted value of each observation
  for (i in 1:n) {
    y_hat[i] = beta0 + eps_s[site[i]]; // fixed + random effects                        
    y_obs[i] = normal_rng(y_hat[i], sd_obs);  // add likelihood noise
  }
}
