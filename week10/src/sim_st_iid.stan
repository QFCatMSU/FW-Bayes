/*
  simulate an iid spatio-temporal random field with an exponential kernel 
*/
data {
  int<lower=1> n_sites;
  int<lower=1> n_year;
  matrix[n_sites, n_sites] dist_sites;
  real<lower=0> gp_theta;
  real<lower=0> gp_sigma;
}
generated quantities {
  vector[n_sites] mu_zeros;
  matrix[n_sites, n_sites] SIGMA;
  real<lower=0> gp_sigma_sq;
  matrix[n_sites, n_year] eps_st;  
  // transformations
  gp_sigma_sq = pow(gp_sigma, 2.0);
  // cov matrix between sites
  SIGMA = gp_sigma_sq * exp(-dist_sites / gp_theta);
  for (k in 1:n_sites) {
    mu_zeros[k] = 0;
  }
  for(t in 1:n_year){
    eps_st[,t] = multi_normal_rng(mu_zeros, SIGMA);
  }
}
