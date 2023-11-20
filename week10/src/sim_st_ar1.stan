/*
  simulate an ar(1) spatio-temporal random field with an exponential kernel 
*/
data {
  int<lower=1> n_sites;
  int<lower=1> n_year;
  matrix[n_sites, n_sites] dist_sites;
  real<lower=0> gp_theta;
  real<lower=0> gp_sigma;
  real<lower = - 1, upper = 1> rho; 
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
 
  //initialize spatial re's first slice
  eps_st[,1] = multi_normal_rng(mu_zeros, SIGMA);
  
  //autoregress spatial-temporal re's
  for (t in 2:n_year) {
    eps_st[,t] = rho * eps_st[ ,t-1] + sqrt(1 - rho^2)*multi_normal_rng(mu_zeros, SIGMA);
  }
}
