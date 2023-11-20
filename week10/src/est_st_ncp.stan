data {
  int n_site;
  int n_year; 
  int n;
  array[n] int y_obs; 
  array[n] int site;
  array[n] int year;
  matrix[n_site, n_site] dist_sites;
}
parameters {
  real beta0;
  real<lower=0> gp_theta;
  real<lower=0> gp_sigma;
  real<lower=-1, upper =1> rho; // rho must be in the range [-1, 1]
  matrix[n_site, n_year] z_st;  // non-centering terms 
}
transformed parameters {
  vector<lower=0>[n_site] mu_zeros;
  matrix<lower=0>[n_site, n_site] SIGMA;
  vector[n] y_hat;
  real<lower=0> gp_sigma_sq;
  matrix [n_site, n_site] L;
  matrix[n_site, n_year] eps_st;  
  matrix[n_site, n_year] delta_st;  

  // transformations
  gp_sigma_sq = pow(gp_sigma, 2.0);

  // cov matrix between sites
  SIGMA = gp_sigma_sq * exp(-dist_sites / gp_theta);
  L = cholesky_decompose(SIGMA);

  eps_st = L*z_st; 

  // reconstruct eps_st to deal with stationarity term 
  delta_st[,1] = eps_st[,1]; 
  for(t in 2:n_year){
    delta_st[,t] = rho * delta_st[,t-1] + sqrt(1 - rho^2)*eps_st[,t]; // LOOK HERE
  }
  
  for (k in 1:n_site) {
    mu_zeros[k] = 0;
  }

  // calculate predicted value of each observation
  for (i in 1:n) {
    y_hat[i] = exp(beta0 + delta_st[site[i], year[i]]); // note use of delta_st
  }

}
model {
  // priors:
  gp_theta ~ std_normal();
  gp_sigma ~ std_normal();
  beta0 ~ std_normal();
  to_row_vector(z_st) ~ std_normal(); 
  // rho ~ uniform(-1,1)

  // likelihood 
  y_obs ~ poisson(y_hat);
}
