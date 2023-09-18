data {
  int<lower=0> n_year; 
  array[n_year] int counts; // clam counts
  vector[n_year] year; // year covariate 
}
transformed data{
  vector[n_year] year_std;
  real sd_year = sd(year); 
  real mu_year = mean(year); 
  year_std = (year - mu_year) / sd_year; 
}
parameters {
  real b0; 
  real b1;
  real b2; 
  real b3;  
}
transformed parameters {
  vector[n_year] lambdas;
  for(t in 1:n_year){
    lambdas[t] = exp(b0 + b1*year_std[t] + b2*year_std[t]^2 + b3*year_std[t]^3); 
  }
}
model {
  // priors
  b0 ~ normal(0, 10);  
  b1 ~ normal(0, 10);
  b2 ~ normal(0, 10); 
  b3 ~ normal(0, 10); 
  
  // likelihood 
  counts ~ poisson(lambdas);  
}
generated quantities {
  // calculate model diagnostics: 
  vector[n_year] y_rep; 
  vector[n_year] resid_obs;
  vector[n_year] resid_sim;  
  for(i in 1:n_year){
    y_rep[i] = poisson_rng(lambdas[i]); 
    resid_obs[i] = (counts[i] - lambdas[i]) / sqrt(lambdas[i]); 
    resid_sim[i] = (y_rep[i] - lambdas[i]) / sqrt(lambdas[i]); 
  }
}
