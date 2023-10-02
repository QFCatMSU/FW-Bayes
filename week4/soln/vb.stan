data {
  int<lower=0> n_obs;  
  vector[n_obs] ages; 
  vector[n_obs] lengths; 
}
parameters {
  real<lower = 0> vbk; 
  real<lower = 0> linf; 
  real t0; 
  real<lower = 0> sd_len; 
}
transformed parameters {
  vector[n_obs] l_pred;
  for(i in 1:n_obs){
    l_pred[i] =  linf * (1 - exp(-vbk * (ages[i] - t0))); 
  }
}
model {
  vbk ~ normal(0.2, 1); 
  linf ~ normal(50, 50); 
  t0 ~ normal(0, 5); 
  sd_len ~ exponential(0.25); 
  lengths ~ normal(l_pred, sd_len); 
}
generated quantities {
   // I will leave it to you to do posterior predictive checks 
}

  