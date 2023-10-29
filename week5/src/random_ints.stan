data {
  int<lower=0> n;
  vector[n] weights; 
  array[n] int population; 
  vector[n] lengths; 
}
parameters {
  real mu_alpha; 
  vector[3] alpha_j;
  real<lower=0> sd_alpha; 
  real b1;  
  real<lower=0> sd_obs;
}
model {
  vector[n] w_preds; 
  mu_alpha ~ normal(0, 100);             // hyper prior
  sd_alpha ~ normal(0, 100);             // hyper prior
  alpha_j ~ normal(mu_alpha, sd_alpha);  // hyper distribution
  b1 ~ normal(0, 10);                    // slope 
  sd_obs ~ normal(0, 25);                // likelihood error term 
  for(i in 1:n){
    w_preds[i] = alpha_j[population[i]] + b1*lengths[i];
  }
  weights ~ normal(w_preds, sd_obs);     // likelihood 
}

