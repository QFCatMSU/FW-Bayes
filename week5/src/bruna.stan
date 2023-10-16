data {
  int<lower=0> n;
  vector[n] weights; 
  array[n] int population; 
  vector[n] lengths; 
}
parameters {
  real mu_alpha; 
  vector[3] eps_j;  // deviations off of the global mu alpha
  real<lower=0> sd_alpha; 
  real b1;  
  real<lower=0> sd_obs;
}
model {
  vector[n] w_preds; 
  mu_alpha ~ normal(0, 100);             // hyper prior
  sd_alpha ~ normal(0, 100);             // hyper prior
  eps_j ~ normal(0, sd_alpha);         // hyper distribution --NOTE ZERO HERE
  b1 ~ normal(0, 10);                    // slope 
  sd_obs ~ normal(0, 25);                // likelihood error term 
  for(i in 1:n){
    w_preds[i] = mu_alpha +              // global intercept  -- LOOK HERE
                 b1*lengths[i] +         // slope  
                 eps_j[population[i]]; // random effect     -- LOOK HERE
  }
  weights ~ normal(w_preds, sd_obs);     // likelihood 
}


