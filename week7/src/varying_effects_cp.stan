data {
  int n;                               // total number of observations
  int n_lakes;                         // number of lakes
  int J;                               // number of predictors + intercept
  array[n] int id;                     // lake id vector
  matrix[n, J] X_ij;                   // matrix of predictors
  array[n] real y;                     // observed y data
} 
parameters { 
  array[n_lakes] vector[J] betas_lake; // intercept and slope for each lake 
  vector<lower=0>[J] sigma_lake;       // sd for intercept and slope among lakes
  vector[J] betas;                     // "global" intercept and slope 
  corr_matrix[J] OMEGA;                // J by J correlation matrix
  real<lower=0> sigma;                 // likelihood/population sigma
}
model {
  vector[n] mu;
  betas ~ normal(0, 25);           // hyper prior
  OMEGA ~ lkj_corr(2);             // hyper prior
  sigma_lake ~ exponential(0.01);  // hyper prior
  sigma ~ exponential(0.01);       // prior
  betas_lake ~ multi_normal(betas, quad_form_diag(OMEGA, sigma_lake)); // hyperdistribution
  // calculate predictions | priors and hyperpriors
  for(i in 1:n) {
    mu[i] = X_ij[i] * (betas_lake[id[i]]); // * is matrix multiplication in this context
  }
  // likelihood to contrast predictions vs. data:
  y ~ normal(mu, sigma);
}
generated quantities {
  array[n_lakes] vector[J] betas_lake_rep; 
  array[n] real y_rep; 
  // simulate random effects 
  for(l in 1:n_lakes){
    betas_lake_rep[l][1:2] = multi_normal_rng(betas, quad_form_diag(OMEGA, sigma_lake));   
  }
  // simulate likelihood noise  
  for(i in 1:n) {
    y_rep[i] = normal_rng(X_ij[i] * (betas_lake_rep[id[i]]), sigma);            
  } 
}
