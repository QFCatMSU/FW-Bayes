data {
  int<lower=0> n_obs; 
  vector[n_obs] C; // kills per wolf 
  vector[n_obs] N; // moose density
}
parameters {
  real<lower = 0, upper = 1> a; 
  real<lower = 0> h; 
  real<lower = 0> sigma;  
}
transformed parameters {
  vector[n_obs] C_pred;
  for(i in 1:n_obs){
    C_pred[i] = (a*N[i]) / (1 + a*h*N[i]); 
  }
}
model {
  h ~ gamma(0.01, 0.01); 
  sigma ~ normal(0, 10)T[0.001,]; 
  C ~ normal(C_pred, sigma); 
}
