data {
 int nt; 
 vector[nt] y_t;
 int nt_pred; 
}
parameters {
  real<lower=0> sd_o;       // obs
  real<lower=0> sd_p;       // process
  vector[nt] x_t;           // latent state 
  real x0;                  // initial conditions
  real alpha;               // drift
}
model{
  sd_o ~ normal(0, 0.5); 
  sd_p ~ normal(0, 0.5); 
  x0 ~ normal(0, 5); 
  alpha ~ normal(0, 0.25); 
  x_t[1] ~ normal(x0, sd_p); 
  for(i in 2:nt){
    x_t[i] ~ normal(x_t[i-1] + alpha, sd_p); 
  }
  y_t ~ normal(x_t, sd_o); 
}
generated quantities {
  vector[nt + nt_pred] pred_y_t;                           // predicted y_t
  vector[nt + nt_pred] pred_x_t;                           // predicted x_t
  for(i in 1:nt){
    pred_x_t[i] = x_t[i]; 
  } 
  for(i in 1:nt_pred){
    pred_x_t[nt + i] = normal_rng(pred_x_t[nt + i-1] + alpha, sd_p); 
  }
  for(i in 1:(nt+nt_pred)){
    pred_y_t[i] = normal_rng(pred_x_t[i], sd_o); 
  }
}
