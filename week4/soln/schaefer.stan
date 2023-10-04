data {
  int<lower=0> n_year; 
  array[n_year] real catches; // harvests
  array[n_year] real log_cpue; // log(catch per unit effort)
}
parameters{
  real<lower=-1> logK; 
  real<lower=-5> logr; 
  real<lower=1, upper =10> iq; // inverse q (catchability)
  real<lower = 50> itau2; 
}
transformed parameters{
  real tau2; 
  real q; 
  real K; 
  real r;
  
  tau2 = 1 / itau2; 
  q = 1 / iq;
  
  K = exp(logK); 
  r = exp(logr); 
}
model {
  array[n_year] real B; 
  array[n_year] real ypred;
  real temp; 
  
 // priors
 logr ~ normal(-1.38, 0.51); 
 iq ~ gamma(0.001, 0.001); 
 itau2 ~ gamma(1.708603, 0.008613854); 
 // implicitly logK ~ uniform(-1, inf)
 
 B[1] = K; 
 ypred[1] = log(B[1]) + log(q); 
 
 for (i in 2:n_year) {
    temp = B[i - 1] + r * B[i - 1] * (1 - B[i - 1] / K) - catches[i-1];
    if (temp < .001) { 
      B[i] = 1 / (2 - temp / .001);
    } else {
      B[i] = temp;
    }
    ypred[i] = log(B[i]) + log(q);
  }
  log_cpue ~ normal(ypred, sqrt(tau2)); 
}
generated quantities {
   array[n_year] real B_med;          // median B | parameter draws
   array[n_year] real ypred2;         // deterministic log cpue 
   array[n_year] real log_cpue_pred;  // predicted log cpue 
   real MSY;                          // maximum sustainable yield 
   real EMSY;                         // effort that gets you MSY
   
   // analytical solutions for equilibrium MSY, EMSY:
   MSY = r*K/4; 
   EMSY = r/(2*q); 
   
   // initialize 
   B_med[1] = K;
   ypred2[1] = log(B_med[1]) + log(q);

   // draw a random deviate for time slice 1: 
   log_cpue_pred[1] = normal_rng(ypred2[1], sqrt(tau2)); 

   // Now do that for remaining time slices:
   for(i in 2:n_year){
    B_med[i] = (B_med[i - 1] + r * B_med[i - 1] * (1 - B_med[i - 1] / K) - catches[i-1]);
    ypred2[i] = log(B_med[i]) + log(q); 
    log_cpue_pred[i] = normal_rng(ypred2[i], sqrt(tau2)); 
  } 
}
