// this is a comment 
data {
  // data section, where we read in a
  int<lower=0> n; // number of observations 
  array[n] int y; // vector of binary observations
}
parameters {
  // parameters, where we declare *estimated* parameters
  real<lower=0, upper=1> p;
}
model {
  // model section, where we define prior(s) and likelihood(s) 
  p ~ beta(1, 1); // prior
  y ~ bernoulli(p); // likelihood  

  // NOTE: stan sets up the ln posterior in the background for you:
  // ln(posterior) = ln(priors)+ln(l)
}


