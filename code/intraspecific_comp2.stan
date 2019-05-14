data {
  int<lower=1> N;                   // sample size
  int<lower=1> S;                   // number of species 
  int<lower=1, upper=S> species[N]; // species index
  vector<lower=0>[N] x;             // density of competitors 
  vector<lower=0>[N] y;             // response 
}
parameters {
  vector<lower=0>[S] r;             // growth rate for each species
  vector<lower=0>[S] K;             // carrying capacity for each species
  real<lower=0> sigma;              // std dev 
}
transformed parameters{
  vector[N] mu;                     // linear predictor

  for ( i in 1:N) { 
    mu[i] = exp(r[species[i]]*( 1 - (x[i] + 1)/K[species[i]] ) );
  }
}
model {
  // priors
  r ~ cauchy(0,2); 
  sigma ~ cauchy(0,1);
  K ~ cauchy(0,5);
  
  // likelihood
  log(y) ~ normal(log(mu), sigma);
}
generated quantities{ 
  vector[N] y_hat; 
  
  for(i in 1:N)
    y_hat[i] = exp(normal_rng(log(mu[i]), sigma));
}
