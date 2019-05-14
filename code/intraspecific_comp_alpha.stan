data {
  int<lower=1> N;                   // sample size
  int<lower=1> S;                   // number of species 
  int<lower=1, upper=S> species[N]; // species index
  vector<lower=0>[N] x;             // density of competitors 
  vector<lower=0>[N] y;             // response 
}
parameters {
  vector<lower=0>[S] lambda;        // gain for each species
  vector<lower=-0.1>[S] alpha;         // intra comp 
  real<lower=0> sigma;              // std dev 
}
transformed parameters{
  vector[N] mu;                     // linear predictor

  for ( i in 1:N) { 
    mu[i] = lambda[species[i]]/(1 + alpha[species[i]]*x[i]);
  }
}
model {
  // priors
  lambda ~ cauchy(0, 10); 
  sigma ~ cauchy(0,2);
  alpha ~ cauchy(0,1); 
  
  // likelihood
  log(y) ~ normal(log(mu), sigma);
}
generated quantities{ 
  vector[N] y_hat; 
  
  for(i in 1:N)
    y_hat[i] = exp(normal_rng(log(mu[i]), sigma));
}
