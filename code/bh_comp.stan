data {
  int<lower=1> N;                   // sample size
  int<lower=1> S;                   // number of species 
  int<lower=1, upper=S> focal[N];  // focal species index
  int<lower=1, upper=S> comp[N];   // background species index
  vector<lower=0>[N] x;             // density of competitors 
  vector<lower=0>[N] y;             // response (seeds produced)
}
parameters {
  vector<lower=0>[S] lambda;        // gain for each species
  vector<lower=0>[S*S] alpha_vec;   // comp 
  real<lower=0> sigma;              // std dev 
}
transformed parameters{
  vector[N] mu;                     // linear predictor
  matrix<lower=0>[S,S] alpha;
  
  alpha = to_matrix(alpha_vec, S, S ); 
  
  for ( i in 1:N) { 
    mu[i] = lambda[focal[i]]/(1 + alpha[focal[i],comp[i]]*x[i]);
  }
}
model {
  // priors
  lambda ~ cauchy(0, 5); 
  sigma ~ cauchy(0,2);
  alpha_vec ~ cauchy(0,1); 
  
  // likelihood
  log(y) ~ normal(log(mu), sigma);
}
generated quantities{ 
  vector[N] y_hat; 
  
  for(i in 1:N)
    y_hat[i] = exp(normal_rng(log(mu[i]), sigma));
}
