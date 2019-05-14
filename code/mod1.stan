data {
  int<lower=1> N;                   // sample size
  int<lower=1> K1;                  // K1 predictors
  int<lower=1> K2;                  // K2 predictors 
  vector<lower=0>[N] x;                // density of competitors 
  vector<lower=0>[N] y;             // response 
  vector<lower=1>[N] d1;               // focal species phenology 
  vector<lower=0>[N] m;             // focal species seed mass 
  vector<lower=0, upper=1>[N] D;    // phenology overlap between focal and competitor 
  matrix[N,K1] X11;                 // trait matrix for focal species
  matrix[N,K1] X12;                 // trait matrix for competitor 
  matrix[N,K2] X21;                 // trait matrix for focal for niche overlap  
  matrix[N,K2] X22;                 // trait matrix for competitor for niche overlap 
}
parameters {
  vector[K1] Beta1;                 // trait effects on u
  vector<lower=0>[K2] Beta2;        // trait effects on niche overlap 
  real<lower=0> sig;                // individual error 
}
transformed parameters{
  vector[N] u1;                     // gain rate species 1 
  vector[N] u2;                     // gain rate species 2 

  vector<lower=0, upper=1>[N] phi;  // niche overlap
  vector[N] mu;                     // linear predictor

  u1 = exp( X11 * Beta1); 
  u2 = exp( X12 * Beta1); 
  
  for( i in 1:N){ 
    {
      vector[K2] X21_scaled; 
      vector[K2] X22_scaled; 
      real dist; 
      
      X21_scaled = to_vector(X21[i,1:K2]) .* Beta2; // elementwise multiplication 
      X22_scaled = to_vector(X22[i,1:K2]) .* Beta2;  
      
      dist = distance(X21_scaled, X22_scaled); 
      phi[i] = inv_logit( 1/dist ); 
    }
  }
  
  
  for ( i in 1:N) { 
    mu[i] = (1/m[i])*u1[i]*d1[i] / (1 + (u2[i]/u1[i])*phi[i]*D[i]*x[i]);
  } 
  
}
model {
  // model calculations
  // priors
  sig ~ cauchy(0, 1); 
  Beta1 ~ normal(0, 1);   
  Beta2 ~ normal(0, 1);              
  
  // likelihood
  y ~ lognormal(mu, sig);
}
