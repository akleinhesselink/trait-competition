data {
  int<lower=1> N;                   // sample size
  int<lower=1> S;                   // number of species 
  int<lower=1> P;                   // number of interspecific species pairs 
  int<lower=1,upper=S> s1[N];    // s1 index  
  int<lower=1,upper=S> s2[N];    // s2 index 
  int<lower=0,upper=P> p[N];     // vector of unique species pairs
  vector<lower=0>[N] x;             // density of competitors 
  vector<lower=0>[N] y;             // response 
  vector<lower=1>[N] d1;            // focal species phenology 
  vector<lower=0>[N] m;             // focal species seed mass 
  vector<lower=0, upper=1>[N] D;    // phenology overlap between focal and competition
  int<lower=0,upper=1> intra[N];   // flag intraspecific competition 
}
parameters {
  vector<lower=0>[S] u;             // gain for each species
  vector<lower=0, upper=1>[P] phi;  // pairwise niche overlap 
  real<lower=0> sig;                // std dev.
}
transformed parameters{
  vector[N] mu;                     // linear predictor

  for ( i in 1:N) { 
    real u1; 
    real u2; 
    real PHI; 
    
    u1 = u[s1[i]];
    u2 = u[s2[i]];
    
    if( intra[i]==1 ){ 
      mu[i] = (1/m[i])*u1*d1[i] / (1 + x[i]); // intraspecific 
    }else{
      PHI = phi[ p[i] ]; 
      mu[i] = (1/m[i])*u1*d1[i] / (1 + (u2/u1)*PHI*D[i]*x[i]); //interspecif
    }

  }
  
}
model {
  // priors
  u ~ gamma(1, 4); 
  phi ~ beta(1, 1); 
  sig ~ cauchy(0,1); 
  // likelihood
  y ~ lognormal(mu, sig);
}
generated quantities{ 
  
}
  
