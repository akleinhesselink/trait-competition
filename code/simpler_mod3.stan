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
  real <lower = 0> sig; 
}
transformed parameters{
  vector[N] mu;                     // linear predictor

  for ( i in 1:N) { 
    real u2; 
    real PHI; 
    
    //u1 = u[s1[i]];
    u2 = u[s2[i]];
    
    if( intra[i]==1 ){ 
      u_pred[i] = (m[i]/d1[i])*(1 + x[i])*y[i]; // intraspecific 
    }else{
      real a; 
      real b; 
      // Use quadratic equation to solve for u1, u2 and phi are estimated
      PHI = phi[ p[i] ]; 
      
      a = d1[i]/m[i];
      b = u2*PHI*D[i]*x[i];
      
      u_pred[i] = (y + sqrt(4*a*b*y + y^2))/(2*a);
    }

  }
  
}
model {
  // priors
  phi ~ beta(1, 1); 
  
  // likelihood
  u ~ normal(u_sp, 1); 
}
generated quantities{ 
  
}
  
