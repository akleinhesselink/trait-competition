functions {
  real[] competition(real t,
                  real[] y,
                  real[] theta,
                  real[] x_r,
                  int[] x_i) {
    real dydt[2];
    dydt[1] = y[1]*(theta[1]/(1 + (theta[2]/theta[1])*theta[3]*y[2]));
    dydt[2] = y[2]*(theta[2]/(1 + (theta[1]/theta[2])*theta[3]*y[1]));
    
    return dydt;
  }
}
data {
  int<lower=1> T;
  real y0[2];
  real t0;
  real ts[T];
  real theta[3];
  real sigma[2];
}
transformed data {
  real x_r[0];
  int x_i[0];
}
model {
  real y_hat[T,2];
  y_hat = integrate_ode_rk45(competition, y0, t0, ts, theta, x_r, x_i);
  y_hat[T] ~ normal(y_hat[T], sigma);
}
