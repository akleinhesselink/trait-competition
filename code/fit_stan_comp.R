rm(list = ls())

library(rstan)

options(mc.cores = parallel::detectCores())
rstan_options(auto_write = T)

load( file = 'output/simulated_data.rda')

my_dat <- 
  list(N = nrow(exp_grid), 
       s1 = exp_grid$s1, 
       s2 = exp_grid$s2, 
       p = as.numeric(exp_grid$p),
       y = exp_grid$y, 
       x = exp_grid$x, 
       d1 = exp_grid$d1, 
       D = exp_grid$D, 
       m = exp_grid$m, 
       intra = as.numeric( exp_grid$intra )) 

my_dat$p[ is.na( my_dat$p) ] <- 0

my_dat$S <- length(unique(my_dat$s1))
my_dat$P <- max(unique(my_dat$p))

out <- stan(file = 'code/simpler_mod2.stan', 
            pars = c('u', 'phi', 'mu'), 
            data = my_dat, 
            chains = 4, 
            iter = 2000, 
            cores = 4, 
            thin = 10)

test <- summary(out)

rn <- row.names(test$summary)

mu_est <- test$summary[grep('^mu', rn), 1]

plot(my_dat$y, mu_est, type = 'n')

for(i in 1:my_dat$N){ 
  points(my_dat$y[i], mu_est[i], col = my_dat$s1[i] )  
}

abline(0,1)

u_est <- test$summary[grep('^u', rn), 1]

u_est2 <- test$summary[grep('^u', rn), 6]

plot( u$u1, u_est) 
abline(0, 1)

phi_est <- test$summary[ grep('phi', rn), 1]

plot( phi$phi, phi_est )
abline(0,1)


my_dat <- 
  list( N = nrow( exp_grid), 
      K1 = ncol(X11), 
      K2 = ncol(X22), 
      x = exp_grid$x, 
      y = exp_grid$y,
      d1 = exp_grid$d1, 
      D = exp_grid$D, 
      m = exp_grid$m, 
      X11 = X11, 
      X12 = X12, 
      X21 = X21, 
      X22 = X22 ) 

# out <- stan(file = 'code/mod1.stan',
#             data = my_dat, 
#             chains = 1,
#             iter = 100, init = 0.1)


