rm(list = ls())

library(rstan)
library(tidyverse)

load('output/stan_dat.rda')

out <- stan(file = 'code/simpler_mod2.stan', 
            pars = c('u', 'phi', 'y_hat'), 
            data = stan_dat, 
            chains = 4, 
            iter = 2000, 
            cores = 4, 
            thin = 5)


y_hat <- summary(out, 'y_hat')

plot(stan_dat$y, y_hat$summary[, 1] )
abline(0,1)

u <- summary(out, 'u')

u_mean <- u$summary[, 1] 

growth_rates <- data.frame( sp1 = levels(inter_pairs$sp1), u = u_mean)

growth_rates %>% 
  ggplot( aes( x = sp1, y = u)) + geom_point() + coord_flip()

phi <- summary(out, 'phi')
phi_mean <- phi$summary[, 1]

overlap <- data.frame( pair = levels( factor(inter_pairs$pair) ), phi = phi_mean)

overlap %>%
  left_join(inter_pairs, by = 'pair') %>% 
  mutate( pair_label = factor(pair, levels = unique(pair[order(phi)]), ordered = T)) %>% 
  distinct(pair_label, phi) %>% 
  filter( str_detect(pair_label, 'AMME') ) %>%
  ggplot( aes( x = pair_label, y = phi)) + 
  geom_point() + 
  coord_flip()

overlap %>%
  left_join(inter_pairs, by = 'pair') %>% 
  mutate( pair_label = factor(pair, levels = unique(pair[order(phi)]), ordered = T)) %>% 
  distinct(pair_label, phi) %>% 
  filter( str_detect(pair_label, 'PLER') ) %>%
  ggplot( aes( x = pair_label, y = phi)) + 
  geom_point() + 
  coord_flip()

