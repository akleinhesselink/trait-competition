rm(list = ls())

library(sedgwicktraits)
library(sedgwickspecies)

load('data/sedgwick_demographic_parms.rda')

traits <- 
  sedgwicktraits %>% 
  left_join(sedgwick_plants, by = c('species' = 'calflora_taxon')) %>%
  filter( dataset == 'TAPIOCA')

tdist <- dist( traits$`leaf_size(cm2)` , upper = T)
tdiff <- outer( traits$`leaf_size(cm2)`, traits$`leaf_size(cm2)`,  '-')

tdist
tdiff

traits$lambda <- as.numeric (parms$lambda[ traits$prior_code  ] )
traits$s <- as.numeric( parms$s[traits$prior_code])
traits$g <- as.numeric ( parms$g[traits$prior_code] )

demo_traits <- 
  traits %>% 
  filter( !is.na(prior_code)) %>% 
  select( prior_code, lambda, g, s,  `leaf_size(cm2)`:`turgor_loss_point(MPa)`, -seed_mass_data_source, -max_height_data_source, -dataset, -notes) 

demo_traits %>% 
  filter( prior_code != 'CEME') %>% 
  gather( trait, value, `leaf_size(cm2)`:`turgor_loss_point(MPa)`) %>% 
  ggplot( aes( x = value, y = lambda) ) + 
  geom_point() + 
  facet_wrap( ~ trait , scales = 'free') 

demo_traits %>% 
  filter( prior_code != 'CEME') %>% 
  mutate( lambda_by_weight = lambda*`seed_mass(g)`) %>%
  gather( trait, value, `leaf_size(cm2)`:`turgor_loss_point(MPa)`) %>% 
  ggplot( aes( x = value, y = lambda_by_weight) ) + 
  geom_point() + 
  facet_wrap( ~ trait , scales = 'free') 

demo_traits %>% 
  filter( prior_code != 'CEME' ) %>%
  mutate( r0 = lambda*(1 - (1-g)*(1-s)) ) %>% 
  gather( trait, value, `leaf_size(cm2)`:`turgor_loss_point(MPa)`) %>% 
  ggplot( aes( x = value, y = r0) ) + 
  geom_point() + 
  geom_smooth(se = F, method = 'lm') + 
  facet_wrap( ~ trait , scales = 'free') 



