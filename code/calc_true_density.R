rm(list = l())
library(tidyverse)
library(sedgwickspecies)

demo_pars <- read_csv('data/demo_pars.csv')
dat <- read.csv(file = 'data/old_data/2013_Neighbor_survey_prepped.csv')

germination <- 
  demo_pars %>% 
  select( focal, g) %>% 
  rename( 'USDA_symbol' = focal )

dat$density <- as.numeric(str_extract( dat$density, '\\d+')  )
dat$density[is.na(dat$density)] <- 0
dat$background <- as.character(dat$background)
dat$target <- as.character(dat$target)
dat$background[is.na(dat$background)] <- 'lambda'

#
dat$background[ dat$background == 'NAJA' & dat$sampled == 1]  <- 'NAAT'

codes <- 
  sedgwick_plants %>% 
  ungroup() %>% 
  select(USDA_symbol, prior_code) %>%
  filter( !is.na(prior_code))

dat <- 
  dat %>% 
  left_join(codes, by = c('target' = 'prior_code')) %>% 
  mutate( target = USDA_symbol)  %>%
  select( -USDA_symbol) %>% 
  left_join (codes, by = c('background' = 'prior_code')) %>% 
  mutate( background = USDA_symbol) %>% 
  select( -USDA_symbol)

seed_mass <- 
  sedgwicktraits::sedgwicktraits %>% 
  group_by( USDA_symbol) %>% 
  filter( dataset == 'TAPIOCA') %>% 
  distinct(USDA_symbol, seed_mass)

test <- 
  dat %>% 
  left_join(seed_mass, by = c('background' = 'USDA_symbol')) %>% 
  left_join(germination, by = c('background' = 'USDA_symbol')) %>% 
  mutate( seeds_per_m = (density*(1/seed_mass)), 
          seeds_per_neighborhood = seeds_per_m*((pi*7^2)/10000), 
          expected_neighbors = seeds_per_neighborhood*g)

test %>% 
  ggplot( aes( x = expected_neighbors, y = neighbours_number)) + 
  geom_point() + 
  stat_summary(fun.y = 'mean', geom = 'line', color = 'red') + 
  geom_abline(aes(intercept = 0, slope = 1)) + 
  facet_wrap( ~ background, scale = 'free')

test %>% 
  group_by( background) %>%
  summarise(my_g = mean( neighbours_number/seeds_per_neighborhood, na.rm = T), g = mean(g, na.rm = T))

  
