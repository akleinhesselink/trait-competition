rm(list = ls()) 
library(tidyverse)

d <- read_csv(file = 'data/old_data/2013_Neighbor_survey_prepped.csv')

lambdas <- 
  d %>% 
  filter( neighbours_number == 0  |  str_detect( 'LAMBDA', density_plot_name ), sampled == 1) %>% 
  group_by( target) %>% 
  summarise( lambda =  mean ( seeds, na.rm = T), lambda_sd = sd(  seeds, na.rm = T))

test <- 
  d %>% 
  ungroup() %>% 
  filter( neighbours_number > 0 , sampled == 1) %>% 
  left_join(lambdas, by = 'target') %>% 
  mutate( y_prime = lambda/seeds  - 1)  %>% 
  select( sort_order, plot_number, target, y_prime, background, neighbours_number, y_prime )

sp_dat <- 
  test %>% 
  filter( target != 'SIGA' , target != 'CLBO' ) %>%
  split( . , .$target )

comp_coeffs <- lapply( sp_dat, FUN = function(x, ... ) lm( data = x, ... ), 
                       formula = y_prime ~ -1 + neighbours_number:background)

A <- as.matrix( do.call( rbind, lapply(comp_coeffs, coef) ) )

colnames(A) <- str_extract( colnames(A), pattern = '[A-Z]{4}')

A <- A[ , colnames(A) %in% rownames(A) ]

lambdas <- lambdas[ lambdas$target %in% colnames(A), ]

lambda <- lambdas$lambda

names( lambda ) <- lambdas$target

parms <- list( lambda = lambda, A = A )

rates <- read_csv('data/old_data/species_rates.csv') %>%
  mutate( lambda = iconv(lambda, 'UTF-8', 'UTF-8', sub = '-')) %>% 
  mutate( g = iconv (g , 'UTF-8', 'UTF-8', sub = '-')) %>% 
  mutate( s = iconv (s, 'UTF-8', 'UTF-8', sub = '-')) %>% 
  separate(lambda, c('lambda', 'lambda_se'), sep = '-') %>% 
  separate(g, c('g', 'g_se'), sep = '-') %>% 
  separate(s, c('s', 's_se'), sep = '-') %>% 
  separate( species, sep = '\\s', c('genus', 'species'))

rates <- 
  rates %>% 
  mutate( target =  toupper( paste0( str_sub( genus, 1, 2), str_sub(species, 1,2) )))  %>% 
  select( -species, -genus)

rates$target[ rates$target == 'SACO' ]  <- 'SACA'

rates <- rates[ rates$target %in% names(parms$lambda) , ]

parms$g <- rates$g
parms$s <- rates$s

names( parms$g ) <- rates$target
names(parms$s) <- rates$target

save(parms, file = 'data/sedgwick_demographic_parms.rda')

