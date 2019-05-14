rm(list = ls())
library(sedgwickspecies)
library(sedgwicktraits)
library(tidyverse)

dat <- read.csv(file = 'data/old_data/2013_Neighbor_survey_prepped.csv')

dat$density <- as.numeric(str_extract( dat$density, '\\d+')  )
dat$density[is.na(dat$density)] <- 0
dat$background <- as.character(dat$background)
dat$target <- as.character(dat$target)
dat$background[is.na(dat$background)] <- 'lambda'

#
dat$background[ dat$background == 'NAJA' & dat$sampled == 1]  <- 'NAAT'

intraspecific <- 
  dat %>% 
  filter( sampled == 1) %>% 
  mutate( background = ifelse(background == 'lambda', target, background)) %>% 
  filter( target == background) 

intraspecific %>%
  group_by( target) %>% 
  mutate( geom_lambda = exp(mean(log(seeds_uncorrected[neighbours_number == 0 ])) ), 
          mean_lambda = mean(seeds_uncorrected[neighbours_number == 0])) %>% 
  mutate( geom_expected = geom_lambda/(1 + neighbours_number), 
          mean_expected = mean_lambda/(1 + neighbours_number)) %>% 
  ggplot( aes( x = neighbours_number,  y = seeds_uncorrected )) + 
  geom_point() + 
  geom_line(aes( y = geom_expected )) +
  geom_line(aes( y = mean_expected), linetype = 2) + 
  facet_wrap(~target, scales = 'free')

intraspecific %>%
  group_by( target) %>% 
  mutate( geom_lambda = exp(mean(log(seeds_uncorrected[neighbours_number == 0 ])) ), 
          mean_lambda = mean(seeds_uncorrected[neighbours_number == 0])) %>% 
  mutate( geom_expected = geom_lambda/(1 + neighbours_number), 
          mean_expected = mean_lambda/(1 + neighbours_number)) %>% 
  ggplot( aes( x = neighbours_number,  y = seeds_uncorrected )) + 
  geom_point() + 
  geom_line(aes( y = geom_expected )) +
  geom_line(aes( y = mean_expected), linetype = 2) + 
  facet_wrap(~target, scales = 'free') + 
  scale_y_log10()

interspecific <- 
  dat %>% 
  filter( sampled == 1) %>% 
  mutate( background = ifelse(background == 'lambda', target, background))

stan_dat <- 
  list( species_labels = factor(intraspecific$target), 
      x = intraspecific$neighbours_number, 
      y = intraspecific$seeds_uncorrected )

stan_dat$species <- as.numeric( stan_dat$species_labels)
stan_dat$S <- max(stan_dat$species)
stan_dat$N <- length(stan_dat$y)

stan_dat <-
  list( focal_labels = factor(interspecific$target), 
        comp_labels = factor(interspecific$background, levels = levels(factor(interspecific$target))), 
        x = interspecific$neighbours_number, 
        y = interspecific$seeds_uncorrected)

stan_dat$focal <- as.numeric(stan_dat$focal_labels)
stan_dat$comp <- as.numeric(stan_dat$comp_labels)
stan_dat$S <- max(stan_dat$focal)
stan_dat$N <- length(stan_dat$y)


library(rstan)
out <- stan(file = 'code/intraspecific_comp.stan', 
            data = stan_dat, 
            chains = 4, 
            cores = 4)

outr <- stan(file = 'code/bh_comp.stan', 
            data = stan_dat, 
            chains = 4, 
            cores = 4)

stan_alphas <- summary(outr, 'alpha' )$summary
stan_lambdas <- summary(outr, 'lambda')$summary

plot_diagnostics <- function(out, stan_dat ){ 
  
  lambda <- data.frame( summary(out, 'lambda')$summary)
  lambda$sp <- levels(stan_dat$species_labels)

  gg_l <- 
    lambda %>%
    ggplot( aes( x = sp, y = mean, ymin = X2.5., ymax = X97.5.)) + 
    geom_point() + 
    geom_errorbar() + 
    coord_flip()

  y_hat <- data.frame( summary(out, 'y_hat')$summary )
  y_hat <- data.frame( y_hat, obs = stan_dat$y)
  mu <- data.frame( summary(out, 'mu')$summary )
  mu <- data.frame( mu, obs = stan_dat$y)

  gg_y_hat <- 
    y_hat %>% 
    ggplot( aes( y = obs, x = X50., xmin = X2.5., xmax = X97.5.)) + 
    geom_errorbarh(color = 'blue') + 
    geom_point() + 
    geom_abline(aes( intercept = 0, slope = 1))
  
  gg_mu <- 
    mu %>% 
    ggplot( aes( y = obs, x = mean )) + 
    geom_point() + 
    scale_y_log10() + 
    scale_x_log10() + 
    geom_abline(aes( intercept = 0, slope = 1))
  
  return(list(gg_l, gg_y_hat, gg_mu))
}

plot_diagnostics(out, stan_dat)
plot_diagnostics(outr, stan_dat)



gg_dens <- 
  function(out, stan_dat){ 
  
  y_hat <- summary(out, 'y_hat')$summary
  
  df <- data.frame(y_hat, obs = stan_dat$y, x = stan_dat$x,  sp  = stan_dat$species_labels )
  
  df %>% 
    ggplot( aes( x = x, y = obs )) + 
    geom_point() + 
    facet_wrap(~sp, scales = 'free') +
    geom_ribbon( aes( ymin = X25., ymax = X75.), color = NA, fill = 'red', alpha = 0.2) + 
    scale_y_log10()
}


gg_dens(outr, stan_dat)
gg_dens(out, stan_dat)

lambda <- data.frame( sp = levels( stan_dat$focal_labels), lambda = stan_lambdas[, 1])

lambda <- lambda %>% 
  left_join(sp_names, by = c('sp' = 'prior_code'))

lambda %>% 
  left_join(demo_pars, by = c('USDA_symbol' = 'focal')) %>% 
  ggplot( aes( x = lambda.x, y = lambda.y)) + 
  geom_point() + 
  geom_abline( aes( intercept = 0, slope = 1))






stan_alphas_mean <- matrix( stan_alphas[, 1], stan_dat$S, stan_dat$S)
stan_alphas_low <- matrix( stan_alphas[, 4], stan_dat$S, stan_dat$S)
stan_alphas_hi <- matrix( stan_alphas[, 8], stan_dat$S, stan_dat$S)

alphas <- data.frame( sp = levels( stan_dat$focal_labels), stan_alphas_mean )
names( alphas)[-1] <- levels(stan_dat$focal_labels)

demo_pars <- read.csv('data/demo_pars.csv')

intras <- diag ( as.matrix( demo_pars[ , -c(1:4)] ) )
demo_pars$intra <- intras

sp_names <- 
  sedgwick_plants %>% 
  ungroup() %>% 
  select( USDA_symbol, prior_code ) %>% 
  filter( !is.na(prior_code)) %>% 
  distinct()

alphas <- 
  alphas %>% 
  left_join(sp_names, by = c('sp' = 'prior_code')) 

alphas <- alphas[, -1]

USDA_symbols <- sp_names$USDA_symbol
names(USDA_symbols) <- sp_names$prior_code

names(alphas)[1:21] <- USDA_symbols[ names(alphas)[1:21] ] 

alphas <- 
  alphas %>% 
  gather( comp, alpha, -USDA_symbol) %>% 
  rename( 'focal' = USDA_symbol) 

demo_pars <- 
  demo_pars %>% 
  select(-intra) %>% 
  gather( comp, alpha, AMME:URLI5) %>% 
  left_join( alphas, by = c('focal', 'comp')) 

demo_pars %>% 
  filter( focal == comp) %>% 
  ggplot(aes( x = alpha.x, y = alpha.y)) + 
  geom_point()  + 
  geom_abline(aes(intercept = 0 , slope = 1))

demo_pars %>% 
  ggplot(aes( x = alpha.x, y = alpha.y)) + 
  geom_point()  + 
  geom_abline(aes(intercept = 0 , slope = 1)) + 
  facet_wrap(~focal)




dat <- 
  dat %>% 
  left_join(sp_names, by = c('target' = 'prior_code' )) %>% 
  mutate( target = USDA_symbol ) %>% 
  left_join( sp_names, by = c('background' = 'prior_code')) %>% 
  mutate( background = USDA_symbol.y) %>% 
  select( -USDA_symbol.x, -USDA_symbol.y)

focal_traits <- 
  sedgwicktraits %>% 
  left_join(sedgwick_plants, by = c('species' = 'calflora_binomial')) %>% 
  filter( dataset == 'TAPIOCA') %>% 
  ungroup() %>%
  filter( !is.na(phenology))

demo_pars %>% 
  left_join(focal_traits ,  by = c('focal' = 'USDA_symbol')) %>% 
  ggplot( aes( x = rooting_depth, y = intra)) + 
  geom_point()

demo_pars %>% 
  left_join(focal_traits ,  by = c('focal' = 'USDA_symbol')) %>% 
  ggplot( aes( x = projected_area, y = intra)) + 
  geom_vline(aes(xintercept = 153))+ 
  scale_x_log10( ) + 
  geom_smooth(se = F) + 
  geom_point()

dat %>% 
  group_by( plot_number, target) %>% 
  filter( sampled == 1) %>% 
  summarise( n = n() ) %>% 
  ungroup() %>% 
  summarise( max(n), mean(n), min(n))

test <- 
  dat %>% 
  filter( background == target) %>% 
  group_by( plot_number ) %>% 
  filter( sampled == 1) %>% 
  summarise( n= n())

demo_pars %>% 
  left_join(focal_traits ,  by = c('focal' = 'USDA_symbol')) %>% 
  ggplot( aes( x = phenology, y = intra)) + 
  geom_point()

dat <- 
  dat %>% 
  left_join(focal_traits, by = c('target' = 'USDA_symbol')) %>% 
  left_join(focal_traits, by = c('background' = 'USDA_symbol'))



dat %>% 
  filter(!is.na(background) | target == background) %>% 
  ggplot( aes( x = density, y = sampled )) + 
  geom_point() + 
  geom_smooth(method = 'glm', formula = sampled ~ neighbours_number, method.args = list( family ='binomial')) + 
  facet_wrap( ~ target ) + 
  geom_point(data = sampled_per_d , color = 'red')

sampled_per_d <- 
  dat %>% 
  group_by( target, neighbours_number ) %>% 
  summarise( sampled = mean(sampled)) 



growth <- 
  dat %>% 
  filter( sampled == 1 ) %>% 
  mutate( seed_weight = seeds_uncorrected * seed_mass.x) %>% 
  filter( !is.na(seed_weight)) 

plot( growth$seed_weight, growth$weight)

growth_phenology <- 
  growth %>% 
  rename( 'final' = seed_weight, 
          'start' = seed_mass.x) %>% 
  gather( time, mass, c(start, final)) %>% 
  mutate( day = phenology.x) %>% 
  mutate( day = ifelse(time == 'start', 0, day) ) %>% 
  distinct( plot_number, grid_position, target, background, neighbours_number, day, mass )

growth_phenology %>% 
  ggplot( aes( x = day, y = log(mass), color = log(1 + neighbours_number))) + 
  stat_summary(geom = 'point', fun.y = 'mean', aes( group = paste0(target, neighbours_number)), alpha = 0.5) + 
  stat_summary(geom = 'line', fun.y = 'mean', aes( group = paste0(target, neighbours_number)), alpha = 0.5) + 
  facet_wrap(~target) + 
  scale_color_gradient(low = 'blue', high = 'red')
  
growth_rate <- 
  growth %>% 
  mutate( r = (log(seed_weight) - log(seed_mass.x))/phenology.x) %>% 
  mutate( background = ifelse(is.na(background), target, background)) 

growth_rate %>% 
  ggplot(aes( x = neighbours_number, y = r)) + geom_point() + 
  facet_wrap(~target)

intra_models <- 
  growth_rate %>% 
  filter( background == target ) %>% 
  group_by( target ) %>% 
  mutate( x1 = 1/(1 + neighbours_number), x2 = - log(1 + neighbours_number)) %>% 
  do( m1 = lm(data = .,  r ~ x1 ), m2 = lm( data = . , r ~ x2))

division_model <- intra_models$m1  
names(division_model) <- intra_models$target
division_model

log_model <- intra_models$m2
names(log_model) <- intra_models$target

intra_effects1 <- do.call( rbind, lapply( division_model, coef) )
intra_effects2 <- do.call(rbind, lapply( log_model, coef))

slopes <- data.frame(target = row.names(intra_effects1),  intra_effects1, intra_effects2 )

plot(slopes$x1, slopes$x2)

plot(slopes$x1, slopes$X.Intercept.)
plot(slopes$x2, slopes$X.Intercept..1)

intercepts <- 
  slopes %>%
  select( target, X.Intercept., X.Intercept..1) %>% 
  gather( type, intrcpt, c(X.Intercept., X.Intercept..1))

slopes <- 
  slopes %>% 
  select( target, x1, x2 ) %>% 
  gather( type, slope, c(x1, x2)) 

slopes %>% 
  ggplot(aes( x = target, color = type,  y = slope)) + 
  geom_point() + 
  coord_flip()

sp_names <- 
  sedgwick_plants %>% 
  ungroup() %>% 
  select( calflora_binomial, USDA_symbol) %>% 
  distinct()

traits <- 
  sedgwicktraits %>% 
  filter( dataset == 'TAPIOCA') %>% 
  left_join(sp_names, by = c('species' = 'calflora_binomial')) 

slopes %>% 
  left_join(traits, by = c('target' = 'USDA_symbol')) %>% 
  ggplot( aes( x = phenology, y = slope, color = type)) + 
  geom_point() + 
  geom_smooth(method = 'lm', se = F)

slopes %>% 
  left_join(traits, by = c('target' = 'USDA_symbol')) %>% 
  ggplot( aes( x = rooting_depth, y = slope, color = type)) + 
  geom_point() + 
  geom_smooth(method = 'lm', se = F)

intra_slope_df <- 
  slopes %>% 
  left_join(traits, by = c('target' = 'USDA_symbol')) 
  
m1 <- lm( data = intra_slope_df[intra_slope_df$type == 'x1', ], slope ~ rooting_depth + phenology + max_height + relative_spread + seed_mass + SRL)
summary(m1)



growth %>%
  mutate( rel_growth = log(seed_weight) - log(seed_mass.x )) %>% 
  mutate( rel_growth_rate = rel_growth/phenology.x) %>% 
  mutate( type = ifelse( target == background, 'intra', 'inter')) %>% 
  ggplot( aes( x = neighbours_number, y= rel_growth_rate )) + 
  geom_point() + 
  geom_smooth(method = 'lm', aes(group = type, color = type)) + 
  facet_wrap(~target)

dat %>% 
  ggplot( aes( x = neighbours_number, y = seeds)) + 
  geom_point() + 
  scale_y_log10()+ 
  facet_wrap(~target)


dat %>%
  filter( neighbours_number == 0) %>%
  ggplot( aes( x = seeds ) ) + 
  geom_histogram() + 
  facet_wrap(~target, scales = 'free')


sedgwicktraits %>% 
  filter( dataset == 'TAPIOCA') %>% 
  ggplot(aes( x = species, y = phenology)) + 
  geom_point() + 
  coord_flip()


sedgwicktraits %>% 
  filter( dataset == 'TAPIOCA') %>% 
  ggplot(aes( x = species, y = seed_mass)) + 
  geom_point() + 
  coord_flip()

PLER_DAT <- 
  dat %>% 
  filter( target == 'PLER', sampled == 1)





y <- PLER_DAT$seeds[PLER_DAT$neighbours_number == 0] 
library(fitdistrplus)

yrange <- 0:max( y)

y_int <- round(y)
y_norm <- fitdist(y, distr = 'norm')
y_lnorm <- fitdist(y, distr = 'lnorm')
y_pois <- fitdist(y_int, distr = 'pois')
y_gamm <- fitdist(y, distr = 'gamma')

dist_name <- 'norm'

fit_dist <- function( y, dist_name = 'norm', my_col = 'red'){ 
  
  yrange <- 0:max(y)
  
  fit <- fitdist(y, distr = dist_name) 
  
  my_dist <- eval(parse( text = paste0('d',dist_name)) )
  
  center <- fit$estimate[1]  
  if(length(fit$estimate) == 2){ 
    scale <- fit$estimate[2]
    qs <- my_dist(yrange, center, scale)
  }else if(length(fit$estimate) == 1){ 
    
    qs <- my_dist(yrange, center)  
  }

  points(yrange, qs, type = 'l', col = my_col)
}

hist(y, freq = F)
fit_dist(y, 'norm')
fit_dist(y_int, 'pois', 'blue')
fit_dist(y, 'gamma', 'orange')
fit_dist(y, 'lnorm', 'green')

y <- dat$seeds[  dat$target == 'HECO' & dat$sampled == 1]
y_int <- round(y)
hist(y, freq = F)
fit_dist(y, 'norm')
fit_dist(y_int, 'pois', 'blue')
fit_dist(y, 'gamma', 'orange')
fit_dist(y, 'lnorm', 'green')

fitdist(y, 'gamma')


