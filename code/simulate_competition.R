rm(list = ls())
library(boot)
# simulate competition 

nsp <- 20 
nt1 <- 3  # number of traits influencing u
nt2 <- 3  # number of traits influencing phi

# trait covariates 
X1 <- matrix(NA, nsp, nt1)  # traits for uptake 
X1[] <- rnorm(nsp*nt1)

X1 <- cbind(1, X1) # add intercept

X2 <- matrix(NA, nsp, nt2)  # traits for overlap 
X2[] <- rnorm(nsp*nt1)

# phenology 
d <-  rnbinom( nsp , size = 60 , prob = 0.5) 

# demo pars 
m <- rlnorm(nsp, 1, 1)/1000  # seed mass (g)

g <- rbeta(nsp, 5, 1) # germination

s <- rbeta(nsp, 4, 1) # seed survival 

# Pairwise differences  

all_pairs <- expand.grid( 1:nsp, 1:nsp)
D <- apply( all_pairs, 1, function(x, d) min(1, d[x[1]]/d[x[2]]), d = d)

# Set parameters ----------------------- # 

# vector of coefficients that determines net resource uptake rate
# from the functional trait values in X1
Beta1 <- c(-4, 0, -0.3, 0.3)

# Resource gain rates determined by traits 
u <- exp( X1 %*% Beta1 ) # constrain to be positive 


# vector of coefficients that weight pairwise differences 
# between species. Positive values only. 
Beta2 <- c(0, 0.3, 2) 

# Calc pairwise niche overlap based on trait difference with 
# weights for each trait given by Beta2 
X2_scaled <- sweep(X2, 2, Beta2 , '*')
trait_dist <- dist(X2_scaled)
phi <- inv.logit( - log (dist(X2_scaled) ))

plot(trait_dist, phi)

all_pairs$intra <- NA
all_pairs$intra[all_pairs$Var1 == all_pairs$Var2] <- T
all_pairs$intra[is.na(all_pairs$intra)] <- F

pair <- apply( all_pairs[, 1:2], 1, function(x) sort(c(x)) )
inter_pairs <- pair[, !all_pairs$intra ]

all_pairs$phi <- NA
all_pairs$phi <- as.numeric( as.matrix(phi) )
all_pairs$phi[all_pairs$intra] <- 1

all_pairs$p <- apply( all_pairs[, 1:2], 1, function(x) paste( sort(x), collapse = '_'))

all_pairs$p <- factor( all_pairs$p, levels = unique( all_pairs$p ))

# Simulate response to competition 
sim_dat  <- 
  data.frame(s1 = all_pairs[,1],
             s2 = all_pairs[,2], 
            u1 = u[ all_pairs[, 1]], 
            u2 = u[ all_pairs[, 2]], 
            d1 = d[all_pairs[,1]], 
            D = D, 
            phi = all_pairs$phi, 
            m = m[all_pairs[,1]],
            p = all_pairs$p,
            intra = all_pairs$intra)

# Merge with experimental grid of densities 
exp_grid <- expand.grid( s1 = 1:nsp, s2 = 1:nsp, x = c(1, 2, 4, 8, 16, 32) - 1, block = 1:4 )

# Next step removes duplicate lambda plots 
exp_grid$s2[ exp_grid$x==0] <- exp_grid$s1[ exp_grid$x == 0 ] 
exp_grid <- unique( exp_grid) 
# 

exp_grid <- merge( exp_grid, sim_dat, by = c('s1', 's2'))


calc_mu <- function(x, u1, u2, phi, d1, D, m){ 
  
  (1/m)*(u1*d1)/(1 + (u2/u1)*phi*D*x )
}

calc_mu <- Vectorize(calc_mu)

exp_grid$mu <- calc_mu(exp_grid$x, exp_grid$u1, exp_grid$u2, exp_grid$phi, exp_grid$d1, exp_grid$D, exp_grid$m)

# simulate response with error from log normal #
exp_grid$y <- rpois(nrow(exp_grid), exp_grid$mu)

# -----------------------------------------  #
plot(exp_grid$x, exp_grid$y, type = 'n')

for( i in unique(exp_grid$s1)){ 
  y <- exp_grid[exp_grid$s1 == i & exp_grid$s2 != i, 'y']
  x <- exp_grid[exp_grid$s1 == i & exp_grid$s2 != i, 'x']
  
  points(x,  y , col = i)
}

# format for fitting the response via the traits
exp_grid <- exp_grid[ , c('mu', 'y', 'x', 's1','s2', 'u1', 'u2', 'p', 'phi', 'd1', 'D', 'm', 'block', 'intra') ]

X11 <- X1[ exp_grid$s1, ]
X12 <- X1[ exp_grid$s2, ]
X21 <- X2[exp_grid$s1, ]
X22 <- X2[exp_grid$s2, ]

u <- unique( exp_grid[, c('s1', 'u1')])
u <- u[order(u$s1), ]

inter_pairs <- unique( exp_grid[ !exp_grid$intra, c('s1', 's2', 'p') ] )
inter_pairs$p <- factor( inter_pairs$p )

exp_grid <- merge( exp_grid[, -grep('^p$', names(exp_grid))], inter_pairs, by = c('s1', 's2'), all.x = T)

phi <- unique(exp_grid[!exp_grid$intra, c('p', 'phi')]) 
phi <- phi[order(phi$p), ]

save(u, phi, Beta1, Beta2, exp_grid, X11, X12, X21, X22, file = 'output/simulated_data.rda' )


