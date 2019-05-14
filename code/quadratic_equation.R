
quadratic_equation <- function( y, a, b ) { 
  
  v1 <- (y + sqrt(4*a*b*y + y^2) )/(2*a)
  v2 <- (y - sqrt(4*a*b*y + y^2) )/(2*a) 
  
  c(v1,v2)
}


a <- 1 
y <- 10 
b <- 0.5 

test <- data.frame( do.call(cbind, my_dat) )
test <- test[!test$intra, ]

a <- test$d1/test$m
u2 <- u$u1[test$s2]
u1 <- u$u1[test$s1]

PHI <- phi$phi[ test$p]
D <- test$D
X <- test$x
b <- u2*PHI*D*X
y <- test$y

data.frame( y = y, a = a, b_part = PHI*D*X  )

u_solve <- matrix( NA, nrow = length(y), ncol = 2)

for( i in 1:length(y)){ 
  u_solve[i, ] <- quadratic_equation(y = y[i], a = a[i] , b = b[i])
}

matplot( u_solve )
points(u1, col = 'blue')

all(u_solve[,1]>0 )
all(u_solve[, 2] < 0)

