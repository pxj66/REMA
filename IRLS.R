##########################################################################
####package
setwd("./Simulation")
source("Ck.R")


####function
psi <- function(x, k){
  sapply(x, function(t){
    y <- numeric(0)
    if(abs(t) <= k)
      y <- t
    else
      y <- k * sign(t)
    return(y)
  })
}

####parameter and data
beta2 = 0
beta3 = 0
cp = 1
a = 0.39
b = 0.39

n <- 200
x <- rnorm(n)

e2 <- rt(n, 2)
e3 <- rt(n, 2)
m <- beta2 + a * x + e2
y <- beta3 + cp * x + b * m + e3
X <- matrix(c(rep(1, n), x), nrow = n)
Z <- matrix(c(rep(1, n), x, m), nrow = n)

V <- solve(t(X) %*% X)
U <- solve(t(Z) %*% Z)
##########################################################################
ls1 <- lm(m ~ x)
ls2 <- lm(y ~ x + m)
CkR <- Ck(ls1)
k1 <- CkR[1]
tau1 <- CkR[2]
CkR <- Ck(ls2)
k2 <- CkR[1]
tau2 <- CkR[2]

coeff1 <- coef(ls1)
for (s in 1:100) {
  oldcoeff1 <- coeff1
  E1 <- m - coeff1[1] - coeff1[2] * x
  W1 <- psi(E1, k1) / E1
  wls1 <- lm(m ~ x, weights = W1)
  coeff1 <- coef(wls1)
  if(sum(abs(coeff1 - oldcoeff1)) < 1e-5) 
    break
}
wls1
summary(wls1)
se.a <- sqrt(tau1*V[2,2])
se.b <- sqrt(tau2*U[3,3])


