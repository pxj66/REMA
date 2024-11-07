#work space
setwd("./Simulation")

#package and source
library(MASS)
library(quantreg)
source("Ck.R")

#huber loss weighted function
psi<-function(x, k){
  sapply(x, function(x){
    y = numeric()
    if (abs(x)>k)
      y = k*sign(x)
    else 
      y = x
    return(y)
  }) 
}

set.seed(1234)
N<-1000

#result matrix: ls, lad, mr
ls<-matrix(0,N,4)
colnames(ls)<-c('a','se.a','b','se.b')

lad<-matrix(0,N,4)
colnames(lad)<-colnames(ls)

mr<-matrix(0,N,4)
colnames(mr)<-colnames(ls)

#sample size and effect
beta2 = 0
beta3 = 0
cp = 1
a = 0
b = 0.39

n <- 200
x <- rnorm(n)

for (i in 1:N) {
  #data
  e2 <- c(rnorm(0.9*n), rnorm(0.1*n,0,100))
  e3 <- c(rnorm(0.9*n), rnorm(0.1*n,0,100))
  m <- beta2 + a * x + e2
  y <- beta3 + cp * x + b * m + e3
  X <- matrix(c(rep(1, n), x), nrow = n)
  Z <- matrix(c(rep(1, n), x, m), nrow = n)
  
  V <- solve(t(X) %*% X)
  U <- solve(t(Z) %*% Z)
  
  #ls
  ls1 <-lm(m ~ x)
  ls2 <- lm(y ~ x + m)
  ls[i,] <- c(summary(ls1)$coefficients[2, 1:2],summary(ls2)$coefficients[3, 1:2])
  
  #lad
  lad1 <- rq(m ~ x)
  lad2 <- rq(y~x+m)
  lad[i,] <- c(summary.rq(lad1,se="iid")$coefficients[2, 1:2],
             summary.rq(lad2,se="iid")$coefficients[3, 1:2])
  
  #mr
  #STEP 1 Choose optimal tuning constant
  CkR <- Ck(lad1)
  k1 <- CkR[1]
  tau1 <- CkR[2]
  CkR <- Ck(lad2)
  k2 <- CkR[1]
  tau2 <- CkR[2]
  
  #STEP 2 IRLS
  coeff1 <- coef(lad1)
  for (s in 1:100) {
    oldcoeff1 <- coeff1
    E1 <- m - coeff1[1] - coeff1[2] * x
    W1 <- psi(E1, k1) / E1
    wls1 <- lm(m ~ x, weights = W1)
    coeff1 <- coef(wls1)
    if(sum(abs(coeff1 - oldcoeff1)) < 1e-5) 
      break
  }
  mr[i,c("a","se.a")]<-c(coeff1[2], sqrt(tau1 * V[2,2]))
  
  coeff2 <- coef(lad2)
  for (s in 1:100) {
    E2 = y - coeff2[1] - coeff2[2] * x - coeff2[3] * m
    oldcoeff2 <- coeff2
    W2 = psi(E2, k2) / E2
    wls2 <- lm(y ~ x + m, weights = W2)
    coeff2 <- coef(wls2)
    if(sum(abs(oldcoeff2 - coeff2)) < 1e-5) break
  }
  mr[i, c("b", "se.b")] <- c(coeff2[3], sqrt(tau2 * U[3, 3]))
  cat('iteration = ', iter <- i, '\n')
}

##

save(ls, lad, mr, file = 'a=0,b=0.39_n=200_Err=Mixture.RData')


ls.se <- sqrt((ls[ , 'a'] * ls[ , 'se.b'])^2 + (ls[ , 'b'] * ls[ , 'se.a'])^2)
ls.ci <- matrix(0, N, 2)
ls.ci[ , 1] <- ls[ , 'a'] * ls[ , 'b'] - qnorm(0.975) * ls.se
ls.ci[ , 2] <- ls[ , 'a'] * ls[ , 'b'] + qnorm(0.975) * ls.se
mean(ls.ci[ , 1] > 0 | ls.ci[ , 2] < 0)

lad.se <- sqrt((lad[ , 'a'] * lad[ , 'se.b'])^2 + (lad[ , 'b'] * lad[ , 'se.a'])^2)
lad.ci <- matrix(0, N, 2)
lad.ci[ , 1] <- lad[ , 'a'] * lad[ , 'b'] - qnorm(0.975) * lad.se
lad.ci[,2] <- lad[ , 'a'] * lad[ , 'b'] + qnorm(0.975) * lad.se
mean(lad.ci[ , 1] > 0 | lad.ci[ , 2] < 0)

mr.se <- sqrt((mr[ , 'a'] * mr[ , 'se.b'])^2 + (mr[ , 'b'] * mr[ , 'se.a'])^2)
mr.ci <- matrix(0, N, 2)
mr.ci[ , 1] <- mr[ , 'a'] * mr[ , 'b'] - qnorm(0.975) * mr.se
mr.ci[,2] <- mr[ , 'a'] * mr[ , 'b'] + qnorm(0.975) * mr.se
mean(mr.ci[ , 1] > 0 | mr.ci[ , 2] < 0)

