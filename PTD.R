setwd("./Mental Health Violent Offending")
Base <- load("./ICPSR_29961/DS0001/29961-0001-Data.rda")
Base <- eval(parse(text = Base))

## Base$S0SRO6   # total offending variety proportion - past 6 months

## Base$S0BSISOM # mental health: somatization
## Base$S0BSIOC  # mental health: obsessive-compulsive
## Base$S0BSIIS  # mental health: interpersonal sensitivity
## Base$S0BSIDEP # mental health: depression
## Base$S0BSIANX # mental health: Anxiety
## Base$S0BSIHOS # mental health: hostility
## Base$S0BSIPHB # mental health: phobic anxiety
## Base$S0BSIPAR # mental health: paranoid ideation
## Base$S0BSIPSY # mental health: Psychoticism
## Base$S0BSIGSI # mental health: global severity index

## Base$S0EXPVIC # experienced victimization

## Base$S0GANG1  # 6mo before JV - member of a gang
Base$S0PAKNOW[is.na(Base$S0PAKNOW)] <- 0 # parent knowledge
Base$S0PARMNT[is.na(Base$S0PARMNT)] <- 0 # Parent monitoring
Base$S0MAWARM[is.na(Base$S0MAWARM)] <- 0 # parent warmth: mother 
Base$S0MAHOTL[is.na(Base$S0MAHOTL)] <- 0 # parent hostility: mother
Base$S0PAWARM[is.na(Base$S0PAWARM)] <- 0 # parent warmth: father
Base$S0PAHOTL[is.na(Base$S0PAHOTL)] <- 0 # parent hostility: father
## Base$S0ROUT   # Unsupervised routine activities
## Base$S0SUBUS2 # Subuse past 6 months how often used alcohol
## Base$S0SUBUS7 # subuse past 6 months num times used  marijuana

## Base$S0AGE    # DEM age
## Base$S0ETHN_R # DEM ethnicity, factor
## Base$S0SGEND  # DEM gender, factor
Famstr <- as.numeric(Base$S0FAMSTR)
Famstr[Famstr %in% c(7, 9, 10, 11)] <- 99 # non biological parents
Famstr[Famstr %in% c(2, 3, 4, 5, 6, 8, 12, 13)] <- 0 # single biological parents
Famstr[Famstr == 1] <- -99 # two biological parents
Famstr <- factor(Famstr, levels = c( -99, 0, 99), 
                 labels = c("two bio pat", "sing bio pat", "non bio pat"))
# DEM family structure, factor

MEV <- data.frame(offend = Base$S0SRO6, health = Base$S0BSIGSI, 
                  expvic = Base$S0EXPVIC, gang = as.numeric(Base$S0GANG1) - 1, 
                  patmnt = Base$S0PARMNT + Base$S0PAKNOW, 
                  patwarm = Base$S0MAWARM + Base$S0PAWARM,
                  pathotl = Base$S0PAHOTL + Base$S0MAHOTL, 
                  rout = Base$S0ROUT, alcohol = as.numeric(Base$S0SUBUS2) - 1, 
                  marijuana = as.numeric(Base$S0SUBUS7), age = Base$S0AGE, 
                  ethn = Base$S0ETHN_R, 
                  gend = as.numeric(Base$S0SGEND) - 1, 
                  famstr = Famstr)
MEV <- MEV[rowSums(is.na(MEV)) <= 2, ]
MEV <- MEV[!is.na(MEV$offend), ]
MEV <- MEV[!is.na(MEV$health), ]
MEV <- MEV[!is.na(MEV$gang), ]
MEV <- MEV[!is.na(MEV$alcohol), ]
MEV <- MEV[!is.na(MEV$marijuana), ]
MEV <- MEV[!MEV$patwarm == 0, ]
MEV <- MEV[!MEV$pathotl == 0, ]
MEV <- MEV[!MEV$patmnt == 0, ]

str(MEV) # check the number of NA
## colSums(is.na(MEV))
## table(rowSums(is.na(MEV)))

#install.packages("fastDummies")
library(fastDummies)
MEV1 <- dummy_cols(MEV, select_columns = c("ethn", "famstr"), 
                   remove_first_dummy = T,
                   remove_selected_columns = T)
str(MEV1)

ls1 <- lm(expvic~.-offend, data = MEV1)
ls2 <- lm(offend~., data = MEV1)

library(quantreg)
lad1 <- rq(expvic~.-offend, data = MEV1)
lad2 <- rq(offend~., data = MEV1)

## Huber estimate
source("./Simulation/Ck.R")
CkR <- Ck(lad1)
k1 = CkR[1]
tau1 = CkR[2]
CkR = Ck(lad2)
k2 = CkR[1]
tau2 = CkR[2]

n = nrow(MEV1)
X <- as.matrix(subset(MEV1, select = - c(expvic, offend)))
X <- cbind(matrix(rep(1, n), nrow = n), X)
Z <- as.matrix(subset(MEV1, select = - offend))
Z <- cbind(matrix(rep(1, n), nrow = n), Z)

V <- solve(t(X) %*% X)
U <- solve(t(Z) %*% Z)

coeff1 <- coef(lad1)
res1 <- resid(lad1)
for (s in 1:100) {
  oldcoeff1 <- coeff1
  E1 <- res1
  W1 <- psi(E1, k1) / E1
  wls1 <- lm(expvic ~. - offend, weights = W1, data = MEV1)
  coeff1 <- coef(wls1)
  res1 <- resid(wls1)
  if(sum(abs(coeff1 - oldcoeff1))<1e-5) 
    break
}
c(coeff1[2], sqrt(tau1 * V[2,2]))

coeff2 <- coef(lad2)
res2 <- resid(lad2)
for (s in 1:100) {
  E2 = res2
  oldcoeff2 <- coeff2
  W2 = psi(E2, k2) / E2
  wls2 <- lm(offend ~ ., weights = W2, data = MEV1)
  coeff2 <- coef(wls2)
  res2 <- resid(wls2)
  if(sum(abs(oldcoeff2 - coeff2))<1e-5) break
}
c(coeff2[3], sqrt(tau2 * U[3,3]))
coeff1[2] * coeff2[3]

## Inference
### M1: Sobel test
a.ls <- summary(ls1)$coefficients[2, 1]
a.se.ls <- summary(ls1)$coefficients[2, 2]
b.ls <- summary(ls2)$coefficients[3, 1]
b.se.ls <- summary(ls2)$coefficients[3, 2]
ls.se <- sqrt(a.ls^2 * b.se.ls^2 + b.ls^2 * a.se.ls^2)

a.lad <- summary.rq(lad1, se = "iid")$coefficients[2, 1]
a.se.lad <- summary.rq(lad1, se = "iid")$coefficients[2, 2]
b.lad <- summary.rq(lad2, se = "iid")$coefficients[3, 1]
b.se.lad <- summary.rq(lad2, se = "iid")$coefficients[3, 2]
lad.se <- sqrt(a.lad^2 * b.se.lad^2 + b.lad^2 * a.se.lad^2)

a.h <- coeff1[2]
a.se.h <- sqrt(tau1 * V[2,2])
b.h <- coeff2[3]
b.se.h <- sqrt(tau2 * U[3,3])
h.se <- sqrt(a.h^2 * b.se.h^2 + b.h^2 * a.se.h^2)

stest.ls.1 <- a.ls * b.ls - qnorm(0.975) * ls.se
stest.ls.2 <- a.ls * b.ls + qnorm(0.975) * ls.se
len.ls.stes <- stest.ls.2 - stest.ls.1
a.ls * b.ls
c(stest.ls.1, stest.ls.2)
len.ls.stes

stest.lad.1 <- a.lad * b.lad - qnorm(0.975) * lad.se
stest.lad.2 <- a.lad * b.lad + qnorm(0.975) * lad.se
len.lad.stes <- stest.lad.2 - stest.lad.1
a.lad * b.lad
c(stest.lad.1, stest.lad.2)
len.lad.stes

stest.h.1 <- a.h * b.h - qnorm(0.975) * h.se
stest.h.2 <- a.h * b.h + qnorm(0.975) * h.se
len.h.stes <- stest.h.2 - stest.h.1
a.h * b.h
c(stest.h.1, stest.h.2)
len.h.stes


### M2: PRCT and BCa
library(boot)
ls.boot <- function(d, i){
  ls1 <- lm(expvic~.-offend, data = d[i, ])
  ls2 <- lm(offend~., data = d[i, ])
  coef(ls1)["health"]*coef(ls2)["expvic"]
}

lad.boot <- function(d, i){
  lad1 <- rq(expvic~.-offend, data = d[i, ])
  lad2 <- rq(offend~., data = d[i, ])
  coef(lad1)["health"]*coef(lad2)["expvic"]
}

h.boot <- function(d, i){
  # lad1 <- rq(expvic~.-offend, data = d[i, ])
  # lad2 <- rq(offend~., data = d[i, ])
  ls1 <- lm(expvic~.-offend, data = d[i, ])
  ls2 <- lm(offend~., data = d[i, ])
  CkR <- Ck(lad1)
  k1 = CkR[1]
  tau1 = CkR[2]
  CkR = Ck(lad2)
  k2 = CkR[1]
  tau2 = CkR[2]
  
  coeff1 <- coef(ls1)
  res1 <- resid(ls1)
  for (s in 1:20) {
    oldcoeff1 <- coeff1
    E1 <- res1
    W1 <- psi(E1, k1) / E1
    wls1 <- lm(expvic ~. - offend, weights = W1, data = d[i, ])
    coeff1 <- coef(wls1)
    res1 <- resid(wls1)
    if(sum(abs(coeff1 - oldcoeff1))<1e-5) break
  }
  
  coeff2 <- coef(ls2)
  res2 <- resid(ls2)
  for (s in 1:20) {
    E2 = res2
    oldcoeff2 <- coeff2
    W2 = psi(E2, k2) / E2
    wls2 <- lm(offend ~ ., weights = W2, data = d[i, ])
    coeff2 <- coef(wls2)
    res2 <- resid(wls2)
    if(sum(abs(oldcoeff2 - coeff2))<1e-5) break
  }
  coeff1[2] * coeff2[3]
}

ls.out <- boot(MEV1, ls.boot, R = 2000, stype = "i")
boot.ci(ls.out, conf = 0.95, type = "all")

lad.out <- boot(MEV1, lad.boot, R = 2000, stype = "i")
boot.ci(lad.out, conf = 0.95, type = "all")

h.out <- boot(MEV1, h.boot, R = 2000, stype = "i")
boot.ci(h.out, conf = 0.95, type = "all")

#install.packages("moments")
library(moments)
skewness(resid(ls1))
kurtosis(resid(ls1))
qqnorm(resid(ls1), main = "y - m & x")
qqline(resid(ls1), col = "red")

ks.test(resid(ls1), "pnorm", mean = mean(resid(ls1)), sd = sd(resid(ls1)))

qqnorm(resid(ls2), main = "m - x")
qqline(resid(ls2), col = "red")
skewness(resid(ls2))
kurtosis(resid(ls2))
ks.test(resid(ls2), "pnorm", mean = mean(resid(ls2)), sd = sd(resid(ls2)))

par(mfrow = c(1,2))

















