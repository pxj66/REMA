setwd("./Real Data Analysis/Automaticity")
data0 <- read.csv("action_planning_PH_dataset.csv")
data1 <- data.frame(pa = data0$PA,
                    auto = (data0$automaticity_1 + data0$automaticity_2 + 
                              data0$automaticity_3 + data0$automaticity_4)/4,
                    plan = (data0$planning_frequency + data0$planning_manner + 
                              data0$planning_place+ data0$planning_time)/4,
                    bmi = data0$BMI,
                    ill = data0$Injury_illness,
                    sex = data0$Sex)

ls.m<-lm(auto ~. - pa, data = data1)
ls.y<-lm(pa~., data = data1)

library(quantreg)
lad.m<-rq(auto ~. - pa, data = data1)
lad.y<-rq(pa ~., data = data1)

source("./Simulation/Ck.R")
CkR <- Ck(ls.m)
k1 = CkR[1]
tau1 = CkR[2]
CkR = Ck(ls.y)
k2 = CkR[1]
tau2 = CkR[2]

n = nrow(data1)
X <- as.matrix(subset(data1, select = - c(auto, pa)))
X <- cbind(matrix(rep(1, n), nrow = n), X)
Z <- as.matrix(subset(data1, select = - pa))
Z <- cbind(matrix(rep(1, n), nrow = n), Z)

V <- solve(t(X) %*% X)
U <- solve(t(Z) %*% Z)

coeff1 <- coef(ls.m)
res1 <- resid(ls.m)
for (s in 1:20) {
  oldcoeff1 <- coeff1
  E1 <- res1
  W1 <- psi(E1, k1) / E1
  wls1 <- lm(auto ~. - pa, weights = W1, data = data1)
  coeff1 <- coef(wls1)
  res1 <- resid(wls1)
  if(sum(abs(coeff1 - oldcoeff1))<1e-5) 
    break
}
c(coeff1[2], sqrt(tau1 * V[2,2]))

coeff2 <- coef(ls.y)
res2 <- resid(ls.y)
for (s in 1:20) {
  E2 = res2
  oldcoeff2 <- coeff2
  W2 = psi(E2, k2) / E2
  wls2 <- lm(pa ~ ., weights = W2, data = data1)
  coeff2 <- coef(wls2)
  res2 <- resid(wls2)
  if(sum(abs(oldcoeff2 - coeff2))<1e-5) break
}
c(coeff2[2], sqrt(tau2 * U[2,2]))
coeff1["plan"] * coeff2["auto"]

################################################################################
coef(ls.m)["plan"] * coef(ls.y)["auto"]
summary(ls.m)$coefficient[2,2]
summary(ls.y)$coefficient[2,2]

coef(lad.m)["plan"] * coef(lad.y)["auto"]
summary.rq(lad.m, se = "iid")$coefficient[2,2]
summary.rq(lad.y, se = "iid")$coefficient[2,2]

coeff1["plan"] * coeff2["auto"]
c(coeff1[2],sqrt(tau1 * V[2,2]))
c(coeff2[3],sqrt(tau2 * U[2,2]))

##M1 Sobel test
library(RMediation)
medci(mu.x = coef(ls.m)["plan"], mu.y = coef(ls.y)["auto"], 
      se.x = summary(ls.m)$coefficient[2,2],
      se.y = summary(ls.y)$coefficient[2,2],
      type="all")$`Asymptotic Normal`
medci(mu.x = coef(lad.m)["plan"], mu.y = coef(lad.y)["auto"], 
      se.x = summary.rq(lad.m,se="iid")$coefficient[2,2],
      se.y = summary.rq(lad.y,se="iid")$coefficient[2,2],
      type = "all")$`Asymptotic Normal`
medci(mu.x = coeff1[2], mu.y = coeff2[3], se.x = sqrt(tau1*V[2,2]),
      se.y = sqrt(tau2*U[2,2]), type = "all")$`Asymptotic Normal`

##M2 PRCT and BCa
library(boot)
ls.fit<-function(d,i){
  ls.m<-lm(auto ~. - pa, data = d[i,])
  ls.y<-lm(pa ~., data = d[i,])
  coef(ls.m)["plan"]*coef(ls.y)["auto"]
}

lad.fit <- function(d, i){
  lad.m <- rq(auto ~. - pa, data = d[i, ])
  lad.y <- rq(pa ~ ., data = d[i, ])
  coef(lad.m)["plan"]*coef(lad.y)["auto"]
}

mr.fit <- function(d, i){
  ls.m<-lm(auto ~. - pa, data = d[i,])
  ls.y<-lm(pa ~. , data = d[i,])
  CkR <- Ck(ls.m)
  k1 <- CkR[1]
  CkR <- Ck(ls.y)
  k2 <- CkR[1]
  
  coeff1 <- coef(ls.m)
  res1 <- resid(ls.m)
  for (s in 1:5) {
    oldcoeff1 <- coeff1
    E1 <- res1
    W1 <- psi(E1, k1) / E1
    wls1 <- lm(auto ~. - pa, weights = W1, data = d[i, ])
    coeff1 <- coef(wls1)
    res1 <- resid(wls1)
    if(sum(abs(coeff1 - oldcoeff1))<1e-5) 
      break
  }
  
  coeff2 <- coef(ls.y)
  res2 <- resid(ls.y)
  for (s in 1:5) {
    E2 = res2
    oldcoeff2 <- coeff2
    W2 = psi(E2, k2) / E2
    wls2 <- lm(pa ~ ., weights = W2, data = d[i, ])
    coeff2 <- coef(wls2)
    res2 <- resid(wls2)
    if(sum(abs(oldcoeff2 - coeff2))<1e-5) break
  }
  coeff1["plan"]*coeff2["auto"]
}

ls.boot<-boot(data1, ls.fit, R = 2000)
boot.ci(ls.boot)

lad.boot<-boot(data1, lad.fit, R = 2000)
boot.ci(lad.boot)

mr.boot<-boot(data1, mr.fit, R = 2000)
boot.ci(mr.boot)










