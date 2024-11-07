######################################################################
## Huber loss, psi2 and psi.dot function
psi <- function(x,k){
  sapply(x, function(t){
    y <- numeric(0)
    if(abs(t) <= k)
      y <- t
    else
      y <- k * sign(t)
    return(y)
  })
}

######################################################################
## Choose the optimal tuning constant
Ck <- function(lad1){
  resid1 <- resid(lad1)
  n <- length(resid1)
  K <- round(3 * median(abs(resid1)) / 0.6745, 2)
  S <- seq(0.2, K, 0.01)
  i <- 1
  tau <- numeric(0)
  for (k in S) {
    B <- mean(abs(resid1) <= k)
    Sig2 <- (sum((resid1[abs(resid1) <= k])^2) + sum(abs(resid1) > k) * k^2)/ n
    tau[i] <- Sig2 / B^2
    i <- i + 1
  }
  return(c(S[which.min(tau)], tau[which.min(tau)]))
}

