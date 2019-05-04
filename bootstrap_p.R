
boot_p_inc <- function (Spec) {
  nT <- Spec[1]
  Spec <- Spec[-1]
  Sobs <- sum(Spec > 0)
  Q1 <- sum(Spec == 1)
  Q2 <- sum(Spec == 2)
  Q0.hat <- ifelse(Q2 == 0, (nT - 1)/nT * Q1 * (Q1 - 1)/2, 
                   (nT - 1)/nT * Q1^2/2/Q2)
  A <- ifelse(Q1 > 0, nT * Q0.hat/(nT * Q0.hat + Q1), 1)
  a <- Q1/nT * A
  b <- sum(Spec/nT * (1 - Spec/nT)^nT)
  if (Q0.hat == 0) {
    w <- 0
    if (sum(Spec > 0) == 1) {
      warning("This site has only one species. Estimation is not robust.")
    }
  }
  else {
    w <- a/b
  }
  Prob.hat <- Spec/nT * (1 - w * (1 - Spec/nT)^nT)
  Prob.hat.Unse <- rep(a/ceiling(Q0.hat), ceiling(Q0.hat))
  return(c(Prob.hat, Prob.hat.Unse))
}


###########################################
#' Estimating species relative abundance 
#' Chao, A., Wang, Y. T., and Jost, L. (2013). 
#' Entropy and the species accumulation curve: a novel estimator of entropy via discovery rates of new species.
#' Appendix S2:  Estimating species relative abundance :one parameter
#' boot_p_abu(x) is a function of estimating detected species relative abundance.
#' @param x a vector of species abundance frequency
#' @return a numerical vector
boot_p_abu <- function(Spec)
{
  Sobs <- sum(Spec > 0)   #observed species
  n <- sum(Spec)        #sample size
  f1 <- sum(Spec == 1)   #singleton 
  f2 <- sum(Spec == 2)   #doubleton
  f0.hat <- ifelse(f2 == 0, (n - 1) / n * f1 * (f1 - 1) / 2, (n - 1) / n * f1 ^ 2/ 2 / f2)  #estimation of unseen species via Chao1
  A <- ifelse(f1>0, n*f0.hat/(n*f0.hat+f1), 1)
  a <- f1/n*A
  b <- sum(Spec / n * (1 - Spec / n) ^ n)
  if(f0.hat==0){
    w <- 0
    if(sum(Spec>0)==1){
      warning("This site has only one species. Estimation is not robust.")
    }
  }else{
    w <- a / b      	#adjusted factor for rare species in the sample
  }
  Prob.hat <- Spec / n * (1 - w * (1 - Spec / n) ^ n)					#estimation of relative abundance of observed species in the sample
  Prob.hat.Unse <- rep(a/ceiling(f0.hat), ceiling(f0.hat))  	#estimation of relative abundance of unseen species in the sample
  return(c(Prob.hat, Prob.hat.Unse))		  							#Output: a vector of estimated relative abundance
}




