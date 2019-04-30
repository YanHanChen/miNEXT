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