"lmompdq4" <- function(para, paracheck=TRUE) {
  z <- list(lambdas=rep(NA, 5), ratios=rep(NA, 5),
            trim=0, leftrim=0, rightrim=0,
            source="lmompdq4")
  if(paracheck == TRUE) {
    if(! are.parpdq4.valid(para)) return()
  }
  U <- para$para[1]
  A <- para$para[2]
  K <- para$para[3]
  if(K > 0.98) {
    warning("kappa > 0.98, later alpha results could be unreliable, ",
            "if alpha back computed by lmompdq4()")
  }
  z$lambdas[1] <- U
  z$lambdas[c(3,5)] <- z$ratios[c(3,5)] <- 0

  if(K > 0) {
    L2 <- A*(1-K^2)*atanh(K)/K
    T4 <- -(1/4) + (5/(4*K)) * (1/K - 1/atanh(K))
  } else if(abs(K) < sqrt(.Machine$double.eps)) {
    L2 <- A
    T4 <- 1/6
  } else {
    L2 <- A*(1+K^2)*atan(K)/K
    T4 <- -(1/4) - (5/(4*K)) * (1/K - 1/atan(K))
  }
  z$lambdas[2] <- L2
  z$ratios[4]  <- T4
  z$ratios[2]  <- L2 / U
  z$lambdas[4] <- L2 * T4
  return(z)
}
