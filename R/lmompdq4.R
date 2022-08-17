"lmompdq4" <- function(para, paracheck=TRUE) {
  zz <- list(lambdas=rep(NA, 5), ratios=rep(NA, 5),
            trim=0, leftrim=0, rightrim=0, ifail=0, ifailtext="",
            source="lmompdq4")
  if(paracheck == TRUE) {
    if(! are.parpdq4.valid(para)) return()
  }
  U <- para$para[1]
  A <- para$para[2]
  K <- para$para[3]
  if(K > 0.99) {
    zz$ifail <- 1
    zz$ifailtext <- paste0("kappa > 0.99, later alpha results could be unreliable, ",
                           "if alpha back computed by lmompdq4(parpdq4())")
  }
  zz$lambdas[1] <- U
  zz$lambdas[c(3,5)] <- zz$ratios[c(3,5)] <- 0

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
  zz$lambdas[2] <- L2
  zz$ratios[4]  <- T4
  zz$ratios[2]  <- L2 / U
  zz$lambdas[4] <- L2 * T4
  return(zz)
}
