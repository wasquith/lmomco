"lmompdq3" <- function(para, paracheck=TRUE) {
  z <- list(lambdas=rep(NA, 4), ratios=rep(NA, 4),
            trim=0, leftrim=0, rightrim=0,
            source="lmompdq3")
  if(paracheck == TRUE) {
    if(! are.parpdq3.valid(para)) return()
  }
  U <- para$para[1]
  A <- para$para[2]
  K <- para$para[3]

  if(K == -1) K <- -1 + .Machine$double.eps
  if(K == +1) K <- +1 - .Machine$double.eps

  TAU3 <- 1/K - 1/atanh(K)
  if(is.nan(TAU3)) TAU3 <- 0
  TAU4 <- ((5*TAU3 / K) - 1) / 4
  if(is.nan(TAU4)) TAU4 <- 1/6
  LAM2 <- A*(1-K^2) / (1-K*TAU3)
  LAM1 <- U + A*((1+K)*log(1+K) - (1-K)*log(1-K) - K*log(4))
  z$lambdas[1] <- LAM1
  z$lambdas[2] <- LAM2
  z$ratios[3]  <- TAU3
  z$ratios[4]  <- TAU4
  z$ratios[2]  <- LAM2 / LAM1
  z$lambdas[3] <- LAM2 * TAU3
  z$lambdas[4] <- LAM2 * TAU4
  return(z)
}
