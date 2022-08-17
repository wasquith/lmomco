"parpdq3" <- function(lmom, checklmom=TRUE) {
  para <- rep(NA, 3)
  names(para) <- c("xi", "alpha", "kappa")
  if(length(lmom$L1) == 0) {
    lmom <- lmorph(lmom)
  }
  if(checklmom & ! are.lmom.valid(lmom)) {
    warning("L-moments are invalid")
    return()
  }

  LAM1 <- lmom$L1
  LAM2 <- lmom$L2
  TAU3 <- lmom$TAU3

  fn <- function(K) {
     val <- 1/K - 1/atanh(K)
     if(is.nan(val)) val <- 0
     return(val - TAU3)
  }
  rt <- NULL
  try(rt <- uniroot(fn, interval=c(-1, 1)))
  if(is.null(rt)) {
    para[3] <- K <- 0
  } else {
    para[3] <- K <- rt$root
  }
  para[2] <- LAM2 * (1 - K*TAU3) / (1 - K^2)
  para[2] <- exp(log(LAM2) + log(1 - K*TAU3) - log(1 - K^2))
  para[1] <- LAM1 - para[2]*((1+K)*log(1+K) - (1-K)*log(1-K) - K*log(4))

  zz <- list(para=para, type="pdq3", ifail=0, ifailtext="", source="parpdq3")

  if(abs(para[3]) > 0.98) {
    zz$ifail <- 1
    zz$ifailtext <- "|kappa| > 0.98, alpha (yes alpha) results could be unreliable"
  }
  return(zz)
}
