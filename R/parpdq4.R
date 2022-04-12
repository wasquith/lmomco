"parpdq4" <- function(lmom, checklmom=TRUE, snapt4uplimit=TRUE) {
  para <- rep(NA, 3)
  names(para) <- c("xi", "alpha", "kappa")
  if(length(lmom$L1) == 0) {
    lmom <- lmorph(lmom)
  }
  if(checklmom & ! are.lmom.valid(lmom)) {
    warning("L-moments are invalid")
    return()
  }

  # By having neginf not hideously deep in magnitude, we greatly reduce the iteration counts
  # for uniroot() when kappa < 0, so we have a truncation on the depth here but notice that
  # we are still in a few parts per million
  neginf  <- -.Machine$double.xmax^(1/64)
  # print(-(1/4) - (5/(4*neginf)) * (1/neginf - 1/atan(neginf)), 16)
  smallTAU4 <- -0.2499878576145593

  bigTAU4 <- 0.874 # see demo in parpdq4.Rd

  para[1] <- lmom$L1
  LAM2 <- lmom$L2
  TAU4 <- lmom$TAU4
  if(is.null(TAU4) || is.na(TAU4)) {
     warning("The fourth L-moment ratio is undefined")
     return()
  }
  if(TAU4 > bigTAU4) {
     warning("The tau4 is too big for the algorithm (tau4 <= ", bigTAU4, ")")
     if(! snapt4uplimit) return()
     warning("reducing tau4 to the upper margin (tau4 <- ", bigTAU4, ")")
     TAU4 <- bigTAU4
  }

  if(TAU4 < smallTAU4) TAU4 <- smallTAU4 + sqrt(.Machine$double.eps)
  #print(TAU4, 16) # -0.2499878427133981

  if(1/6 <= TAU4 & TAU4 <= 1) {
    fn <- function(K) {
       val <- -(1/4) + (5/(4*K)) * (1/K - 1/atanh(K))
       if(is.nan(val)) val <- 1/6
       #points(K, val, pch=21, bg=4)
       return(val - TAU4)
    }
    rt <- NULL
    try(rt <- uniroot(fn, interval=c(0,1)))
    para[3] <- rt$root
    para[2] <- LAM2 * para[3] / ((1-para[3]^2)*atanh(para[3]))
    if(is.nan(para[2])) para[2] <- LAM2 # logistic limit
  } else if(TAU4 < 1/6) {
    fn <- function(K) {
       val <- -(1/4) - (5/(4*K)) * (1/K - 1/atan(K))
       if(is.nan(val)) val <- 1/6
       #print(c(K, TAU4, val))
       #points(K, val, pch=21, bg=2)
       return(val - TAU4)
    }
    rt <- NULL # could use -.Machine$double.xmax^0.25 or even further but, 1E5 is
    try(rt <- uniroot(fn, interval=c(neginf, 0)), silent=TRUE) # tau4=-0.249992
    if(is.null(rt)) {
      message("rooting for solution to Tau4 failed, ",
              "kappa too small, setting kappa to lower limit")
      para[3] <- neginf
    } else {
      para[3] <- rt$root#; print(rt$root)
    }
    para[2] <- LAM2 * para[3] / ((1+para[3]^2)*atan(para[3]))
    if(is.nan(para[2])) {
      print(para[3])
    }
  } else {
    warning("what to do about TAU4 == 1?") # well are.lmom.valid() would catch
    return(NULL)
  }
  if(para[3] > 0.98) {
    warning("kappa > 0.98, alpha (yes alpha) results could be unreliable")
  }
  zz <- list(para=para, type="pdq4", source="parpdq4")
  return(zz)
}
