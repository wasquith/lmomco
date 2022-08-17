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

  ifail = ""
  # By having neginf not hideously deep in magnitude, we greatly reduce the iteration counts
  # for uniroot() when kappa < 0, so we have a truncation on the depth here but notice that
  # we are still in a few parts per million
  neginf  <- -.Machine$double.xmax^(1/64)
  # print(-(1/4) - (5/(4*neginf)) * (1/neginf - 1/atan(neginf)), 16)
  smallTAU4 <- -0.2499878576145593

  bigTAU4 <- 0.866 # see code described in parpdq4.Rd and lmompdq4.Rd
  bigTAU4 <- 0.845 # see code described in parpdq4.Rd and lmompdq4.Rd

  para[1] <- lmom$L1
  LAM2 <- lmom$L2
  TAU4 <- lmom$TAU4
  if(is.null(TAU4) || is.na(TAU4)) {
     warning("The fourth L-moment ratio (tau4) is undefined")
     return()
  }
  if(TAU4 > bigTAU4) {
     txt <- paste0("tau4 is deemed too large (tau4 <= ", bigTAU4, ")\n")
     if(snapt4uplimit) {
       warning(txt, "  reducing to upper algorithm margin (tau4 <- ", bigTAU4, ")")
       TAU4 <- bigTAU4
     } else {
       warning(txt, "  not reducing to upper algorithm margin (tau4 <- ", bigTAU4, ")")
     }
  }

  if(TAU4 < smallTAU4) TAU4 <- smallTAU4 + sqrt(.Machine$double.eps)
  #print(TAU4, 16) # -0.2499878427133981

  if(1/6 <= TAU4 & TAU4 < 1) {
    fn <- function(K) {
       val <- -(1/4) + (5/(4*K)) * (1/K - 1/atanh(K))
       if(is.nan(val)) val <- 1/6
       return(val - TAU4)
    }
    rt <- NULL
    try(rt <- uniroot(fn, interval=c(0, 1), tol=.Machine$double.eps^0.5))
    para[3] <- rt$root
    para[2] <- LAM2 * para[3] / ( (1 - para[3]^2) * atanh(para[3]) )
    if(is.nan(para[2])) para[2] <- LAM2 # logistic limit
  } else if(TAU4 < 1/6) {
    fn <- function(K) {
       val <- -(1/4) - (5/(4*K)) * (1/K - 1/atan(K))
       if(is.nan(val)) val <- 1/6
       return(val - TAU4)
    }
    rt <- NULL # could use -.Machine$double.xmax^0.25 or even further but, 1E5 is
    try(rt <- uniroot(fn, interval=c(neginf, 0), tol=.Machine$double.eps^0.5), silent=TRUE) # tau4=-0.249992
    if(is.null(rt)) {
      message("rooting for solution to Tau4 failed, ",
              "kappa too small, setting kappa to lower limit")
      para[3] <- neginf
    } else {
      para[3] <- rt$root#; print(rt$root)
    }
    para[2] <- LAM2 * para[3] / ( (1 + para[3]^2) * atan(para[3]) )
    if(is.nan(para[2])) {
      ifail <-
      return(list(para=para, type="pdq4", ifail=ifail,
                  ifailtext=paste0("alpha is NaN with a kappa=", para[2], "on L2=", LAM2),
                  source="parpdq4"))
    }
  } else {
    return(list(para=para, type="pdq4", ifail=1,
                ifailtext="TAU4 >= 1 is incompatible",
                source="parpdq4"))
  }
  if(para[3] > 0.99) {
    return(list(para=para, type="pdq4", ifail=ifail,
                ifailtext="kappa > 0.99, alpha (yes alpha) results could be unreliable",
                source="parpdq4"))
  }
  return(list(para=para, type="pdq4", ifail=0, ifailtext="", source="parpdq4"))
}
