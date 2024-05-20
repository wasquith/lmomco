"pargdd" <-
function(lmom, checklmom=TRUE, symgdd=FALSE, init.para=NULL,
               silent=FALSE, trace=FALSE, control=list(abstol=0.001, maxit=1000), ...) {

  zz <- list(type="gdd", para=c(NA,NA,NA,NA, as.numeric(symgdd)), source="pargdd")
  if(length(lmom$L1) == 1) { # convert to named L-moments
    lmom <- lmorph(lmom)     # nondestructive conversion!
  }
  if(checklmom & ! are.lmom.valid(lmom)) {
    warning("L-moments are invalid")
    return()
  }

  L1 <- lmom$lambdas[1]
  L2 <- lmom$lambdas[2]
  L3 <- lmom$lambdas[3]
  L4 <- lmom$lambdas[4]

  small <- 1E-4
  ofunc <- function(par) {
    para <- exp(par)
    para <- list(para=para, type="gdd")
    if(symgdd) {
      para$para[3] <- para$para[1]
      para$para[4] <- para$para[2]
      para$para[5] <- 1
    } else {
      para$para[5] <- NA
    }
    if(! are.pargdd.valid(para)) return(Inf); # message(para$para)
    tlmr <- lmomgdd(para, nmom=4, paracheck=FALSE)
    mu <- (para$para[1]/para$para[2]) - (para$para[3]/para$para[4])
    if( abs(tlmr$lambdas[1] - mu) > small) return(Inf)
    err <- (tlmr$lambdas[1] - L1)^2 + (tlmr$lambdas[2] - L2)^2 +
           (tlmr$lambdas[3] - L3)^2 + (tlmr$lambdas[4] - L4)^2
    if(trace) {
      #print(lmom$lambdas)
      #print(tlmr$lambdas)
      message("TRACE ofunc : para ", paste(c(round(para$para, digits=6), err), collapse=", ", sep=""))
    }
    return(sqrt(err))
  }

  integrate(function(f, n=1) { qgamma(f, A1, rate=B1)^n - qgamma(f, A2, rate=B2)^n }, 0, 1, n=1)$value ???
  # E[X^n] = sum(sapply(0:n) function(k) choose(n,k)*(-1)^k*E[X1^(n-k)]*E[X2^k] )
  mk <- function(n) {
    txt <- ""
    for(k in 0:n) {
      a <- choose(n,k)
      b <- (-1)^k
      b <- ifelse(b == -1, "-", "+")
      c <- paste0("E[X1^", n-k, "]")
      d <- paste0("E[X2^",   k, "]")
      txt <- paste0(txt, " ", paste0(b, a,c,d))
    }
    print(txt)
  }
  pmfunc <- function(par, pm=NA) {
    para <- exp(par); A1 <- para[1]; B1 <- para[2]; A2 <- para[3]; B2 <- para[4]
    mugdd <-      A1/B1 - A2/B2  # arithmetic mean
    sdgdd <- sqrt(A1/B1 + A2/B2) # standard deviation
    skgdd <- 2*(A1*B2^3 - A2*B1^3) / (A1*B2^2 + A2*B1^2)^(3/2) # revise to 3rd moment, not skew
    #ktgdd <- ???? # need to derive
    err <- (mugdd - pm[1])^2 + (sdgdd - pm[2])^2 + (skgdd - pm[3])^2 + (skgdd - pm[3])^2
    return(err)
  }

  if(is.null(init.para)) {
    #pm <- pmoms(init.para, checklmom=FALSE) # when testing starts, init.para sneak way to get x into function
    #rt <- NULL
    #try(rt <- optim(c(0,0,0,0), pmfunc, pm=pm))
    #print(rt)
    #init.para <- exp(rt$par)
    #stop()
    init.para <- c(10, 3, 2, 0.2)
    if(symgdd) init.para <- init.para[c(1:2, 1:2)]
  } else {
    if(length(init.para) == 5) {
      if(! is.na(init.para[5]) & init.para[5] == 1) symgdd <- TRUE
      init.para <- init.para[1:4]
    }
    if(symgdd) init.para <- init.para[c(1:2, 1:2)]
  }
  init.para <- log(init.para)
  if(trace) message("TRACE : log init.para = ", paste(init.para, collapse=", ", sep=""))
  rt <- NULL
  try(rt <- optim(init.para, ofunc, control=control, ...), silent=silent)
  if(is.null(rt)) {
    zz$optim <- rt
    return(zz)
  }
  para <- exp(rt$par)
  if(symgdd) {
    para[3] <- para[1]
    para[4] <- para[2]
    para[5] <- 1
  }
  zz$para <- para
  zz$optim <- rt
  return(zz)
}

