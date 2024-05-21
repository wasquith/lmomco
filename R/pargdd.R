"pargdd" <-
function(lmom, checklmom=TRUE, symgdd=FALSE, init.para=NULL, snap.tau4=FALSE,
               silent=FALSE, trace=FALSE, control=list(abstol=0.0001, maxit=1000), ...) {

  zz <- list(type="gdd", para=c(NA,NA,NA,NA, as.numeric(symgdd)), source="pargdd")
  if(length(lmom$L1) == 1) { # convert to named L-moments
    lmom <- lmorph(lmom)     # nondestructive conversion!
  }
  if(checklmom & ! are.lmom.valid(lmom)) {
    warning("L-moments are invalid")
    return()
  }

  message <- ""

  nudge  <- 0.001
  pe3t3t4 <- lmrdia()$pe3
  pe3t4   <- approx(pe3t3t4[,1], pe3t3t4[,2], xout=lmom$ratios[3], rule=2)$y
  if(snap.tau4 & lmom$ratios[4] <= pe3t4) {
    lmom$ratios[ 4] <- pe3t4 + nudge
    lmom$lambdas[4] <- lmom$lambdas[2] * lmom$ratios[4]
    message <- "T4 was below PE3, snapped it to a nudge above, "
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
    mu <- (para$para[1] / para$para[2]) - (para$para[3] / para$para[4])
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

  # E[X^n] = sum(sapply(0:n) function(k) choose(n,k)*(-1)^k*E[X1^(n-k)]*E[X2^k] )
  #mk <- function(n) {
  #  txt <- ""
  #  for(k in 0:n) {
  #    a <- choose(n,k)
  #    b <- (-1)^k
  #    b <- ifelse(b == -1, "-", "+")
  #    c <- paste0("E[X1^", n-k, "]")
  #    d <- paste0("E[X2^",   k, "]")
  #    txt <- paste0(txt, " ", paste0(b, a,c,d))
  #  }
  #  print(txt)
  #}

  if(is.null(init.para)) {
    F1 <- runif(5000); F2 <- runif(5000)
    mcfunc <- function(par, lmr=NA) {
      p <- exp(par); if(symgdd) p <- c(p[1:2], p[1:2])
      tlmr <- pwm2lmom( pwm( qgamma(F1, p[1], p[2]) -  qgamma(F2, p[3], p[4]) ) )
      if(! are.lmom.valid(tlmr)) return(Inf)
      pe3t4 <- approx(pe3t3t4[,1], pe3t3t4[,2], xout=tlmr$ratios[3], rule=2)$y
      if(tlmr$ratios[4] <= pe3t4) return(Inf) # below the Pearson III
      sum((lmr$lambdas[1:4] - tlmr$lambdas[1:4])^2)
    }
    try(sara <- optim(c(0, 0, 0, 0), mcfunc, lmr=lmom), silent=FALSE )
    if(is.null(sara)) {
      init.para <- c(10, 3, 2, 0.2)
    } else {
      init.para <- exp(sara$par)
      message <- paste0(message, "initial parameters from Monte Carlo qgamma()x2")
    }
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
  zz$para      <- para
  zz$message   <- message
  zz$optim     <- rt
  zz$init.para <- exp(init.para)
  return(zz)
}
