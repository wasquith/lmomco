"x2pars" <- function(x, verbose=TRUE, ...) {
      if(verbose) message("Working on L-moments")
      lmr <- lmr2par(x, ...)
      if(verbose) message("Working on Maximum Likelihood")
      mle <- mle2par(x, para.int=lmr, ...)
      if(verbose) message("Working on Maximum Product of Spacings")
      mps <- mps2par(x, para.int=lmr, ...)
      zz <- list(lmr=lmr, mle=mle, mps=mps)
      return(zz)
}
