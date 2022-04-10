"quapdq4" <- function(f, para, paracheck=TRUE) {
  if(! check.fs(f)) return()
  if(paracheck == TRUE) {
    if(! are.parpdq4.valid(para)) return()
  }
  U <- para$para[1]
  A <- para$para[2]
  K <- para$para[3]
  if(K >= 1-sqrt(.Machine$double.eps)) K <- 1-sqrt(.Machine$double.eps)
  LF  <- log(f/(1-f))
  KF2 <- K*(2*f - 1)
  if(abs(K) < sqrt(.Machine$double.eps)) { # logistic
    x <- U + A*LF # https://en.wikipedia.org/wiki/Logistic_distribution
  } else if(K > 0) {
    x <- U + A*(LF - 2*K*atanh(KF2))
  } else {
    x <- U + A*(LF + 2*K*atan( KF2))
  }
  names(x) <- NULL
  return(x)
}
