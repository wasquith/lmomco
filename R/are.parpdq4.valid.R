"are.parpdq4.valid" <- function(para, nowarn=FALSE) {
  if(! is.pdq4(para)) return(FALSE)
  if(any(is.na(para$para))) return(FALSE)

  op <- options()
  if(nowarn == TRUE) options(warn=-1)
  GO <- TRUE

  ALPHA <- para$para[2]
  KAPPA <- para$para[3]
  if(ALPHA <= 0) {
    warning("Parameter ALPHA is not > 0, invalid")
    GO <- FALSE
  }
  if(KAPPA > 1) {
    warning("Parameter KAPPA is not < 1, invalid")
    GO <- FALSE
  }
  options(op)
  if(GO) return(TRUE)
  return(FALSE)
}
