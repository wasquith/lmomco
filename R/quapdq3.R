"quapdq3" <- function(f, para, paracheck=TRUE) {
  if(! check.fs(f)) return()
  if(paracheck == TRUE) {
    if(! are.parpdq3.valid(para)) return()
  }
  U <- para$para[1]
  A <- para$para[2]
  K <- para$para[3]

  f[f <     .Machine$double.eps] <-     .Machine$double.eps
  f[f > 1 - .Machine$double.eps] <- 1 - .Machine$double.eps

  if(K >=       1-sqrt(.Machine$double.eps))  K <-       1-sqrt(.Machine$double.eps)
  if(K <= -1 * (1-sqrt(.Machine$double.eps))) K <- -1 * (1-sqrt(.Machine$double.eps))
  LF <- log(f/(1-f))
  KK <- K * log( (1 - K*(2*f-1))^2 / (4*f*(1-f)))
  x <- U + A*(LF + KK)
  names(x) <- NULL
  return(x)
}
