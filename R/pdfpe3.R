"pdfpe3" <-
function(x,para) {
  if(! are.parpe3.valid(para)) return()

  names(para$para) <- NULL
  MU    <- para$para[1] # location
  SIGMA <- para$para[2] # scale
  GAMMA <- para$para[3] # shape

  SMALL <- sqrt(.Machine$double.eps)
  if(abs(GAMMA) <= SMALL) return(dnorm(x, mean=MU, sd=SIGMA))
  ALPHA <- 4/GAMMA^2
  BETA  <-   (1/2) * SIGMA * abs(GAMMA)
  XI    <- MU - 2  * SIGMA/GAMMA
  Y <- sign(GAMMA) * (x - XI)
  f <- dgamma(Y/BETA, ALPHA)/BETA
  names(f) <- NULL
  f[! is.finite(f)] <- NA
  f[is.na(f)] <- 0 # decision Dec. 2015
  return(f)
}
