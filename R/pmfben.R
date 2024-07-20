"pmfben" <- function(d, para=list(para=c(1, 10)), ...) {
  m <- para$para[1]
  b <- para$para[2]
  if(b != 10) {
    warning("pmfben is only configured for base10 at the moment")
    return(NULL)
  }
  d <- floor(d)
  s <- (1*b^(m-1))
  e <- as.integer(paste(rep(9, m), collapse=""))
  f <- logb(1 + 1/d, base=b)
  f[d < s] <- 0
  f[d > e] <- 0
  f
}
