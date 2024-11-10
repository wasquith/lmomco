"cdfben" <- function(d, para=list(para=c(1, 10)), ...) {
  m <- para$para[1]
  if(m > 9) {
    warning("expect failure as first parameter is an integer greater than 9")
  }
  b <- para$para[2]
  if(b != 10) {
    warning("pmfben is only configured for base10 at the moment in the second parameter")
    return(NULL)
  }
  d <- floor(d)
  s <- (1*b^(m-1))
  e <- as.integer(paste(rep(9, m), collapse=""))
  f <- sapply(d, function(t) { sum(logb(1 + 1/(s:t), base=b)) })
  f[d < s] <- 0
  f[d > e] <- 1
  return(f)
}
