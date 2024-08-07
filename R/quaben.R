"quaben" <- function(f, para=list(para=c(1, 10)), ...) {
  m <- para$para[1]
  b <- para$para[2]
  if(b != 10) {
    warning("pmfben is only configured for base10 at the moment")
    return(NULL)
  }
  s <- (1*b^(m-1))
  e <- as.integer(paste(rep(9, m), collapse=""))
  x <- cut(f, breaks=c(0, sapply(s:e, function(k) sum(pmfben(s:k, para=para)))),
                          labels=s:e, include.lowest=TRUE)
  as.integer( as.character(x) )
}
