#https://www.math.kit.edu/stoch/~klar/seite/veroeffentlichungen/media/note-vg-revision.pdf
"cdfgdd" <- function(x, para, paracheck=TRUE, silent=TRUE, ...) {
  if(paracheck) if(! are.pargdd.valid(para)) return()

  A1 <- para$para[1];  B1 <- para$para[2]
  A2 <- para$para[3];  B2 <- para$para[4]
  if(length(para$para) == 5) {
    if(! is.na(para$para[5]) & para$para[5] == 1) {
      A2 <- A1
      B2 <- B1
    }
  }

  d  <- B2^A2 / exp(lgamma(A1) + lgamma(A2))
  hi <- Inf

  pgdd <- function(t, x=NA) {
    t^(A2-1) * exp(-B2*t + lgamma(A1) + pgamma(B1*(t+x), A1, log.p=TRUE, lower.tail=TRUE))
  }
  cdf <- lapply(x, function(x) {
    if(x == -Inf) return(list(value=NA))
    if(x == +Inf) return(list(value=NA))
    int <- NULL
    try(int <- integrate(pgdd, pmax(0, -x), hi, x=x, ...), silent=silent)
    if(is.null(int)) return(list(value=NA))
    return(int)
  })
  n <- length(cdf)
  f <- d * sapply(seq_len(n), function(i) cdf[[i]]$value)
  f[x == -Inf] <- 0
  f[x == +Inf] <- 1
  return(f)
}
