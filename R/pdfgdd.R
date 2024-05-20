"pdfgdd" <- # https://arxiv.org/pdf/2309.09516v2
function(x, para, paracheck=TRUE, silent=TRUE, ...) {
  if(paracheck) if(! are.pargdd.valid(para)) return()

  A1 <- para$para[1];  B1 <- para$para[2]
  A2 <- para$para[3];  B2 <- para$para[4]
  if(length(para$para) == 5) {
    if(! is.na(para$para[5]) & para$para[5] == 1) {
      A2 <- A1
      B2 <- B1
    }
  }

  f <- rep(NA, length(x))

  c <- (B1^A1 * B2^A2) / exp( lgamma(A1) + lgamma(A2) )

  A1m1 <- A1 - 1; A2m1 <- A2 - 1; BB   <- -(B1 + B2)

  wnt <- x > 0
  d <- rep(NA, length(x))
  d[  wnt] <- exp( B2*x[  wnt])
  d[! wnt] <- exp(-B1*x[! wnt])

  f <- sapply(seq_len(length(x)), function(i) {
    int <- NULL
    if(x[i] > 0) {
      try(int <- integrate(function(z) z^A1m1*(z-x[i])^A2m1*exp(BB*z), +x[i], Inf, ...), silent=silent)
    } else {
      try(int <- integrate(function(z) z^A2m1*(z+x[i])^A1m1*exp(BB*z), -x[i], Inf, ...), silent=silent)
    }
    if(is.null(int)) return(NA)
    return(int$value) })

  f <- c * d * f
  f[f < 0] <- 0
  return(f)
}


