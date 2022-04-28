"pdfgpa" <-
function(x, para, paracheck=TRUE) {
  if(paracheck) {
    if(! are.pargpa.valid(para)) return()
  }
  U <- para$para[1]
  A <- para$para[2]
  K <- para$para[3]
  Y <- (x - U) / A
  ZERO <- sqrt(.Machine$double.eps)
  if(abs(K) > ZERO) {
    suppressWarnings(Y <- -(1/K) * log(1 - K*Y)) # traps NaN
  }
  f <- 1/A * exp(-(1-K)*Y)

  # conditions below will silently fix the NaN
  if(K > 0) {
    f[x < U      ] <- 0
    f[x > U + A/K] <- 0
  } else {
    f[x < U      ] <- 0
  }
  names(f) <- NULL
  return(f)
}

# Disabled version with bad support for X, reported by
# Christophe Dutang during comparisons of GPA/GPD implementations
# in R and various packages in Spring 2022.
#"pdfgpa" <-
#function(x,para) {
#    if(! are.pargpa.valid(para)) return()
#    XI <- para$para[1]
#    A  <- para$para[2]
#    K  <- para$para[3]
#
#    Y <- (x - XI)/A
#    if(K == 0) {
#        f <- A^(-1) * exp(-Y)
#    } else {
#        ARG <- 1-K*Y
#        Y <- suppressWarnings( -log(ARG)/K )
#        f <- A^(-1) * exp(-(1-K)*Y)
#    }
#
#    names(f) <- NULL
#    f[! is.finite(f)] <- NA
#    f[is.na(f)] <- 0 # decision Dec. 2015
#    return(f)
#}


