"lmoms" <-
function(x, nmom=5, no.stop=FALSE, vecit=FALSE) {
  n <- length(x)

  if(nmom > n) {
    if(no.stop) return(NULL)
    stop("More L-moments requested by parameter 'nmom' than data points available in 'x'")
  }

  if(length(unique(x)) == 1) {
    if(no.stop) return(NULL)
    stop("all values are equal--Lmoments can not be computed")
  }

  z <- TLmoms(x,nmom=nmom)
  z$source <- "lmoms"
  if(! vecit) return(z)
  if(nmom == 1) {
     z <- z$lambdas[1]
  } else if(nmom == 2) {
     z <- c(z$lambdas[1], z$lambdas[2])
  } else {
     z <- c(z$lambdas[1], z$lambdas[2], z$ratios[3:nmom])
  }
  attr(z, which="trim")     <- NULL
  attr(z, which="rightrim") <- NULL
  attr(z, which="leftrim")  <- NULL
  attr(z, which="source")   <- "lmoms"
  return(z)
}
