"lmomgdd" <-
function(para, nmom=6, paracheck=TRUE, silent=TRUE, ...) {
  if(paracheck) if(! are.pargdd.valid(para)) return()

  A1 <- para$para[1];  B1 <- para$para[2]
  A2 <- para$para[3];  B2 <- para$para[4]
  if(length(para$para) == 5) {
    if(! is.na(para$para[5]) & para$para[5] == 1) {
      A2 <- A1
      B2 <- B1
    }
  }

  lower <- -Inf
  upper <- +Inf
  #lower <- quagdd(0.0001, para, paracheck=FALSE)
  #upper <- quagdd(1-0.0001, para, paracheck=FALSE)

  afunc <- function(x, r=NA, j=NA) {
     FFx <- cdfgdd(x, para, paracheck=FALSE)
     FFx^(r-j-1) * (1-FFx)^(j+1)
  }
  LAM <- vector(mode="numeric", nmom)
  LAM[1] <- A1/B1 - A2/B2
  for(r in 2:nmom) {
    Lr <- sapply(0:(r-2), function(j) {
             a <- (-1)^j * choose(r-2, j) * choose(r, j+1)
             int <- NULL
             try(int <- integrate(afunc, lower=lower, upper=upper, r=r, j=j, ...), silent=silent)
             if(is.null(int)) return(NA)
             return(a * int$value) })
    LAM[r] <- (1/r) * sum(Lr)
  }
  RAT <- rep(NA, nmom)
  if(nmom > 1) {
    RAT[2] <- LAM[2] / LAM[1]
    if(nmom > 2) {
      for(r in 3:nmom) RAT[r] <- LAM[r] / LAM[2]
    }
  }
  zz <- list(lambdas=LAM, ratios=RAT, trim=0, leftrim=0, rightrim=0, source="lmomgdd")
  return(zz)
}
