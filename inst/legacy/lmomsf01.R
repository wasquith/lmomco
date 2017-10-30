"lmomsf01" <-
function(x, f=NULL, a=0, nmom=5, subdivisionscale=50,
         sort=TRUE, constantails=TRUE, efs=1e-3) {

  if(nmom < 1) {
    warning("Number of L-moments is less than 1")
    return()
  }
  if(length(unique(x)) == 1) {
    warning("all values are equal--Lmoments can not be computed")
    return()
  }

  if(sort) {
    ix <- sort(x, index.return=TRUE)$ix
     x <-  x[ix]
    if(! is.null(f)) f <- f[ix]
  }

  if(is.null(f)) {
    f <- pp(x, a=a)
  }
  n <- length(x)

  if(nmom > n) {
    warning("More Lmoments requested by parameter 'nmom' than data points available in 'x'")
    return()
  }
  if(length(x) != length(f)) {
    warning("Number of 'x' values is not equal to number of 'f' values")
    return()
  }
  if(n < 5) {
    warning("Too few 'x' values, need > 5")
    return()
  }

  # Nonexceedances
  f1p0 <- f[1]   # left edge
  f1p1 <- f[2]   # one in from left edge
  fnm1 <- f[n-1] # one in from right edge
  fnm0 <- f[n]   # right edge

  # Values
  x1p0 <- x[1]   # left edge
  x1p1 <- x[2]   # one in from left edge
  xnm1 <- x[n-1] # one in from right edge
  xnm0 <- x[n]   # right edge

  # slopes right and left expressed as a positive
  delfr <- fnm0 - fnm1
  delfl <- f1p1 - f1p0
  sr <- ifelse(delfr > efs, (xnm0 - xnm1) / delfr, 0)
  sl <- ifelse(delfl > efs, (x1p1 - x1p0) / delfl, 0)

  "XofF.func" <- function(F) {
     if(constantails) return(approx(f, y=x, xout=F, rule=2)$y)
     X <- vector(mode="numeric", length=length(F))
       tf  <- F > fnm0 # right tail
     X[tf] <- xnm0 + sr * (F[tf] - fnm0)
       tf  <- F < f1p0 # left tail
     X[tf] <- x1p0 - sl * (f1p0 - F[tf])
       tf  <- (F <= fnm0 & F >= f1p0) # middle
     X[tf] <- approx(f, y=x, xout=F[tf], rule=2)$y
     return(X)
  }

  "PstarR.func" <- function(r,F) {
     Pr <- vector(mode="numeric", length=length(F))
     for(i in 1:length(F)) {
       tmp <- sapply(0:r,
         function(k) {
             return((-1)^(r-k)*choose(r,k)*choose(r+k,k)*F[i]^k)
          } )
       Pr[i] <- sum(tmp)
     }
     return(Pr)
  }

  "lam.func" <- function(F,r) {
     xf <- XofF.func(F)
     Ps <- PstarR.func(r-1,F)
     return(xf*Ps)
  }

  L <- vector(mode="numeric",length=nmom)
  R <- vector(mode="numeric",length=nmom)
  for(r in 1:nmom) {
     my.moment <- NA
     try(my.moment <- integrate(lam.func, 0, 1, r=r,
                       subdivisions=n*subdivisionscale)$value)
     L[r] <- my.moment
  }
  if(nmom >= 2) {
    R[2] <- L[2]/L[1]
  }
  if(nmom >= 3) {
    for(r in 3:nmom) {
      R[r] <- L[r]/L[2]
    }
  }
  R[1] <- NA

  z <- list(lambdas = L, ratios = R,
            trim=NULL, leftrim=NULL, rightrim=NULL,
            source="lmomsf")
  return(z)
}
