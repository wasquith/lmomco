"TLmom" <-
function(x, order=1, trim=NULL, leftrim=NULL, rightrim=NULL, sortdata=TRUE) {
  r <- order
  n <- length(x)

  if(r > n) {
    stop("TLmoment 'order' greater than data points available in 'x'")
  }

  if(! is.null(trim) && (! is.null(leftrim) || ! is.null(rightrim))) {
    warning("Ambiguous trimming arguments")
    return()
  }

  if(length(unique(x)) == 1) stop("all values are equal--TLmoment can not be computed")

  t1 <- NULL
  t2 <- NULL
  if(is.null(trim) && is.null(leftrim) && is.null(rightrim)) {
    trim <- 0
  }

  if(length(trim) == 1 && trim >= 0) {
    t1 <- trim
    t2 <- trim
    leftrim <- NULL
    rightrim <- NULL
  }
  else {
    if(length(leftrim)  == 1 && leftrim  >= 0) t1 <- leftrim
    if(length(rightrim) == 1 && rightrim >= 0) t2 <- rightrim
    if(is.null(leftrim) ) {  leftrim <- 0; t1 <- 0 }
    if(is.null(rightrim)) { rightrim <- 0; t2 <- 0 }
  }

  if(is.null(t1) || is.null(t2)) {
    warning("Ambiguous asymmetrical trimming values--use explicit leftrim and rightrim arguments")
    return()
  }

  if(sortdata == TRUE) x <- sort(x)
  denom  <- lchoose(n, r+t1+t2)
  lambda <- sum(sapply(seq(t1+1, n-t2), function(i) {
                        wk <- sapply(seq(0, r-1), function(k) {
                 (-1)^k * exp(lchoose(r-1, k) + lchoose(i-1, r+t1-1-k) + lchoose(n-i, t2+k) - denom)
                         })
                   wk * x[i] }) )
  lambda <- lambda/r

  z <- list(lambda = lambda, order = r, trim = trim, leftrim = t1, rightrim = t2)
  return(z)
}

