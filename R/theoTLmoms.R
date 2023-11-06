"theoTLmoms" <-
function(para, nmom=5, trim=NULL, leftrim=NULL, rightrim=NULL,
               minF=0, maxF=1, quafunc=NULL,
               nsim=50000, fold=5,
               silent=TRUE, verbose=FALSE, ...) {
  if(nmom < 1) {
    warning("Number of TL-moments requested is less than 1")
    return()
  }
  if(! is.null(trim) && trim < 0) {
    warning("Trimming value is less than 0")
    return()
  }
  if(! is.null(leftrim) && leftrim < 0) {
    warning("Left rimming value is less than 0")
    return()
  }
  if(! is.null(rightrim) && rightrim < 0) {
    warning("Right trimming value is less than 0")
    return()
  }

  if(is.null(trim) && is.null(leftrim) && is.null(rightrim)) {
    trim <- 0
  }

  # t1 is the left-hand trimming and t2 is the right-hand trimming
  t1  <- NULL
  t2  <- NULL
  if(length(trim) == 1 && trim >= 0) {
    t1 <- t2 <- trim
    leftrim  <- rightrim <- NULL
  }
  else {
    trim <- NULL
    if(length(leftrim)  == 1 && leftrim  >= 0) t1 <- leftrim
    if(length(rightrim) == 1 && rightrim >= 0) t2 <- rightrim
    if(is.null(leftrim) ) { leftrim <- 0; t1 <- 0 }
    if(is.null(rightrim)) { rightrim <- 0; t2 <- 0 }
  }

  if(is.null(t1) || is.null(t2)) {
    warning("Ambiguous asymmetrical trimming values--use explicit leftrim ",
            "and rightrim arguments")
    return(NULL)
  }

  VR <- NULL # data frame built from verbose=TRUE on the numerical integration
  MC <- rep(NA, length=nmom)
  zz <- list(lambdas=rep(NA, nmom), ratios=rep(NA, nmom), trim=trim,
             leftrim=leftrim, rightrim=rightrim, monte_carlo=MC, source="theoTLmoms")

  L <- vector(mode="numeric", length=nmom)
  R <- vector(mode="numeric", length=nmom)
  for(r in seq_len(nmom)) { # for each  order of moment
    sum <- 0
    for(k in seq(0, r-1, by=1) ) {
      tmp <- (-1)^k*choose(r-1,k)
      tmp <- tmp * exp( lgamma(r+t1+t2+1) - lgamma(r+t1-k-1+1) - lgamma(t2+k+1) )
      # Quantile function X(F), which will require numerical integration
      XofF <- NULL
      if(is.null(quafunc)) {
         XofF <- function(FF) par2qua(FF, para, paracheck=FALSE) * FF^(r+t1-k-1) * (1-FF)^(t2+k)
      } else {
         XofF <- function(FF) quafunc(FF, para)                  * FF^(r+t1-k-1) * (1-FF)^(t2+k)
      }
      # Perform the numerical integration
      int <- NULL; used_monte <- FALSE
      try( int <- integrate(XofF, minF, maxF), silent=silent )
      # Error in integrate(XofF, minF, maxF) : the integral is probably divergent
      if(is.null(int)) { # Perform Monte Carlo integration
        used_monte <- TRUE
        folds <- rep(NA, fold )
        for(f in seq_len(fold )) folds[f] <- mean(XofF(runif(nsim)))
        muc <- mean(     folds               )
        mad <- mean( abs(folds - mean(folds)))
        int <- list(lmr.order=r, k.order=k, value=muc, abs.error=mad,
                    subdivisions="--", message="--", monte_carlo=used_monte)
        if(! silent) {
          print(str(int))
        }
      } else {
        int$lmr.order   <- r
        int$k.order     <- k
        int$monte_carlo <- used_monte
      }
      MC[r] <- used_monte
      sum <- sum + tmp * int$value # Sum up
      if(verbose) { # Handy messages
        ni <- c(int$lmr.order, int$k.order, round(int$value, digits=8), round(int$abs.error, digits=8), int$subdivisions, int$message, int$monte_carlo)
        names(ni) <- c("lmr.order", "k.order", "value", "abs.error", "subdivisions", "message", "monte_carlo")
        VR <- rbind(VR, ni) # if(verbose) print(ni)
      }
    }
    L[r] <- sum / r  # do not forget to divide by the order of the L-moment!
  }

  R[1] <- NA # because the first L-moment ratio is always and forevermore NA
  if(nmom >= 2)                       R[2] <- L[2] / L[1]
  if(nmom >= 3) for(r in seq(3,nmom)) R[r] <- L[r] / L[2]

  if(! any(MC == TRUE)) {
    nsim <- folds <- "not needed"
  }

  zz <- list(lambdas=L, ratios=R, trim=trim,
             leftrim=leftrim, rightrim=rightrim,
             nsim=nsim, folds=folds,
             monte_carlo=MC, source="theoTLmoms")

  VR <- as.data.frame(VR)
  row.names(VR) <- NULL
  if(verbose) zz$integrations <- VR

  return(zz)
}
