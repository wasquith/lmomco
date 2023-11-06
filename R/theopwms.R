"theopwms" <-
function(para, nmom=5, minF=0, maxF=1, quafunc=NULL,
               nsim=50000, fold=5,
               silent=TRUE, verbose=FALSE, ...) {

  if(nmom < 1) {
    warning("Number of PWMs requested is less than 1")
    return()
  }

  VR <- NULL # data frame built from verbose=TRUE on the numerical integration
  MC <- rep(NA, length=nmom)

  B <- vector(mode="numeric", length=nmom)
  for(r in seq(0, nmom-1)) { # for each  order of moment
    XofF <- NULL
    if(is.null(quafunc)) {
      XofF <- function(FF) par2qua(FF, para, paracheck=FALSE) * FF^r
    } else {
      XofF <- function(FF) quafunc(FF, para) * FF^r
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
      int <- list(pwm.order=r, value=muc, abs.error=mad,
                  subdivisions="--", message="--", monte_carlo=used_monte)
      if(! silent) print(str(int))
    } else {
      int$pwm.order   <- r
      int$monte_carlo <- used_monte
    }
    MC[r+1] <- used_monte
    B[r+1] <- int$value
    if(verbose) { # Handy messages
      ni <- c(int$pwm.order, round(int$value, digits=8), round(int$abs.error, digits=8), int$subdivisions, int$message, int$monte_carlo)
      names(ni) <- c("pwm.order", "value", "abs.error", "subdivisions", "message", "monte_carlo")
      VR <- rbind(VR, ni) # if(verbose) print(ni)
    }
  }

  if(! any(MC == TRUE)) {
    nsim <- folds <- "not needed"
  }

  zz <- list(betas=B, nsim=nsim, folds=folds,
             monte_carlo=MC, source="theopwms")

  VR <- as.data.frame(VR)
  row.names(VR) <- NULL
  if(verbose) zz$integrations <- VR

  return(zz)
}
