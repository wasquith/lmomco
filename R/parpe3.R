"parpe3" <-
function(lmom, checklmom=TRUE, useHosking=TRUE, ...) {
    para <- rep(NA, 3)
    names(para) <- c("mu", "sigma", "gamma")

    if(! useHosking) {
      # Coefficients from the embedded example in parpe3.Rd from which refitting to numerically
      # integrated Tau3 from theoLmoms() for a given gamma is used with curation of settings for
      # integrate() and optim() call to estimate other C and D parameters (also using a new SMALL).
      # However, we use Hosking original structure and use Hosking's parameters as initial starting
      # points in a 9D optim()ization. These coefficients result from the script
      # lmomco/inst/newcoes/makePE3coes.R.
      SMALL <- .Machine$double.eps
      C1 <-  0.292135202635314
      C2 <-  0.189791143390214
      C3 <-  0.0440498237092469
      D1 <-  0.36066214333264
      D2 <- -0.595506123861203
      D3 <-  0.252854052270765
      D4 <- -2.788853998752
      D5 <-  2.56107741806592
      D6 <- -0.770667511918909
    } else {
      # METHOD: RATIONAL APPROXIMATION IS USED TO EXPRESS ALPHA, THE SHAPE
      # PARAMETER OF THE GAMMA DISTRIBUTION, AS A FUNCTION OF TAU-3.
      # RELATIVE ACCURACY OF THE APPROXIMATION IS BETTER THAN 3E-5.

      # SMALL IS USED TO TEST WHETHER SKEWNESS IS EFFECTIVELY ZERO
      SMALL <- 1e-6

      # CONSTANTS USED IN MINIMAX APPROXIMATIONS
      C1 <-  0.2906
      C2 <-  0.1882
      C3 <-  0.0442
      D1 <-  0.36067
      D2 <- -0.59567
      D3 <-  0.25361
      D4 <- -2.78861
      D5 <-  2.56096
      D6 <- -0.77045
    }

    PI3    <-    3*pi
    ROOTPI <- sqrt(pi)

    if(length(lmom$L1) == 0) { # convert to named L-moments
      lmom <- lmorph(lmom)     # nondestructive conversion!
    }

    if(checklmom & ! are.lmom.valid(lmom)) {
      warning("L-moments are invalid")
      return()
    }

    L1 <- lmom$L1
    L2 <- lmom$L2
    T3 <- abs(lmom$TAU3)
    if(T3 <= SMALL) {
      # ZERO SKEWNESS
      para[1] <- L1
      para[2] <- L2*ROOTPI
      para[3] <- 0
      return(list(type="pe3", para=para, source="parpe3"))
    }
    if(T3 >= 1/3) {
      TT <- 1 - T3
      ALPHA <- TT * ( D1 + TT * ( D2 + TT * D3) ) /
                      (1 + TT * ( D4 + TT * (D5 + TT * D6) ) )
    }
    else {
      TT <- PI3 * T3^2
      ALPHA <- (1 + C1 * TT) / (TT * ( 1 + TT * (C2 + TT * C3) ) )
    }
    RTALPH  <- sqrt(ALPHA)
    BETA    <- ROOTPI * L2 * exp( lgamma(ALPHA) - lgamma(ALPHA + 0.5) )
    para[1] <- L1
    para[2] <- BETA * RTALPH
    para[3] <- 2 / RTALPH
    if(lmom$TAU3 < 0) para[3] <- -para[3]
    return(list(type="pe3", para=para, source="parpe3"))
}

