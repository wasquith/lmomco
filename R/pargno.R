"pargno" <-
function(lmom, checklmom=TRUE, useHosking=TRUE, ...) {
    para <- rep(NA, 3)
    names(para) <- c("xi", "alpha", "kappa")
    "erf" <- function(x) 2 * pnorm(x * sqrt(2)) - 1

    if(length(lmom$L1) == 1) { # convert to named L-moments
      lmom <- lmorph(lmom)     # nondestructive conversion!
    }
    if(checklmom & ! are.lmom.valid(lmom)) {
       warning("L-moments are invalid")
       return()
    }
    T3 <- lmom$ratios[3]; SGN <- sign(T3)

    if(! useHosking) {
      # Coefficients from the embedded example in pargno.Rd from which refitting to numerically
      # integrated Tau3 from theoLmoms() for a given gamma is used with curation of settings for
      # integrate() and optim() call to estimate other A and B parameters (also using a new SMALL).
      # However, we use Hosking original structure and use Hosking's parameters as  initial starting
      # points in a 7D optim()ization. These coefficients result from the script
      # lmomco/inst/newcoes/makeGNOcoes.R.
      SMALL <- .Machine$double.eps
      A0 <-  2.05861903000273; A1 <- -3.722332886770820
      A2 <-  1.81265195204820; A3 <- -0.142649680677246
      B1 <- -2.02334708603110; B2 <-  1.195993222763030; B3 <- -0.171431551954349

      if(abs(T3) > 0.999) {
        warning("L-skew is too large, truncating L-skew to 0.999")
        T3 <- SGN*0.999
        #para[1] = 0
        #para[2] = -1
        #para[3] = 0
        #return(list(type = 'gno', para = para))
      }
    } else {
      #  SMALL IS USED TO TEST WHETHER SKEWNESS IS EFFECTIVELY ZERO
      SMALL <- 1e-8
      # METHOD: RATIONAL-FUNCTION APPROXIMATION OF K IN TERMS OF TAU-3
      #  COEFFICIENTS OF RATIONAL-FUNCTION APPROXIMATION
      #  A0 IS 0.5*sqrt(3/pi)
      A0 <-  0.20466534e1;   A1 <- -0.36544371e+1;
      A2 <-  0.18396733e+1;  A3 <- -0.20360244;
      B1 <- -0.20182173e+1;  B2 <-  0.12420401e+1;  B3 <- -0.21741801

      if(abs(T3) > 0.95) {
        warning("L-skew is too large, truncating L-skew to 0.95")
        T3 <- SGN*0.95
        #para[1] = 0
        #para[2] = -1
        #para[3] = 0
        #return(list(type = 'gno', para = para))
      }
    }


    if(abs(T3) <= SMALL) {
      para[1] <- lmom$lambdas[1]
      para[2] <- lmom$lambdas[2]*sqrt(pi)
      para[3] <- 0
      return(list(type="gno", para=para))
    } else {
      TT <-  T3 * T3
      K  <- -T3 * (A0 + TT * (A1 + TT * (A2 + TT * A3) ) ) /
                   (1 + TT * (B1 + TT * (B2 + TT * B3) ) )
      E  <- exp(0.5 * K * K)
      A  <- lmom$lambdas[2] * K / (E * erf(0.5 * K))
      XI <- lmom$lambdas[1] + A * (E-1) / K
      para[1] <- XI
      para[2] <- A
      para[3] <- K
      return(list(type="gno", para=para, source="pargno"))
    }
}

