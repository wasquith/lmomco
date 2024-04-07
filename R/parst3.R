"parst3" <-
function(lmom, checklmom=TRUE, ...) {
    para <- rep(NA, 3)
    names(para) <- c("xi", "alpha", "nu")

    if(length(lmom$lambdas) == 0) { # convert to named L-moments
      lmom <- lmorph(lmom)     # nondestructive conversion!
    }

    if(checklmom & ! are.lmom.valid(lmom)) {
      warning("L-moments are invalid")
      return()
    }

    if(length(lmom$ratios) <= 3) {
      warning("Not enough L-moments, need Tau4 to fit Student t (3-parameter)")
      return()
    }

    # theoLmoms(list(para=c(0,1,1.001), type="st3"), nmom=6, verbose=TRUE)$ratios[4]
    # Tau4 with integration is 0.998167


    SMALL.NU <- 1.001  # arrived from manual experiments involving theoLmoms() testing
    LARGE.NU <- 10^5.5 # arrived from manual experiments involving theoLmoms() testing
    NUDGE <- 0.0000002 # to ensure that parst3(vec2lmom(c(10, 2, 0, 0.1226)))$para or
    # otherwise Tau4 near to exactly the normal still uniroots()
    HIGHEST.TAU4 <- 0.998

    TAU4.NORMAL <- 30/pi * atan(sqrt(2)) - 9

    TAU4 <- lmom$ratios[4]
    if(TAU4 >= 0.1226 & TAU4 <= TAU4.NORMAL) TAU4 <- TAU4.NORMAL + NUDGE
    if(TAU4 <  0.1226) {
      warning("TAU4 is less than the lowest shorthand TAU4 '0.1226' for the normal distribution")
      return(NULL)
    }
    if(TAU4 >= HIGHEST.TAU4) TAU4 <- HIGHEST.TAU4

    # This is just a bit expensive to waste CPU on the r != 2,4 L-moments for the distribution.
    # But presumably more accurate and simpler code to uniroot for the Tau4 and then back-
    # substitution for the rest of the parameters in lieu of using optim().
    "ofunc" <- function(nu, tau4=NA) {
       val <- theoTLmoms(list(para=c(0, 1, nu), type="st3"), nmom=4)$ratios[4]
       return(val - tau4)
    }
    rt <- NULL; N <- NA
    try( rt <- uniroot(ofunc, interval=c(SMALL.NU, LARGE.NU), tau4=TAU4), silent=TRUE)
    #print(rt)
    if(is.null(rt)) {
      warning("Could not root solution of Tau4 as a function of Nu")
      return(NULL)
    } else {
      N <- rt$root
    }

    denom <- 2^(6-4*N) * pi * sqrt(N) * exp(lgamma(2*N-2) - 4*lgamma(N/2))
    if(! is.finite(denom) | is.nan(denom)) {
      denom <- 1/sqrt(pi) # from experiments and limiting arguments
    }

    A <- lmom$lambdas[2] / denom
    U <- lmom$lambdas[1]
    para[1:3] <- c(U, A, N)
    return(list(type='st3', para=para, rt=rt, source="parst3"))
}

