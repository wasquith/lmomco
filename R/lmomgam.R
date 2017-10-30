"lmomgam" <-
function(para, ...) {
    z <- list(lambdas=rep(NA,5), ratios= rep(NA,5), source="lmomgam")

    if(length(para$para) == 3) {
       # The upper end point can be problematic for a reasonably robust
       # algorithm for the numerical integrations needed. The while() loop
       # is a scheme to back into the first nonInf upper limit, which is
       # presumably at a very large nonexceedance probability. The L-moment
       # validity check is useful to make a terminal trap on errors.
       eps <- .Machine$double.eps
       xmax <- Inf; maxit <- 100; i <- 0
       while(! is.finite(xmax)) {
         xmax <- quagam(1-eps, para)
         i <- i + 1; eps <- eps^0.5
         if(i > maxit) break
       }
       lmr <- theoLmoms.max.ostat(para=para, cdf=cdfgam, pdf=pdfgam,
                                  lower=0, upper=xmax, ...)
       lmr$source <- "lmomgam"
       if(! are.lmom.valid(lmr)) {
          warning("integrated L-moments test as invalid")
          return(NULL)
       }
       return(lmr)
    }

    # Note that TAU5 and L5 are not available from Hosking's FORTRAN base.

    # CONST IS 1/sqrt(pi)
    CONST <- 0.564189583547756287

    # COEFFICIENTS OF RATIONAL-FUNCTION APPROXIMATIONS
    #   A0 IS 1/sqrt(3*pi)
    #   C0 IS TAU-4 FOR THE NORMAL DISTRIBUTION
    A0 <- 0.32573501
    A1 <- 0.16869150;     A2 <- 0.78327243e-01; A3 <- -0.29120539e-02
    B1 <- 0.46697102;     B2 <- 0.24255406
    C0 <- 0.12260172;
    C1 <- 0.53730130e-01; C2 <- 0.43384378e-01; C3 <- 0.11101277e-01
    D1 <- 0.18324466;     D2 <- 0.20166036e+00
    E1 <- 0.23807576e+01; E2 <- 0.15931792e+01; E3 <- 0.11618371e+00
    F1 <- 0.51533299e+01; F2 <- 0.71425260e+01; F3 <- 0.19745056e+01
    G1 <- 0.21235833e+01; G2 <- 0.41670213e+01; G3 <- 0.31925299e+01
    H1 <- 0.90551443e+01; H2 <- 0.26649995e+02; H3 <- 0.26193668e+02

    if(! are.pargam.valid(para)) return()
    attributes(para$para) <- NULL

    ALPHA <- para$para[1]
    BETA  <- para$para[2]
    z$lambdas[1] <- ALPHA*BETA
    z$lambdas[2] <- BETA*CONST*exp(lgamma(ALPHA + 0.5) -
                                   lgamma(ALPHA))
    if(ALPHA < 1) {
      Z <- ALPHA
      z$ratios[3] <- (((E3*Z+E2)*Z+E1)*Z+1)/(((F3*Z+F2)*Z+F1)*Z+1)
      z$ratios[4] <- (((G3*Z+G2)*Z+G1)*Z+1)/(((H3*Z+H2)*Z+H1)*Z+1)
    }
    else {
      Z <- 1/ALPHA
      z$ratios[3] <- sqrt(Z)*(((A3*Z+A2)*Z+A1)*Z+A0) /
                                    ((B2*Z+B1)*Z+1)
      z$ratios[4] <- (((C3*Z+C2)*Z+C1)*Z+C0)/((D2*Z+D1)*Z+1)
    }
    z$ratios[2]  <- z$lambdas[2]/z$lambdas[1]
    z$lambdas[3] <- z$ratios[3]*z$lambdas[2]
    z$lambdas[4] <- z$ratios[4]*z$lambdas[2]
    return(z)
}

