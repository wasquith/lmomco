"lmomkap" <-
function(para, nmom=5) {
   if(! are.parkap.valid(para)) return()
   attributes(para$para) <- NULL
   z <- list(lambdas=rep(NA, nmom), ratios=rep(NA, nmom),
             trim=0, leftrim=NULL, rightrim=NULL, source = "lmomkap")
   r <- 1:nmom; EU <- 0.577215664901532861 # EU IS EULER'S CONSTANT
   # CLUSTER FOR SMALL G
   c1 <- function(G,H) {   #  - CASE H < 0,   G <- 0
        return(EU+log(-H)+digamma( -r/H)) } # Hosking(1994, [12b])
   c2 <- function(G,H) {   #  - CASE H SMALL, G <- 0
        return(EU+log(r)) }                  # Hosking(1994, [12b])
   c3 <- function(G,H) {   #  - CASE H > 0,   G <- 0
        return(EU+log(+H)+digamma(1+r/H)) }  # Hosking(1994, [12b])
   # CLUSTER FOR NONZERO G
   c4 <- function(G,H) {   #  - CASE H < 0,   G NONZERO
        return(exp(lgamma(1+G)+lgamma(-r/H-G) -
                   lgamma(-r/H)-G*log(-H))) }  # Hosking(1994,[12b])
   c5 <- function(G,H) {   #  - CASE H SMALL, G NONZERO
        HOSKING_EXTRA <- (1-0.5*H*G*(1+G)/r)   # Not in Hosking (1994)
        # This extra from Asquith old lmomkap via Hosking (1996) FORTRAN and
        # confirmed in 2005 FORTRAN bundle of Hosking. Check lmom::lmrkap!
        return(exp(lgamma(1+G)-G*log(r))*HOSKING_EXTRA) }
   c6 <- function(G,H) {   #  - CASE H > 0,   G NONZERO
        return(exp(lgamma(1+G)+lgamma(1+r/H) -
                   lgamma(1+G+r/H)-G*log(H))) } # Hosking(1994,[12b])

    U <- para$para[1]; A <- para$para[2]
    G <- para$para[3]; H <- para$para[4]

    # CALCULATE FUNCTIONS OCCURRING IN THE PWM'S BETA-SUB-R
    # SMALL IS USED TO TEST WHETHER H IS EFFECTIVELY ZERO
    SMALL <- 100*.Machine$double.eps # Testing in March 2017 indicates that we
    # can push this to this instead of say sqrt(.Machine$double.eps) without
    # burden of checking for NaN etc. That testing has made during evaluation
    # of an R idiom code nearly 1/2 the code length as here (see below).
    # So it is not quite known if SMALL is optimal but surely a modern R can do
    # better than the 1E-8 in Hosking original code.
    icase <- ifelse(abs(G) < SMALL, 0, 3) # ensembles split with h, test first
    icase <- icase + ifelse(abs(H) < SMALL, 2, ifelse(H < 0, 1, 3))
    BETA <- switch(icase, c1(G,H), c2(G,H), c3(G,H), c4(G,H), c5(G,H), c6(G,H))
    # The author wanted an excuse to use switch()
    ALAM2 <- BETA[2]-BETA[1] # This is setup to Hosking attacking LMR ratios
    z$lambdas[1] <- ifelse(abs(G) < SMALL, U+A*BETA[1], U+A*(1-BETA[1])/G)
    z$lambdas[2] <- ifelse(abs(G) < SMALL, A*ALAM2,       A*ALAM2/(-G))
    z$ratios[2]  <- z$lambdas[2] / z$lambdas[1]

    Zo <- 1 # Begin processing for the upper L-moment ratios
    for(j in 3:nmom) {  # This is ported from Hosking's FORTRAN, hard to see how
      Zo  <- Zo*(4*j-6)/j     # to convert to R idioms. Hosking is clearly using
      Z   <- Zo*3*(j-1)/(j+1) # linear combinations of PWMs with added division
      SUM <- Zo*(BETA[j]-BETA[1])/ALAM2-Z # to get into L-moment ratios. Testing
      if(j == 3) {          # of a process in R idiom (like above caseN's) shows
        z$ratios[j] <- SUM     # concerns in computing the full PWM ensemble and
        z$lambdas[j] <- z$ratios[j] * z$lambdas[2]
      } else {                 # the swapping to L-moments via pwm2lmom() as h
        for(i in 2:(j-2)) {    # because "large". WHA finds by full PWM ensemble
          Z <- Z*(i+i+1)*(j-i)/((i+i-1)*(j+i)) # that about 36 lines are needed
          SUM <- SUM - Z*z$ratios[i+1] # to the say 62 lines for this function.
        }
        z$ratios[j] <- SUM; z$lambdas[j] <- z$lambdas[2] * z$ratios[j]
      }
    }
    return(z)
}

