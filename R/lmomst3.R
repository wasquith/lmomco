"lmomst3" <-
function(para, ...) {
    if(! are.parst3.valid(para)) return()

    U <- para$para[1]
    A <- para$para[2]
    N <- para$para[3]

    SMALL.NU <- 1.001  # arrived from manual experiments involving theoLmoms() testing
    LARGE.NU <- 10^5.5 # arrived from manual experiments involving theoLmoms() testing

    if(N < SMALL.NU) N <- SMALL.NU
    if(N > LARGE.NU) N <- LARGE.NU

    L1   <- U
    L2   <- 2^(6-4*N)*pi*A*sqrt(N)*exp(lgamma(2*N - 2) - 4*lgamma(N/2))
    TAU3 <- 0
    TAU4 <- NA
    TAU5 <- 0
    "poly6" <- function(tau4) {
        tau4 <- log(tau4)
    	   b  <- -2.351846e-06
        ce <- c( 1.493587e+00,  2.272624e-02, -8.763728e-04, -1.529720e-02,
                -1.112732e-02, -7.823726e-03, -3.269371e-03, -4.751368e-04)
        tau6 <- b + ce[1]*tau4^1 + ce[2]*tau4^2 + ce[3]*tau4^3 + ce[4]*tau4^4 +
                    ce[5]*tau4^5 + ce[6]*tau4^6 + ce[7]*tau4^7 + ce[8]*tau4^8
        return(exp(tau6))
    }

    if(N >= LARGE.NU) { # treat as a normal distribution
       TAU4 <- 30/pi * atan(sqrt(2)) - 9
    } else {
       TAU4 <- (15/2)*exp(lgamma(N) - lgamma(1/2) - lgamma(N - 1/2))
       afunc <- function(x) {
          AA <- 1/sqrt(x) * (1 - x)^(N - (3/2)) * pbeta(x, 0.5, N/2)^2
          return(AA)
       }
       int <- NULL
       try(int <- integrate(afunc, 0, 1), silent=TRUE)
       if(is.null(int)) {
          warning("Bad integration")
       } else {
          TAU4 <- TAU4 * int$value - (3/2)
          if(TAU4 >= 1) TAU4 <- 1
       }
    }
    TAU6 <- poly6(TAU4)
    lam <- c(L1, L2,    TAU3*L2, TAU4*L2, TAU5*L2, TAU6*L2)
    rat <- c(NA, L2/L1, TAU3,    TAU4,    TAU5,    TAU6   )
    names(rat) <- NULL
    names(lam) <- NULL
    zz <- list(lambdas=lam, ratios=rat,
               trim=0, leftrim=NULL, rightrim=NULL,
               source="lmomst3")
    return(zz)
}

