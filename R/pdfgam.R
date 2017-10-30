"pdfgam" <-
function(x,para) {
    #if(! are.pargam.valid(para)) return()
    if(length(para$para) == 2) {
       ALPHA <- para$para[1]
       BETA  <- para$para[2]

       f <- dgamma(x, ALPHA, scale=BETA)

       names(f) <- NULL # likely never to be needed here
       f[! is.finite(f)] <- NA
       f[is.na(f)] <- 0 # decision Dec. 2015
       return(f)
    } else if(length(para$para) == 3) {
       # This is the parameterization of the gamlss package
       # f(y|mu,sigma,nu)=theta^theta*z^theta*nu*e^(-theta*z)/(Gamma(theta)*y)
       # where z =(y/mu)^nu, theta = 1/(sigma^2*abs(nu)^2)
       # for y>0, mu>0, sigma>0 and -Inf>nu>Inf. Note that for nu=0 the distribution is log normal.
       MU <- para$para[1]; SIGMA <- para$para[2]; NU <- para$para[3]
       Z <- (x/MU)^NU
       theta <- 1/(SIGMA^2*abs(NU)^2)
       lGT <- suppressWarnings(lgamma(theta))
       #abs(NU) is from gamlss.dist:dGG sources
       if(! is.finite(lGT) | abs(NU) < 1e-06) {
          # taken from gamlss.dist:dGG sources
          lf <- -log(x) - 0.5*log(2*pi) - log(SIGMA) -
                 (1/(2*SIGMA^2)) * (log(x)-log(MU))^2
       } else {
          #f <- theta^theta * Z^theta * NU * exp(-theta*Z) / (GT * x)
          lf <- theta*log(theta) + theta*log(Z) + log(abs(NU)) -
                theta*Z - lGT - log(x)
           # loglik <- theta*log(theta)+theta*log(z)-theta*z-lgamma(theta)+
           #                                            log(abs(nu))-log(x)
           # above is line from gamlss.dist:dGG sources but notice the abs(nu)
           # wrapped by the log.
       }
       f <- exp(lf)
       names(f) <- NULL # likely never to be needed here
       f[! is.finite(f)] <- NA
       f[is.na(f)] <- 0 # decision Dec. 2015
       return(f)
    } else {
       stop("should not be here in logic flow")
    }
}
