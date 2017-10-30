"quagam" <-
function(f,para,paracheck=TRUE) {
    if(! check.fs(f)) return()
    if(paracheck == TRUE) {
      if(! are.pargam.valid(para)) return()
    }
    if(length(para$para) == 2) {
       ALPHA <- para$para[1]
       BETA  <- para$para[2]

       x <- qgamma(f,ALPHA,scale=BETA)
       names(x) <- NULL
       return(x)
    } else if(length(para$para) == 3) {
       MU <- para$para[1]; SIGMA <- para$para[2]; NU <- para$para[3]
       if(NU < 0) f <- 1-f # Notice the flip in probability outside of both
       # the GG and the log-normal. This is not seen in gamlss.dist::qGG where
       # the 1-f is commented out for the log-Normal (for a nonvectorized nu)
       # only for a 2nd check as used in pdfgam
       theta <- 1/(SIGMA^2*abs(NU)^2)
       lGT <- suppressWarnings(lgamma(theta))
       if(! is.finite(lGT) | abs(NU) < 1e-06) {
         # Call taken from inspection of gamlss.dist::qGG and gamlss.dist::qNO
         z <- qnorm(f, mean=log(MU), sd=SIGMA)
         x <- exp(z)
       } else {
         # Call taken from inspection of gamlss.dist::qGG and gamlss.dist::qGA
         B <- SIGMA*abs(NU)
         z <- qgamma(f, shape=1/B^2, scale=1*B^2)
         x <- MU*z^(1/NU)
       }
       names(x) <- NULL
       return(x)
    } else {
       stop("should not be here in logic flow")
    }
}

