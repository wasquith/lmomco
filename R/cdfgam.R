"cdfgam" <-
function(x,para) {
    if(! are.pargam.valid(para)) return()
    if(length(para$para) == 2) {
       ALPHA <- para$para[1]; BETA  <- para$para[2]
       f <- sapply(1:length(x), function(i) {
                                 if(x[i] <= 0) return(0)
                                 return(pgamma(x[i],ALPHA,scale=BETA)) })
       names(f) <- NULL
       return(f)
    } else if(length(para$para) == 3) {
       MU <- para$para[1]; SIGMA <- para$para[2]; NU <- para$para[3]
       z <- (x/MU)^NU
       theta <- 1/(SIGMA^2*abs(NU)^2) # only for a 2nd check as used in pdfgam
       lGT <- suppressWarnings(lgamma(theta))
       if(! is.finite(lGT) | abs(NU) < 1e-06) {
         # Call taken from inspection of gamlss.dist::pGG and gamlss.dist::pNO
         f <- pnorm(z, mean=log(MU), sd=SIGMA)
       } else {
         # Call taken from inspection of gamlss.dist::pGG and gamlss.dist::pGA
         B <- SIGMA*abs(NU)
         f <- pgamma(z,shape=1/B^2,scale=1*B^2)
       }
       if(NU < 0) f <- 1-f # Notice the flip in probability outside of both
       # the GG and the log-normal. This is not seen in gamlss.dist::pGG where
       # the 1-f does not appear triggered in contrast to code in qGG.
       names(f) <- NULL
       return(f)
    } else {
       stop("should not be here in logic flow")
    }
}
