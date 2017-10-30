"lmom.test.exp" <-
function(data,digits=4) {
   lmom <- lmom.ub(data)
   para <- parexp(lmom)
   cat("EXPONENTIAL DISTRIBUTION PARAMETERS\n")
   str(para)
   lmompara <- lmomexp(para)
   Q50 <- signif(quaexp(0.5,para),digits=digits)
   cat(c("Computed median=",Q50,"\n"),sep="")
   P50 <- signif(cdfexp(Q50,para),digits=digits)
   cat(c("Nonexceedance of computed median=",P50,"\n"),sep="")
   lmom.diff(lmompara,lmom,digits=digits)
}

