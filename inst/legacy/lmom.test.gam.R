"lmom.test.gam" <-
function(data,digits=4) {
   lmom <- lmom.ub(data)
   para <- pargam(lmom)
   cat("GAMMA DISTRIBUTION PARAMETERS\n")
   str(para)
   lmompara <- lmomgam(para)
   Q50 <- signif(quagam(0.5,para),digits=digits)
   cat(c("Computed median=",Q50,"\n"),sep="")
   P50 <- signif(cdfgam(Q50,para),digits=digits)
   cat(c("Nonexceedance of computed median=",P50,"\n"),sep="")
   lmom.diff(lmompara,lmom,digits=digits)
}

