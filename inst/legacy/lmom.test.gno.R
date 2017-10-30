"lmom.test.gno" <-
function(data,digits=4) {
   lmom <- lmom.ub(data)
   para <- pargno(lmom)
   cat("GENERALIZED NORMAL DISTRIBUTION PARAMETERS\n")
   str(para)
   lmompara <- lmomgno(para)
   Q50 <- signif(quagno(0.5,para),digits=digits)
   cat(c("Computed median=",Q50,"\n"),sep="")
   P50 <- signif(cdfgno(Q50,para),digits=digits)
   cat(c("Nonexceedance of computed median=",P50,"\n"),sep="")
   lmom.diff(lmompara,lmom,digits=digits)
}

