"lmom.test.glo" <-
function(data,digits=4) {
   lmom <- lmom.ub(data)
   para <- parglo(lmom)
   cat("GENERALIZED LOGISTIC DISTRIBUTION PARAMETERS\n")
   str(para)
   lmompara <- lmomglo(para)
   Q50 <- signif(quaglo(0.5,para),digits=digits)
   cat(c("Computed median=",Q50,"\n"),sep="")
   P50 <- signif(cdfglo(Q50,para),digits=digits)
   cat(c("Nonexceedance of computed median=",P50,"\n"),sep="")
   lmom.diff(lmompara,lmom,digits=digits)
}

