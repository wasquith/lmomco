"lmom.test.gev" <-
function(data,digits=4) {
   lmom <- lmom.ub(data)
   para <- pargev(lmom)
   cat("GENERALIZED EXTREME VALUE DISTRIBUTION PARAMETERS\n")
   str(para)
   lmompara <- lmomgev(para)
   Q50 <- signif(quagev(0.5,para),digits=digits)
   cat(c("Computed median=",Q50,"\n"),sep="")
   P50 <- signif(cdfgev(Q50,para),digits=digits)
   cat(c("Nonexceedance of computed median=",P50,"\n"),sep="")
   lmom.diff(lmompara,lmom,digits=digits)
}

