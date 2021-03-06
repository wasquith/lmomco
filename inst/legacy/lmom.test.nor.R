"lmom.test.nor" <-
function(data,digits=4) {
   lmom <- lmom.ub(data)
   para <- parnor(lmom)
   cat("NORMAL DISTRIBUTION PARAMETERS\n")
   str(para)
   lmompara <- lmomnor(para)
   Q50 <- signif(quanor(0.5,para),digits=digits)
   cat(c("Computed median=",Q50,"\n"),sep="")
   P50 <- signif(cdfnor(Q50,para),digits=digits)
   cat(c("Nonexceedance of computed median=",P50,"\n"),sep="")
   lmom.diff(lmompara,lmom,digits=digits)
}

