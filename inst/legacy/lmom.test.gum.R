"lmom.test.gum" <-
function(data,digits=4) {
   lmom <- lmom.ub(data)
   para <- pargum(lmom)
   cat("GUMBEL DISTRIBUTION PARAMETERS\n")
   str(para)
   lmompara <- lmomgum(para)
   Q50 <- signif(quagum(0.5,para),digits=digits)
   cat(c("Computed median=",Q50,"\n"),sep="")
   P50 <- signif(cdfgum(Q50,para),digits=digits)
   cat(c("Nonexceedance of computed median=",P50,"\n"),sep="")
   lmom.diff(lmompara,lmom,digits=digits)
}

