"lmom.test.ray" <-
function(data,digits=4) {
   lmom <- lmom.ub(data)
   para <- parray(lmom)
   cat("RAYLEIGH DISTRIBUTION PARAMETERS\n")
   str(para)
   lmompara <- lmomray(para)
   Q50 <- signif(quaray(0.5,para),digits=digits)
   cat(c("Computed median=",Q50,"\n"),sep="")
   P50 <- signif(cdfray(Q50,para),digits=digits)
   cat(c("Nonexceedance of computed median=",P50,"\n"),sep="")
   lmom.diff(lmompara,lmom,digits=digits)
}

