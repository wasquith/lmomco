"lmom.test.wei" <-
function(data,digits=4) {
   lmom <- lmom.ub(data)
   para <- parwei(lmom)
   cat("WEIBULL DISTRIBUTION PARAMETERS\n")
   str(para)
   lmompara <- lmomwei(para)
   Q50 <- signif(quawei(0.5,para),digits=digits)
   cat(c("Computed median=",Q50,"\n"),sep="")
   P50 <- signif(cdfwei(Q50,para),digits=digits)
   cat(c("Nonexceedance of computed median=",P50,"\n"),sep="")
   lmom.diff(lmompara,lmom,digits=digits)
}

