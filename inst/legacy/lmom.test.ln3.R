"lmom.test.ln3" <-
function(data,digits=4) {
   lmom <- lmom.ub(data)
   para <- parln3(lmom)
   cat("LOG-NORMAL (3 parameter) DISTRIBUTION PARAMETERS\n")
   str(para)
   lmompara <- lmomln3(para)
   Q50 <- signif(qualn3(0.5,para),digits=digits)
   cat(c("Computed median=",Q50,"\n"),sep="")
   P50 <- signif(cdfln3(Q50,para),digits=digits)
   cat(c("Nonexceedance of computed median=",P50,"\n"),sep="")
   lmom.diff(lmompara,lmom,digits=digits)
}

