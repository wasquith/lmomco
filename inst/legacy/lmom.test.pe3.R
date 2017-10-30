"lmom.test.pe3" <-
function(data,digits=4) {
   lmom <- lmom.ub(data)
   para <- parpe3(lmom)
   cat("PEARSON TYPE III DISTRIBUTION PARAMETERS\n")
   str(para)
   lmompara <- lmompe3(para)
   Q50 <- signif(quape3(0.5,para),digits=digits)
   cat(c("Computed median=",Q50,"\n"),sep="")
   P50 <- signif(cdfpe3(Q50,para),digits=digits)
   cat(c("Nonexceedance of computed median=",P50,"\n"),sep="")
   lmom.diff(lmompara,lmom,digits=digits)
}

