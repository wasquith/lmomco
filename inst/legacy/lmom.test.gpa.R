"lmom.test.gpa" <-
function(data,digits=4) {
   lmom <- lmom.ub(data)
   para <- pargpa(lmom)
   cat("GENERALIZED PARETO DISTRIBUTION PARAMETERS\n")
   str(para)
   lmompara <- lmomgpa(para)
   Q50 <- signif(quagpa(0.5,para),digits=digits)
   cat(c("Computed median=",Q50,"\n"),sep="")
   P50 <- signif(cdfgpa(Q50,para),digits=digits)
   cat(c("Nonexceedance of computed median=",P50,"\n"),sep="")
   lmom.diff(lmompara,lmom,digits=digits)
}

