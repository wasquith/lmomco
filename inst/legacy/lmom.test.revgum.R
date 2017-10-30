"lmom.test.revgum" <-
function(data,zeta=1,digits=4) {
   lmom <- lmom.ub(data)
   para <- parrevgum(lmom,zeta=zeta)
   cat("REVERSE GUMBEL DISTRIBUTION PARAMETERS\n")
   lmompara <- lmomrevgum(para)
   Q50 <- signif(quarevgum(0.5,para),digits=digits)
   cat(c("Computed median=",Q50,"\n"),sep="")
   P50 <- signif(cdfrevgum(Q50,para),digits=digits)
   cat(c("Nonexceedance of computed median=",P50,"\n"),sep="")
   lmom.diff(lmompara,lmom,digits=digits)
}
