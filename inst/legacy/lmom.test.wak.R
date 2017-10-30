"lmom.test.wak" <-
function(data,digits=4) {
   lmom <- lmom.ub(data)
   para <- parwak(lmom)
   cat("WAKEBY DISTRIBUTION PARAMETERS\n")
   str(para)
   if(para$ifail > 0) {
     warning("early exit")
     return()
   }
   lmompara <- lmomwak(para)
   Q50 <- signif(quawak(0.5,para),digits=digits)
   cat(c("Computed median=",Q50,"\n"),sep="")
   P50 <- signif(cdfwak(Q50,para),digits=digits)
   cat(c("Nonexceedance of computed median=",P50,"\n"),sep="")
   lmom.diff(lmompara,lmom,digits=digits)
}

