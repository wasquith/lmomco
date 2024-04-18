"plotlmrdia" <-
function(lmr=NULL,
      nopoints=FALSE,
      nolines=FALSE,
      nolimits=FALSE,
      noaep4=FALSE,
      nogev=FALSE,
      noglo=FALSE,
      nogno=FALSE,
      nogov=FALSE,
      nogpa=FALSE,
      nope3=FALSE,
      nopdq3=FALSE,
      nowei=TRUE,
      nocau=TRUE,
      noexp=FALSE,
      nonor=FALSE,
      nogum=FALSE,
      noray=FALSE,
      nosla=TRUE,
      nouni=FALSE,
         xlab="L-skew (Tau3), dimensionless",
         ylab="L-kurtosis (Tau4), dimensionless",
         add=FALSE, empty=FALSE,
         autolegend=FALSE, xleg=NULL, yleg=NULL, legendcex=0.9,
         ncol=1, text.width=NULL, lwd.cex=1, expand.names=FALSE,
         ...) {

   entries <- vector(mode = "character")
   Elwd    <- vector(mode = "numeric")
   Epch    <- vector(mode = "numeric")
   Ecol    <- vector(mode = "numeric")
   Elty    <- vector(mode = "numeric")
   Ecex    <- vector(mode = "numeric")
   entryi  <- 0

   popts <- par(lend=2, no.readonly=TRUE)

   if(is.null(lmr)) empty <- TRUE
   if(! add) {
      plot(c(-1,1), c(-0.25,1), xlab = xlab, ylab = ylab, type = "n", ...)
      axis(3, at=axTicks(1), labels=NA, lwd=0, lwd.ticks=1, ...)
      axis(4, at=axTicks(2), labels=NA, lwd=0, lwd.ticks=1, ...)
   }
   if(empty) return(invisible())

   if(! nolimits) {
     lines(lmr$limits, lwd=2*lwd.cex, col=grey(0.4))
     entryi <- entryi + 1
     entries[entryi] <- "Theoretical limits"
     Elwd[entryi] <- 2*lwd.cex
     Ecol[entryi] <- grey(0.4)
     Epch[entryi] <- NA
     Elty[entryi] <- 1
     Ecex[entryi] <- 1
   }
   if(! nolines) {
     if(! noaep4) {
        lines(lmr$aep4, col="red", lwd=1*lwd.cex, lty=4)
        entryi <- entryi + 1
        entries[entryi] <- ifelse(expand.names, "Asymmetrical Exponential Power4 lower bounds",
                                                "AEP4 lower bounds")
        Elwd[entryi] <- 1*lwd.cex
        Ecol[entryi] <- "red"
        Epch[entryi] <- NA
        Elty[entryi] <- 4
        Ecex[entryi] <- 1
     }
     if(! nogev) {
        lines(lmr$gev, col="darkred", lwd=1*lwd.cex, lty=1)
        entryi <- entryi + 1
        entries[entryi] <- ifelse(expand.names, "Generalized Extreme Value", "GEV")
        Elwd[entryi] <- 1*lwd.cex
        Ecol[entryi] <- "darkred"
        Epch[entryi] <- NA
        Elty[entryi] <- 1
        Ecex[entryi] <- 1
     }
     if(! noglo) {
        lines(lmr$glo, col="green", lwd=1*lwd.cex)
        entryi <- entryi + 1
        entries[entryi] <- ifelse(expand.names, "Generalized Logistic", "GLO")
        Elwd[entryi] <- 1*lwd.cex
        Ecol[entryi] <- "green"
        Epch[entryi] <- NA
        Elty[entryi] <- 1
        Ecex[entryi] <- 1
     }
     if(! nogno) {
        lines(lmr$gno, col="blue", lwd=1*lwd.cex, lty=2)
        entryi <- entryi + 1
        entries[entryi] <- ifelse(expand.names, "Generalized Normal", "GNO")
        Elwd[entryi] <- 1*lwd.cex
        Ecol[entryi] <- "blue"
        Epch[entryi] <- NA
        Elty[entryi] <- 2
        Ecex[entryi] <- 1
     }
     if(! nogov) {
        lines(lmr$gov, col="magenta", lwd=1*lwd.cex, lty=2)
        entryi <- entryi + 1
        entries[entryi] <- ifelse(expand.names, "Govindarajulu", "GOV")
        Elwd[entryi] <- 1*lwd.cex
        Ecol[entryi] <- "magenta"
        Epch[entryi] <- NA
        Elty[entryi] <- 2
        Ecex[entryi] <- 1
     }
     if(! nogpa) {
        lines(lmr$gpa, col="blue", lwd=1*lwd.cex)
        entryi <- entryi + 1
        entries[entryi] <- ifelse(expand.names, "Generalized Pareto", "GPA")
        Elwd[entryi] <- 1*lwd.cex
        Ecol[entryi] <- "blue"
        Epch[entryi] <- NA
        Elty[entryi] <- 1
        Ecex[entryi] <- 1
     }
     if(! nope3) {
        lines(lmr$pe3, col="purple", lwd=1*lwd.cex)
        entryi <- entryi + 1
        entries[entryi] <- ifelse(expand.names, "Pearson Type III" , "PE3")
        Elwd[entryi] <- 1*lwd.cex
        Ecol[entryi] <- "purple"
        Epch[entryi] <- NA
        Elty[entryi] <- 1
        Ecex[entryi] <- 1
     }
     if(! nopdq3) {
        lines(lmr$pdq3, col="darkgreen", lwd=1.3*lwd.cex, lty=2)
        entryi <- entryi + 1
        entries[entryi] <- ifelse(expand.names, "Polynomial Quantile-Density3", "PDQ3")
        Elwd[entryi] <- 1.3*lwd.cex
        Ecol[entryi] <- "darkgreen"
        Epch[entryi] <- NA
        Elty[entryi] <- 2
        Ecex[entryi] <- 1
     }
     if(! nowei) {
        lines(lmr$wei, col="darkorange", lwd=1*lwd.cex, lty=1)
        entryi <- entryi + 1
        entries[entryi] <- ifelse(expand.names, "Weibull", "WEI")
        Elwd[entryi] <- 1*lwd.cex
        Ecol[entryi] <- "darkorange"
        Epch[entryi] <- NA
        Elty[entryi] <- 1
        Ecex[entryi] <- 1
     }
   }
   if(! nopoints) {
     if(! noexp) {
        points(lmr$exp, pch=16, col="red", cex=1.5)
        entryi <- entryi + 1
        entries[entryi] <- ifelse(expand.names, "Exponential", "EXP")
        Elwd[entryi] <- NA
        Ecol[entryi] <- "red"
        Epch[entryi] <- 16
        Elty[entryi] <- NA
        Ecex[entryi] <- 1.5
     }
     if(! nonor) {
        points(lmr$nor, pch=15, col="red", cex=1.5)
        entryi <- entryi + 1
        entries[entryi] <- ifelse(expand.names, "Normal", "NOR")
        Elwd[entryi] <- NA
        Ecol[entryi] <- "red"
        Epch[entryi] <- 15
        Elty[entryi] <- NA
        Ecex[entryi] <- 1.5
     }
     if(! nogum) {
        points(lmr$gum, pch=17, col="red", cex=1.5)
        entryi <- entryi + 1
        entries[entryi] <- ifelse(expand.names, "Gumbel", "GUM")
        Elwd[entryi] <- NA
        Ecol[entryi] <- "red"
        Epch[entryi] <- 17
        Elty[entryi] <- NA
        Ecex[entryi] <- 1.5
     }
     if(! noray) {
        points(lmr$ray, pch=18, col="red", cex=1.5)
        entryi <- entryi + 1
        entries[entryi] <- ifelse(expand.names, "Rayleigh", "RAY")
        Elwd[entryi] <- NA
        Ecol[entryi] <- "red"
        Epch[entryi] <- 18
        Elty[entryi] <- NA
        Ecex[entryi] <- 1.5
     }
     if(! nouni) {
        points(lmr$uniform, pch=12, cex=1.25, col="red")
        entryi <- entryi + 1
        entries[entryi] <- ifelse(expand.names, "Uniform", "UNI")
        Elwd[entryi] <- NA
        Ecol[entryi] <- "red"
        Epch[entryi] <- 12
        Elty[entryi] <- NA
        Ecex[entryi] <- 1.25
     }
     if(! nocau) {
        points(lmr$cau, pch=13, col="turquoise4", cex=1.25)
        entryi <- entryi + 1
        entries[entryi] <- ifelse(expand.names, "CAU (TL1)", "Cauchy (TL1)")
        Elwd[entryi] <- NA
        Ecol[entryi] <- "turquoise4"
        Epch[entryi] <- 13
        Elty[entryi] <- NA
        Ecex[entryi] <- 1.25
     }
     if(! nosla) {
        points(lmr$sla, pch=10, cex=1.25, col="turquoise4")
        entryi <- entryi + 1
        entries[entryi] <- ifelse(expand.names, "Slash (TL1)", "SLA (TL1)")
        Elwd[entryi] <- NA
        Ecol[entryi] <- "turquoise4"
        Epch[entryi] <- 10
        Elty[entryi] <- NA
        Ecex[entryi] <- 1.25
     }
   }

   if(autolegend == TRUE & length(entries) > 0) {
     if(is.character(xleg)) {
       lopts <- par(lend=2, no.readonly=TRUE)
       legend(xleg, entries,
              lwd=Elwd,
              col=Ecol,
              pch=Epch,
              lty=Elty,
              pt.cex=Ecex, ncol=ncol, text.width=text.width,
              xjust=0.5, bty="n", cex=legendcex, ...)
       par(lopts)
     } else {
       if(is.null(xleg)) warning("xleg is NULL, but needed")
       if(is.null(yleg)) warning("yleg is NULL, but needed")
       lopts <- par(lend=2, no.readonly=TRUE)
       legend(xleg, yleg, entries,
              lwd=Elwd,
              col=Ecol,
              pch=Epch,
              lty=Elty,
              pt.cex=Ecex, ncol=ncol, text.width=text.width,
              xjust=0.5, bty="n", cex=legendcex, ...)
       par(lopts)
     }
   }

   par(popts)
}
