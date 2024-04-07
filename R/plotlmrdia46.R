"plotlmrdia46" <-
function(lmr=NULL,
      nopoints=FALSE,
      nolines=FALSE,
      noaep4=FALSE,
      nogld=FALSE,
      nopdq4=FALSE,
      nost3=FALSE,
      nosymstable=FALSE,
      notukey=FALSE,
      nonor=FALSE,
      trucate.tau4.to.gtzero=TRUE,
         xlab="L-kurtosis, Tau4, dimensionless",
         ylab="Sixth L-moment ratio, Tau6, dimensionless",
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
      x <- c(-0.25, 1); y <- c(0, 1)
      if(trucate.tau4.to.gtzero) x <- c(0, 1)
      plot(x, y, xlab = xlab, ylab = ylab, type = "n", ...)
      axis(3, at=axTicks(1), labels=NA, lwd=0, lwd.ticks=1, ...)
      axis(4, at=axTicks(2), labels=NA, lwd=0, lwd.ticks=1, ...)
   }
   if(empty) return("empty")

   if(! nolines) {
     if(! noaep4) {
        lines(lmr$aep4, col="red", lty=4, lwd=1*lwd.cex)
        entryi <- entryi + 1
        entries[entryi] <- ifelse(expand.names, "Exponential Power (L-skew = 0)",
                                                "AEP4 (L-skew = 0)")
        Elwd[entryi] <- 1*lwd.cex
        Ecol[entryi] <- "red"
        Epch[entryi] <- NA
        Elty[entryi] <- 4
        Ecex[entryi] <- 1
     }
     if(! nogld) {
        tmp <- lmr$gld[, 2:3]
        if(trucate.tau4.to.gtzero) tmp <- tmp[tmp[,1] >= 0, ]
        lines(tmp, col="purple", lwd=1*lwd.cex)
        entryi <- entryi + 1
        entries[entryi] <- ifelse(expand.names, "Generalized Lambda  (L-scale constant, L-skew = 0)",
                                                "GLD (L-scale constant, L-skew = 0)")
        Elwd[entryi] <- 1*lwd.cex
        Ecol[entryi] <- "purple"
        Epch[entryi] <- NA
        Elty[entryi] <- 1
        Ecex[entryi] <- 1
     }
     if(! nopdq4) {
        tmp <- lmr$pdq4[, 2:3]
        if(trucate.tau4.to.gtzero) tmp <- tmp[tmp[,1] >= 0, ]
        lines(tmp, col="darkgreen", lwd=1*lwd.cex, lty=1)
        entryi <- entryi + 1
        entries[entryi] <- ifelse(expand.names, "Polynomial Quantile Density4", "PDQ4")
        Elwd[entryi] <- 1*lwd.cex
        Ecol[entryi] <- "darkgreen"
        Epch[entryi] <- NA
        Elty[entryi] <- 1
        Ecex[entryi] <- 1
     }
     if(! nost3) {
        lines(lmr$st3, col="blue", lwd=1*lwd.cex, lty=1)
        entryi <- entryi + 1
        entries[entryi] <- ifelse(expand.names, "Student t", "ST3")
        Elwd[entryi] <- 1*lwd.cex
        Ecol[entryi] <- "blue"
        Epch[entryi] <- NA
        Elty[entryi] <- 1
        Ecex[entryi] <- 1
     }
     if(! nosymstable) {
        lines(lmr$symstable[,3:4], col=grey(0.7), lwd=2*lwd.cex, lty=1)
        entryi <- entryi + 1
        entries[entryi] <- ifelse(expand.names, "Symmetric Stable", "Stable")
        Elwd[entryi] <- 2*lwd.cex
        Ecol[entryi] <- grey(0.7)
        Epch[entryi] <- NA
        Elty[entryi] <- 1
        Ecex[entryi] <- 1
     }
     if(! notukey) {
        lines(lmr$tukey[, 3:4], col="purple", lwd=1*lwd.cex, lty=2)
        entryi <- entryi + 1
        entries[entryi] <- ifelse(expand.names, "Tukey Lambda (L-scale not constant)",
                                                "Tukey Lambda")
        Elwd[entryi] <- 1*lwd.cex
        Ecol[entryi] <- "purple"
        Epch[entryi] <- NA
        Elty[entryi] <- 2
        Ecex[entryi] <- 1
     }
   }
   if(! nopoints) {
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
              xjust=0.5, bty="n", cex=legendcex)
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
              xjust=0.5, bty="n", cex=legendcex)
       par(lopts)
     }
   }

   par(popts)
}