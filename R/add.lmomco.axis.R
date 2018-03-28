"add.lmomco.axis" <-
function(side=1, twoside=FALSE,
         side.type=c("NPP", "RI", "SNV"), otherside.type=c("NA", "RI", "SNV", "NPP"),
         alt.lab=NA, alt.other.lab=NA, npp.as.aep=FALSE, case=c("upper", "lower"),
         NPP.control=NULL, RI.control=NULL, SNV.control=NULL, ...) {
   case <- match.arg(case)
   other.side <- switch(as.character(side), "1"=3, "2"=4, "3"=1, "4"=1)

   side.type      <- match.arg(side.type)
   otherside.type <- match.arg(otherside.type)
   if(otherside.type == "NA") otherside.type <- NA
   if(twoside & ! is.na(otherside.type)) twoside <- FALSE
   lims <- par()$usr
   ifelse(side == 1 | side == 3, lims <- lims[3:4], lims <- lims[1:2])
   "my.nonexceeds" <- function(minors=FALSE) {
      if(minors) {
         FF <- c(0.55, 0.65, 0.75, 0.825, 0.850, 0.875,
                0.91, 0.92, 0.93, 0.94, 0.96, 0.97,
                0.9925, 0.996, 0.997, 0.9996, 0.9997, 0.99996, 0.99997)
         FF <- c(sort(1-FF), FF)
         FF <- FF[FF > lims[1]]; FF <- FF[FF < lims[2]]
         return(FF)
      } else {
         FF <- c(0.6, 0.7, 0.8, 0.9, 0.95, 0.98, 0.99, 0.995, 0.998, 0.999,
                 0.9995, 0.9998, 0.9999, 0.99995, 0.99998, 0.99999)
         FF <- c(sort(1-FF), 0.5, FF)
         FF <- FF[FF > lims[1]]; FF <- FF[FF < lims[2]]
         return(FF)
      }
   }

   if(is.null(NPP.control)) {
      txt.npp <- ifelse(case == "upper", "NONEXCEEDANCE PROBABILITY", "Nonexceedance probability")
      txt.epp <- ifelse(case == "upper", "EXCEEDANCE PROBABILITY",    "Exceedance probability")

      the.label       <- ifelse(is.na(alt.lab),       txt.npp, alt.lab)
      the.other.label <- ifelse(is.na(alt.other.lab), txt.epp, alt.other.lab)
      NPP.control <- list(label=the.label,
                          other.label=the.other.label,
                          probs=my.nonexceeds(minors=TRUE),
                          probs.label=my.nonexceeds(minors=FALSE),
                          digits=3, line=3, as.exceed=FALSE)
   }
   if(is.null(RI.control)) {
      txt.ri <- ifelse(case == "upper", "RECURRENCE INTERVAL, IN YEARS", "Recurrence interval, in years")
      the.label <- ifelse(is.na(alt.lab), txt.ri, alt.lab)
      RI.control <- list(label=the.label,
                         Tyear=c(2, 5, 10, 25, 50, 100, 200, 500), line=2)
   }
   if(is.null(SNV.control)) {
      txt.snv <- ifelse(case == "upper", "STANDARD NORMAL VARIATE", "Standard normal variate")
      the.label <- ifelse(is.na(alt.lab), txt.snv, alt.lab)
      SNV.control <- list(label=the.label,
                          begin=-5, end=5, by=0.5, line=2)
   }

   NPPf <- function(side, other.side) {
      dots <- list(...)
      tcl <- ifelse("tcl" %in% names(dots), dots$tcl, par()$tcl)
      NPP <- NPP.control$probs;  NPP.lab <- NPP.control$probs.lab
      if(NPP.control$as.exceed) {
         the.true.NPP.lab <- 1 - NPP.lab
      } else {
         the.true.NPP.lab <-     NPP.lab
      }
      qNPP <- qnorm(NPP); qNPP.lab <- qnorm(NPP.lab)
      if(npp.as.aep) {
        NPP.lab <- format(1-the.true.NPP.lab, nsmall=NPP.control$digits)
      } else {
        NPP.lab <- format(  the.true.NPP.lab, nsmall=NPP.control$digits)
      }
      # By placing tcl last (after ...), its value will trump that potentially in ...

      Axis(qNPP,     at=qNPP,     labels=NA,      side=side, ..., tcl=0.8*tcl)
      Axis(qNPP.lab, at=qNPP.lab, labels=NPP.lab, side=side, ..., tcl=1.3*tcl)
      if(npp.as.aep) {
        mtext(NPP.control$other.label, line=NPP.control$line, side=side)
      } else {
        mtext(NPP.control$label, line=NPP.control$line, side=side)
      }

      if(twoside) {
         Axis(qNPP,     at=qNPP,     labels=NA,      side=other.side, ..., tcl=0.8*tcl)
         Axis(qNPP.lab, at=qNPP.lab, labels=NPP.lab, side=other.side, ..., tcl=1.3*tcl)
      }
   }

   RIf <- function(side, other.side) {
      F <- 1 - 1/RI.control$Tyear; qF <- qnorm(F); labF <- RI.control$Tyear
      Axis(qF, at=qF, labels=labF, side=side, ...)
      if(twoside) {
         Axis(at=qF,  labels=NA, side=other.side, ...)
      }
      mtext(RI.control$label, line=RI.control$line, side=side)
   }

   SNVf <- function(side, other.side) {
      SNV <- NULL
      try( SNV <- seq(SNV.control$begin, SNV.control$end, by=SNV.control$by) )
      if(is.null(SNV)) {
         warning("Poorly constructed SNV.control, trapping, and using alternative")
         SNV <- seq(-5, 5, by=0.5)
      }
      Axis(SNV, at=SNV, side=side, ...)
      mtext(SNV.control$label, line=SNV.control$line, side=side)
      if(twoside) {
         Axis(SNV, at=SNV, side=other.side, ...)
      }
   }

   NULLf <- function() { return("no axis function made") }

   primary.axis <- switch(side.type, NPP=NPPf, RI=RIf, SNV=SNVf, NULLf)
   primary.axis(side, other.side)

   if(! is.na(otherside.type)) {
      secondary.axis <- switch(otherside.type, NPP=NPPf, RI=RIf, SNV=SNVf, NULLf)
      secondary.axis(other.side, side)
   }
}

