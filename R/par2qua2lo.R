"par2qua2lo" <- function(f, para1, para2, xlo1, xlo2,
                            wfunc=NULL, weight=NULL, addouts=FALSE,
                            inf.as.na=TRUE, ...) {

   # Potential early out of the function if the user is asking for the creation of the
   if(! is.null(wfunc) & ! is.function(wfunc) & is.character(wfunc)) { # weight function and not its use
      x1 <- c(xlo1$xin, xlo1$xout); x2 <- c(xlo2$xin, xlo2$xout) # intrasample recombine
      bin <- data.frame(f=pp(c(x1,x2), sort=FALSE), x=c(x1,x2),
                        prob=c(rep(0,length(x1)), rep(1,length(x2))))
      bin <- bin[complete.cases(bin),] # safety feature
      bim <- NULL # this NULL will permit existance testing outside the try()
      opts <- options(warn=-1) # want to silence these messages (tiny sample sizes)
         # Warning: glm.fit: algorithm did not converge
         # Warning: glm.fit: fitted probabilities numerically 0 or 1 occurred
      try(bim <- glm(bin$prob~bin$f, family=binomial), silent=TRUE) # logit is the link
      options(opts)
      if(is.null(bim)) {
         message("could not form the binomial model")
         return(function(f) { return(NA) }) # one why or another
      }
      coe <- coefficients(bim) # get the coefficients and create the weight function
      "wfunc" <- function(f) { z <- exp(coe[1]+coe[2]*f); names(z) <- rep("px2", length(z))
                               attr(z, "coes") <- coe; attr(z, "glm") <- bim
                               attr(z, "data") <- bin; return(z/(z+1)) }
      return(wfunc)
   }

   if(! check.fs(f)) return()
   if(! are.par.valid(para1)) {
      warning("para1 does not appear to hold viable lmomco parameters")
      return()
   }
   if(! are.par.valid(para2)) {
      warning("para2 does not appear to hold viable lmomco parameters")
      return()
   }
   if(xlo1$source != "x2xlo") {
      warning("xlo1 does not appear from x2xlo()")
      return()
   }
   if(xlo2$source != "x2xlo") {
      warning("xlo2 does not appear from x2xlo()")
      return()
   }

   Q1tmp <- par2qua( f2flo(f[f >= xlo1$pp], pp=xlo1$pp), para1, ...)
   Q2tmp <- par2qua( f2flo(f[f >= xlo2$pp], pp=xlo2$pp), para2, ...)
   if(inf.as.na) {
      Q1tmp[! is.finite(Q1tmp)] <- NA
      Q2tmp[! is.finite(Q2tmp)] <- NA
   }

   curve1 <- data.frame(FF=f[f >= xlo1$pp], Q1=Q1tmp)
   curve2 <- data.frame(FF=f[f >= xlo2$pp], Q2=Q2tmp)
   curve12 <- merge(curve1, curve2)

   if(is.null(wfunc)) {
      if(is.null(weight)) {
         t1 <- length(xlo1$xin); t2 <- length(xlo2$xin)
         if(addouts) {
            t1 <- t1+length(xlo1$xout); t2 <- t2+length(xlo2$xout)
         }
         weight <- c(t1, t2) / (t1 + t2)
      } else if(length(weight) == 1) {
         weight <- c((1-weight), weight)
      } else if(length(weight) == 2) {
         if(sum(weight) != 1) {
            warning("sum of the two weights is not unity, going to rescale as such")
            weight <- weight/sum(weight)
         }
      } else {
          warning("weight can not be a vector longer than 2")
          return()
      }
      curve12$Qall <- sapply(1:length(curve12$Q1), function(i) {
                        A <- weight[1]*curve12$Q1[i]
                        B <- weight[2]*curve12$Q2[i]
                        if(is.na(A)) {
                           A <- 0; B <- curve12$Q2[i] # the effective weight of other goes to unity
                        }
                        if(is.na(B)) {
                           B <- 0; A <- curve12$Q1[i] # the effective weight of other goes to unity
                        }
                        return(A+B) })
      zz <- data.frame(f=curve12$FF, quamix=curve12$Qall)
      zz$delta_curve1 <- zz$quamix - curve12$Q1
      zz$delta_curve2 <- zz$quamix - curve12$Q2
      attr(zz, "curve1 weight") <- weight[1]
      attr(zz, "curve2 weight") <- weight[2]
   } else {
      if(! is.function(wfunc)) {
         warning("wfunc is not a function")
         return(NA)
      }
      weights <- wfunc(curve12$FF)
      curve12$Qall <- sapply(1:length(curve12$Q1), function(i) {
                        A <- (1-weights[i])*curve12$Q1[i]
                        B <-     weight[i] *curve12$Q2[i]
                        if(is.na(A)) {
                           A <- 0; B <- curve12$Q2[i] # the effective weight of other goes to unity
                        }
                        if(is.na(B)) {
                           B <- 0; A <- curve12$Q1[i] # the effective weight of other goes to unity
                        }
                        return(A+B) })
      zz <- data.frame(f=curve12$FF, quamix=curve12$Qall)
      zz$delta_curve1 <- zz$quamix - curve12$Q1
      zz$delta_curve2 <- zz$quamix - curve12$Q2
      attr(zz, "curve1 weight") <- (1-weights)
      attr(zz, "curve2 weight") <-    weights
   }
   tmp <- zz; tmp <- tmp[order(tmp$f), ]
   if(any(diff(tmp$quamix) < 0)) warning("result is nonmonotonic increasing")
   return(zz)
}
