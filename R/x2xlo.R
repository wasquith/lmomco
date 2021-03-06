"x2xlo" <-
function(x, leftout=0, a=0, ghost=NULL) {
    pp  <- pp(x, a=a, sort=FALSE)

   xin  <-  x[x >  leftout]
   xlo  <-  x[x <= leftout]
   ppin <- pp[x >  leftout]
   pplo <- pp[x <= leftout]

   if(! is.null(ghost)) {
   	  if(length(x) != length(ghost)) {
   	     warning("Length of x is not the same as the ghosting variable")
   	  }
      ghostin  <-  ghost[x >  leftout]
      ghostlo  <-  ghost[x <= leftout]
   }
   nin <- length(xin); nlo <- length(xlo)
   if(length(xin)  == 0) {
       xin     <- NULL
       ghostin <- NULL
       ppin    <- NULL
       ppthres <- 1
   }
   if(length(xlo)  == 0) {
       xlo     <- NULL
       ghostlo <- NULL
       pplo    <- NULL
       ppthres <- 0
   } else {
       ppthres <- max(pplo)
   }

   z <- list(xin=xin,  ppin=ppin,
             xout=xlo, ppout=pplo,
             pp=ppthres, thres=leftout,
             nin=nin, nlo=nlo, n=nin+nlo,
             source="x2xlo")

   if(! is.null(ghost)) {
      z$ghostin  = ghostin
      z$ghostout = ghostlo
   }
   return(z)
}
