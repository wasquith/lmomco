"f2flo" <-
function(f, pp=NA, xlo=NULL, type=c("ge", "gt")) {
   if(! check.fs(f)) return(FALSE)
   if(! is.null(xlo)) pp <- xlo$pp
   if(is.na(pp)) {
      warning("pp can not be NA")
      return(FALSE)
   } else {
      if(pp < 0 || pp > 1) {
        print("pp argument is not a valid nonexceedance probability")
        return(FALSE)
      }
   }
   type <- match.arg(type)
   if(type == "gt") {
     f <- f[f >  pp] # see identifical subsetting in f2f
   } else {
     f <- f[f >= pp] # see identifical subsetting in f2f
   }
   zs <- (f-pp)/(1-pp)
   if(any(zs < 0) || any(zs > 1)) {
      warning("invalid nonexceedance probability after pp conditioning")
      return(FALSE)
   }
   return(zs)
}
