"f2f" <- function(f, pp=NA, xlo=NULL) {
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
   f  <- f[f >= pp] # see identifical subsetting in f2flo
   return(f)
}