"xlo2qua" <-
function(f, para=NULL, xlo=NULL, augasNA=FALSE,
            retrans=function(x) x, paracheck=TRUE, ...) {
   if(any(is.na(f))) {
     warning("NAs in the nonexceedance probability f, removing them")
   }
   f <- sort(f)
   qua <- par2qua(f2flo(f, xlo=xlo), para=para, paracheck=TRUE, ...)
   n <- length(f); m <- length(qua); #print(qua)
   qua <- retrans(qua);                #print(qua)
   qua[qua <= xlo$thres] <- xlo$thres; #print(qua)
   if(m != n) {
     qua <- c(rep(ifelse(augasNA, "NA", xlo$thres), n-m), sort(qua))
     suppressWarnings(qua <- as.numeric(qua))
     # To silence : In xlo2qua(...) : NAs introduced by coercion
     # Confusing why this showed up and motiviated me to use "NA" and not NA
     # in the qua building and explicity do the numeric conversion as shown
     # as I was trying to isolate what what making R trigger this well known warning.
   }
   return(qua)
}
