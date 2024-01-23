"xlo2qua" <-
function(f, para=NULL, xlo=NULL, augasNA=FALSE, sort=FALSE, fillthres=TRUE,
            retrans=function(x) x, paracheck=TRUE, ...) {
   if(is.null(para)) {
     warning("para argument must be specified")
     return(NULL)
   }
   if(is.null(xlo)) par2qua(f, para=para, paracheck=paracheck, ...)
   if(! sort) {
     ff <- (f - xlo$pp) / (1 - xlo$pp) # operation of the f2flo() w/o the overhead and nuances
     ff[ff < 0] <- NA # retain the presence of the slots for the quantiles
     qua <- sapply(ff, function(f) {
                   ifelse(is.na(f), NA, par2qua(f, para=para, paracheck=paracheck, ...))
            })
     if(! is.null(retrans)) qua <- retrans(qua)
     if(fillthres) qua[qua <= xlo$thres] <- xlo$thres; #print(qua)
     if( augasNA ) qua[qua <= xlo$thres] <- NA
     return(qua)
   }

   if(any(is.na(f))) {
     warning("NAs in the nonexceedance probability f, removing them")
     f <- f[! is.na(f)]
   }
   f <- sort(f)
   qua <- par2qua(f2flo(f, xlo=xlo), para=para, paracheck=paracheck, ...)
   n <- length(f); m <- length(qua);                 #print(qua)
   if(! is.null(retrans)) qua <- retrans(qua);      #print(qua)
   if(fillthres) qua[qua <= xlo$thres] <- xlo$thres; #print(qua)
   if(m != n) {
     qua <- c(rep(ifelse(augasNA, "NA", xlo$thres), n-m), sort(qua))
     suppressWarnings(qua <- as.numeric(qua))
     # To silence : In xlo2qua(...) : NAs introduced by coercion
     # Confusing why this showed up and motivated me to use "NA" and not NA
     # in the qua building and explicitly do the numeric conversion as shown
     # as I was trying to isolate what what making R trigger this well known warning.
   }
   if( augasNA ) qua[qua <= xlo$thres] <- NA
   return(qua)
}
