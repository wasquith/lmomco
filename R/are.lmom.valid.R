"are.lmom.valid" <-
function(lmom, checkt3t4=TRUE) {
   if(is.null(lmom)) return(FALSE)
   if(length(lmom$L1) == 0) { # convert to named L-moments
     lmom <- lmorph(lmom)     # nondestructive conversion!
   }

   # It was a design mistake to use named L-moments. Can not query the length
   # say of lambdas to ascertain which L-moments are being pursued by the
   # user.   Late (2017) additions are finiteness check on the mean and L2. The
   # L2 check should catch other L-moments from being infinite. The finiteness
   # checks do provide protection against say the lmoms() function being called
   # on logarithms of the data and the data having <= zero values.

   # The early return trues are for situations in which the higher moments
   # are simply not available--say from computing the l-moments of a distribution
   if(! is.finite(lmom$L1))                         return(FALSE)
   if(    is.null(lmom$L2)   | is.na(lmom$L2))      return(TRUE )
   if(! is.finite(lmom$L2)   |       lmom$L2 <= 0)  return(FALSE)
   if(    is.null(lmom$TAU3) | is.na(lmom$TAU3))    return(TRUE )
   if(        abs(lmom$TAU3) >= 1)                  return(FALSE)
   if(    is.null(lmom$TAU4) | is.na(lmom$TAU4))    return(TRUE )
   if(            lmom$TAU4  >= 1 )                 return(FALSE)
   if(checkt3t4) { # This test must come before a Tau5 check! (03/07/2017)
      if(lmom$TAU4 < (5*lmom$TAU3^2 - 1)/4)         return(FALSE)
   }
   if(    is.null(lmom$TAU5) | is.na(lmom$TAU5))    return(TRUE )
   if(        abs(lmom$TAU5) >= 1)                  return(FALSE)
   return(TRUE)
}

