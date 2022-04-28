"pdfgev" <-
function(x, para, paracheck=TRUE) {
   if(paracheck) {
     if(! are.pargev.valid(para)) return()
   }
   U <- para$para[1]
   A <- para$para[2]
   K <- para$para[3]

   Y  <- (x - U)/A
   ZERO <- sqrt(.Machine$double.eps)
   if(abs(K) > ZERO) {
     Y <- suppressWarnings(-log(1-K*Y)/K)
   }
   f <- A^(-1) * exp(-(1-K)*Y-exp(-Y))

   if(abs(K) > ZERO) {
     if(K > 0) {
       f[x > U + A/K] <- 0
     } else {
       f[x < U + A/K] <- 0
     }
   }

   names(f) <- NULL
   #f[! is.finite(f)] <- NA
   #f[is.na(f)] <- 0 # decision Dec. 2015
   return(f)
}

# Disabled version with bad support for X, reported by
# Christophe Dutang during comparisons of GEV implementations
# in R and various packages in Spring 2022.
# "pdfgev2" <-
# function(x,para) {
#     if(! are.pargev.valid(para)) return()
#     XI <- para$para[1]
#     A  <- para$para[2]
#     K  <- para$para[3]
#
#     Y  <- (x - XI)/A
#     if(K != 0) {
#        ARG <- 1-K*Y
#        Y <- suppressWarnings(-log(ARG)/K)
#     }
#     f <- A^(-1) * exp(-(1-K)*Y-exp(-Y))
#
#     names(f) <- NULL
#     f[! is.finite(f)] <- NA
#     f[is.na(f)] <- 0 # decision Dec. 2015
#     return(f)
# }
