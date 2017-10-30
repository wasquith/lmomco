"lmom.diff" <-
function(lmomparm, lmomdata, verbose=TRUE, digits=4) {

   if(length(lmomparm$L1) == 0) { # convert to named L-moments
     lmomparm <- lmorph(lmomparm)     # nondestructive conversion!
   }
   if(length(lmomdata$L1) == 0) { # convert to named L-moments
     lmomdata <- lmorph(lmomdata)     # nondestructive conversion!
   }

   L1diff <- NA
   L2diff <- NA
   T3diff <- NA
   T4diff <- NA
   T5diff <- NA

   if(! is.null(lmomparm$L1)) {
     L1diff <- lmomparm$L1 - lmomdata$L1
     L1diff <- signif(L1diff,digits=digits)
   }

   if(! is.null(lmomparm$L2)) {
     L2diff <- lmomparm$L2 - lmomdata$L2
     L2diff <- signif(L2diff,digits=digits)
   }

   if(! is.null(lmomparm$TAU3)) {
     T3diff <- lmomparm$TAU3 - lmomdata$TAU3
     T3diff <- signif(T3diff,digits=digits)
   }

   if(! is.null(lmomparm$TAU4)) {
     T4diff <- lmomparm$TAU4 - lmomdata$TAU4
     T4diff <- signif(T4diff,digits=digits)
   }

   if(! is.null(lmomparm$TAU5)) {
     T5diff <- lmomparm$TAU5 - lmomdata$TAU5
     T5diff <- signif(T5diff,digits=digits)
   }

   z <- data.frame(L1diff = L1diff,
                   L2diff = L2diff,
                   T3diff = T3diff,
                   T4diff = T4diff,
                   T5diff = T5diff)

   if(verbose == TRUE) {
     cat("THE FIVE DIFFERENCES BETWEEN L-MOMENTS OF DISTRIBUTION AND DATA\n")
     print(z)
   } 

   return(z)
}

