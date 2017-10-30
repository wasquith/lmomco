"pargum" <-
function(lmom,checklmom=TRUE,...) {
   # euler <- print(-digamma(1), digits=15)
   euler <- 0.5772156649015323 # Euler's constant
   para <- rep(NA,2)
   names(para) <- c("xi","alpha")
   if(length(lmom$L1) == 1) { # convert to named L-moments
     lmom <- lmorph(lmom)     # nondestructive conversion!
   }
   if(checklmom & ! are.lmom.valid(lmom)) {
     warning("L-moments are invalid")
     return()
   }
   para[2] <- lmom$lambdas[2]/log(2)
   para[1] <- lmom$lambdas[1]-euler*para[2]
   return(list(type = 'gum', para=para, source="pargum"))
}

