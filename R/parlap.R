"parlap" <-
function(lmom, checklmom=TRUE, minvar=TRUE, ...) {

    para <- rep(NA,2)
    names(para) <- c("xi","alpha")
    if(length(lmom$lambdas) == 0) { # convert to named L-moments
      lmom <- lmorph(lmom)     # nondestructive conversion!
    }
    if(checklmom & ! are.lmom.valid(lmom)) {
       warning("L-moments are invalid")
       return()
    }

    L1 <- lmom$lambdas[1]
    L2 <- lmom$lambdas[2]
    L3 <- lmom$lambdas[3]
    L4 <- lmom$lambdas[4]

    if(minvar) {
      para[1] <- L1 - (50/31)*L3 # Hosking (1986), IBM RC12210 # 54860, p. 57.
      para[2] <- 1.4741*L2 - 0.5960*L4
      txt <- ":minvar"
    } else {
      para[1] <- L1
      para[2] <- 4/3 * L2
      txt <- ""
    }
    return(list(type='lap', para=para, source=paste0("parlap",txt)))
}


