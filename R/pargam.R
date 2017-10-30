"pargam" <-
function(lmom, p=c("2", "3"), checklmom=TRUE,...) {

  if(length(lmom$lambdas) == 0) { # convert to named L-moments
     lmom <- lmorph(lmom)         # nondestructive conversion!
  }
  if(checklmom & ! are.lmom.valid(lmom)) {
     warning("L-moments are invalid")
     return()
  }
  LL <- lmom$lambdas
  p <- as.character(p)
  p <- as.numeric(match.arg(p))
  if(p == 2) {
    if(length(LL) < 2) {
       warning("not enough L-moments (need 2)")
       return(NULL)
    }
    para <- rep(NA,2)
    names(para) <- c("alpha","beta")
    # METHOD: RATIONAL APPROXIMATION IS USED TO EXPRESS ALPHA AS A FUNCTION
    # OF L-CV. RELATIVE ACCURACY OF THE  APPROXIMATION IS BETTER THAN 5E-5.
    #
    #  CONSTANTS USED IN MINIMAX APPROXIMATIONS
    #
    A1 <- -0.3080; A2 <- -0.05812; A3 <-  0.01765
    B1 <-  0.7213; B2 <- -0.5947;  B3 <- -2.1817; B4 <- 1.2113
    L1 <- LL[1]; LCV <- LL[2]/LL[1]
    if(LCV >= 0.5) {
      TT <- 1-LCV
      ALPHA <- TT*(B1+TT*B2)/(1+TT*(B3+TT*B4))
    }
    else {
      TT <- pi*LCV^2
      ALPHA <- (1+A1*TT)/(TT*(1+TT*(A2+TT*A3)))
    }
    para[1] <- ALPHA
    para[2] <- L1/ALPHA
    z <- list(type = "gam", para = para, source="pargam")
    if(are.pargam.valid(z)) {
      return(z)
    }
    else {
      warning("Parameters can not be computed likely because ",
              "L1 <= L2 or L2 <= 0")
      return(NULL)
    }
  } else if(p == 3) {
    if(length(LL) < 3) {
       warning("not enough L-moments (need 3)")
       return(NULL)
    }
    COE <- c(+3.196931e+00,
     +7.305355e-03, -5.528250e-06, +1.915532e-09, -2.937845e-13, +1.621107e-17,
     -1.351629e+01, +1.951213e+01, -1.377507e+01, +4.292661e+00, -5.062394e-01 )
    ETA <- sqrt(log(LL[1]/LL[2])) # NOTE L1/L2 not LCV!
    PWR <- c(0, -1,-2,-3,-4,-5, +1,+2,+3,+4,+5)
    lSIG <- sum(COE*ETA^PWR)
    if(is.nan(lSIG)) lSIG <- 8
    SIG <- exp(lSIG)
    para.intA <- c(1,.2) # log(MU), NU
    objfuncA <- function(k) {
       tmp <- vec2par(c(exp(k[1]),SIG,k[2]), type="gam")
       if(is.null(tmp)) return(Inf)
       tLL <- lmomgam(tmp); if(is.null(tLL)) return(Inf)
       return(sum((tLL$lambdas[1:3]-LL[1:3])^2))
    }
    objfuncB <- function(k) {
       tmp <- vec2par(c(exp(k[1]),exp(k[2]),k[3]), type="gam")
       if(is.null(tmp)) return(Inf)
       tLL <- lmomgam(tmp); if(is.null(tLL)) return(Inf)
       return(sum((tLL$lambdas[1:3]-LL[1:3])^2))
    }

    rtA <- NULL
    try(rtA <- optim(para.intA, objfuncA), silent=TRUE)
    if(is.null(rtA)) {
       warning("optim() attempt A is NULL"); return(NULL)
    }
    paraA <- c(exp(rtA$par[1]), SIG, rtA$par[2],
               rtA$value, rtA$counts[1],  rtA$convergence)
    names(paraA) <- c("mu", "sigma", "nu",
                      "optimAvalue", "optimAcounts", "optimAconvergence")
    para.intB <- c(rtA$par[1], lSIG, rtA$par[2])
    rtB <- NULL
    try(rtB <- optim(para.intB, objfuncB), silent=TRUE)
    if(is.null(rtB)) {
       warning("optim() attempt B is NULL"); return(NULL)
    }
    paraB <- vec2par(c(exp(rtB$par[1]), exp(rtB$par[2]), rtB$par[3]), type="gam")
    paraB$paraA <- paraA
    paraB$optim <- rtB
    return(paraB)
  } else {
    stop("should not be here in logic")
  }
}

