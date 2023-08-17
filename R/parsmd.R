"parsmd" <-
function(lmom, checklmom=TRUE, emplims=TRUE, ...) {
    para <- rep(NA, 4)
    names(para) <- c("xi", "a", "b", "q")

    if(length(lmom$L1) == 0) { # convert to named L-moments
      lmom <- lmorph(lmom)     # nondestructive conversion!
    }

    if(checklmom & ! are.lmom.valid(lmom)) {
      warning("L-moments are invalid")
      return()
    }

    # Real world L-moments
    OF   <- lmom$L1 - 1
    MU   <- lmom$L1 - OF

    # Mean zero with L2 of the LCV
    plmom <- list(L1=MU, L2=lmom$L2, L3=lmom$L3, L4=lmom$L4)
    #  print(plmom)
    L1 <- plmom$L1
    L2 <- plmom$L2
    L3 <- plmom$L3
    L4 <- plmom$L4

    uprc <- c(0.16632766,   0.11827093,  0.20727983,   1.21984416,
              2.99333488, -16.28365500, 26.66653095, -19.52994036, 5.44350006)
    lwrc <- c(0.10720630, -0.11626597,   0.78811967,   0.36782311,
              0.51537125, -5.08964587,  11.18233815, -10.20597198, 3.45030275)
    bndtxt <- ""
    if(emplims) {
      Tau3 <- L3 / L2; Tau4 <- L4 / L2
      if(Tau3  <= -0.170) {
        Tau3   <- -0.170
        bndtxt <- "Tau3 <= -0.170, snapped to -0.170; "
      }
      if(Tau3  >=  0.999) {
        Tau3   <-  0.999
        bndtxt <- "Tau3 <= +0.999, snapped to +0.999; "
      }
      if(Tau4  <= 0.104) {  # double assurance knowing that polynomials are to come
        Tau4   <- 0.104
        bndtxt <- "Tau4 <= +0.104, snapped to +0.104; "
      }
      if(Tau4  >=  0.999) { # double assurance knowing that polynomials are to come
        Tau4   <-  0.999
        bndtxt <- "Tau4 <= +0.999, snapped to +0.999; "
      }
      upr <- sum( c( uprc[1], sapply(2:9, function(i) uprc[i] * Tau3^(i-1) ) ) )
      lwr <- sum( c( lwrc[1], sapply(2:9, function(i) lwrc[i] * Tau3^(i-1) ) ) )
      if(Tau4 > upr) {
        L4 <- upr * L2
        bndtxt <- paste0(bndtxt, "Tau4(~Tau3) snapped to upper limit, Tau4=",
                    round(upr, digits=5), " for Tau3=", round(L3/L2, digits=5) )
      }
      #print(c(Tau3, lwr, upr, Tau4))
      if(Tau4 < lwr) {
        L4 <- lwr * L2
        bndtxt <- paste0(bndtxt, "Tau4(~Tau3) snapped to lower limit, Tau4=",
                    round(lwr, digits=5), " for Tau3=", round(L3/L2, digits=5) )
      }
    }

    ofunc <- function(par, L2=NA, L3=NA) {
       B <- exp(par[1]); Q <- exp(par[2])
       IB <- 1/B
       t1 <- exp( lgamma(1*Q - IB) - lgamma(1*Q)  )
       t2 <- exp( lgamma(2*Q - IB) - lgamma(2*Q)  )
       t3 <- exp( lgamma(3*Q - IB) - lgamma(3*Q)  )
       t4 <- exp( lgamma(4*Q - IB) - lgamma(4*Q)  )
       tau3 <- (t1 - 3*t2 +  2*t3       ) / (t1 - t2)
       tau4 <- (t1 - 6*t2 + 10*t3 - 5*t4) / (t1 - t2)
       err <- sqrt((tau3 - L3/L2)^2  + (tau4 - L4/L2)^2)
       #print(c(L3/L2, L4/L2, tau3, tau4, err))
       return(err)
    }

    para <- rep(NA, 4); names(para) <- c("xi", "a", "b", "q")
    para <- list(type="smd", para=para, source="parsmd")
    para.init <- c(1.5, 2.5)
    maxit <- 7
    broken <- NA
    for(i in seq_len(maxit)) {
      rt <- NULL
      try(rt <- optim(log( para.init ), ofunc, L2=L2, L3=L3,
                                         control=list(maxit=1000)), silent=TRUE)
      if(is.null(rt)) next
      if(i < maxit & rt$convergence != 0) next
      B <- exp( rt$par[1] ); Q <- exp( rt$par[2] )
      IB <- 1/B
      t1 <- exp( lgamma(1*Q - IB) - lgamma(1*Q)  )
      t2 <- exp( lgamma(2*Q - IB) - lgamma(2*Q)  )
      A  <- L2 / ( gamma(1+IB) * (t1 - t2) )
      mu <- A * exp(lgamma(1 + IB) + log(t1))
      #print(c(OF, mu-1))
      para <- list(type="smd", para=c(OF - (mu-1), A, B, Q), source="parsmd")
      if(! are.parsmd.valid(para)) break
      lmrsmd <- lmomsmd(para)
      if(! are.lmom.valid(lmrsmd)) break
      #errt1 <- abs(    L1 - lmrsmd$lambdas[1] )
      #errt2 <- abs( L2/L1 - lmrsmd$ratios[ 2] )
      errt3 <- abs( L3/L2 - lmrsmd$ratios[ 3] )
      errt4 <- abs( L4/L2 - lmrsmd$ratios[ 4] )
      #print(c(              errt3, errt4))
      #print(c(errt1, errt2, errt3, errt4))
      if(errt3 < 0.001 & errt4 < 0.001) {
        broken <- TRUE
        break
      } else {
        broken <- FALSE
      }
      para.init <- 10^runif(2, min=-2, max=4)
    }
    para$iter <- i
    para$rt   <- rt
    para$bndtxt <- bndtxt
    ifelse(broken, para$ifail <- 0, para$ifail <- 1)
    return(para)
}
