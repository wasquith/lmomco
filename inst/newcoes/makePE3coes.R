if(Sys.getenv("RSTUDIO") == "1") {
  # Automatic change directory to location of this script when RStudio is running.
  setwd(dirname(rstudioapi::getSourceEditorContext()$path))
}

library(lmomco)
  set.seed(1)
  MU <- 0; SD <- 1      # mean of zero and standard deviation of 1
  GM <- 10^seq(-3, log10(125), by=0.001) # given gamma parameter values
  GM <- sort( c(0, GM)) # with a 152 heuristic limit and also include zero
  LM <- NULL            # data frame built by numerical integration
  for(gm in GM) {
    lmr  <- theoLmoms(list(para=c(MU, SD, gm), type="pe3"), nmom=4,
                      subdivisions=200, rel.tol=.Machine$double.eps^0.50)
    LM   <- rbind(LM, data.frame(mu=0, sd=1, gamma=gm,
                                 L1=lmr$lambdas[1], L2=lmr$lambdas[2],
                                 T3=lmr$ratios[ 3], T4=lmr$ratios[ 4]))
  }

  LM <- LM[abs(LM$T3) < 1,]
  LM <- LM[abs(LM$T4) < 1,]
  plot(LM$gamma, LM$T3, type="l")

  plot(log10(LM$gamma), LM$T4)

  # from Hosking's FORTRAN code in 1990s (C[1-3]; D[1-6])
  CoesA <- c(0.2906,   0.1882,  0.0442)
  DoesA <- c(0.36067, -0.59567, 0.25361, -2.78861, 2.56096, -0.77045)

  Cpara.init <- CoesA
  Dpara.init <- DoesA
  para.init  <- c(Cpara.init, Dpara.init)

  "tau3gamma" <-
  function(tau3, C=Cpara.init, D=Dpara.init,
           SMALL=.Machine$double.eps) {
     sg <- sign(  tau3 )
     gg <- sapply(tau3, function(t3) {
         T3 <- abs(t3)
         if(T3 <= SMALL) return(0)
         if(T3 >= 1/3) {
           TT <- 1 - T3
           ALPHA <- TT * (D[1] + TT * (D[2] + TT * D[3])) /
               (1 + TT * (D[4] + TT * (D[5] + TT * D[6])))
         } else {
           TT <- 3 * pi * T3^2
           ALPHA <-             (1 + C[1] * TT) /
                    (TT * (1 + TT * (C[2] + TT * C[3])))
       }
       suppressWarnings(RTALPH <- sqrt(ALPHA))
       gamma <- 2 / RTALPH
       return(gamma) })
    return(sg * gg)
  }

  "ofunc" <- function(para, t3=NA, tgamma=NA) {
    gg <- tau3gamma(t3, C=para[1:3], D=para[4:9])
    return(sum((gg - tgamma)^2)) # sum of squared error
  }

  rt <- optim(para.init, fn=ofunc, t3=LM$T3, tgamma=LM$gamma,
              control=list(reltol=.Machine$double.eps^0.5, maxit=10000))
  CoesB <- rt$par[1:3]; DoesB <- rt$par[4:9]

  GammaHosk <- LM$gamma_hosking <- tau3gamma(LM$T3)
  GammaMade <- LM$gamma_lmomco  <- tau3gamma(LM$T3, C=rt$par[1:3], D=rt$par[4:9])

  plot(LM$gamma, GammaMade, lwd=0.7, pch=16, cex=0.5, log="xy",
       xlab="Given PE3 gamma", ylab="Gamma by new coefficients")
  abline(0,1, col="red")

  plot(LM$gamma, abs(GammaHosk - GammaMade), lwd=0.7, pch=16, cex=0.5,
       xlab="Given theoretical gamma", ylab="abs(Gamma difference)")
  plot(LM$T3, abs(GammaHosk - GammaMade), lwd=0.7, pch=16, cex=0.5,
       xlab="Given theoretical T3 to gamma", ylab="abs(Gamma difference)")

  message("Hosking Coes: ", paste(CoesA, collapse=", "))
  message("optim() Coes: ", paste(CoesB, collapse=", "))
  message("Hosking Does: ", paste(DoesA, collapse=", "))
  message("optim() Does: ", paste(DoesB, collapse=", "))

  sqrt(mean((GammaHosk-LM$gamma)^2)) # 0.0005561338
  sqrt(mean((GammaMade-LM$gamma)^2)) # 0.0004815001
  100 * (0.0005561338 - 0.0004815001)/abs(0.0004815001)
