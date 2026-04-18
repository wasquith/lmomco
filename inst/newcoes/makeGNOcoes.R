if(Sys.getenv("RSTUDIO") == "1") {
  # Automatic change directory to location of this script when RStudio is running.
  setwd(dirname(rstudioapi::getSourceEditorContext()$path))
}

library(lmomco)
  set.seed(1)
  MU <- 0; SD <- 1      # first parameter == 0 and second parameter == 1
  KK <- 10^seq(-3, log10(5), by=0.001) # given gamma parameter values
  LM <- NULL            # data frame built by numerical integration
  for(kk in KK) {
    lmr  <- theoLmoms(list(para=c(MU, SD, kk), type="gno"), nmom=4,
                      subdivisions=200, rel.tol=.Machine$double.eps^0.50)
    LM   <- rbind(LM, data.frame(xi=0, alpha=1, kappa=kk,
                                 L1=lmr$lambdas[1], L2=lmr$lambdas[2],
                                 T3=lmr$ratios[ 3], T4=lmr$ratios[ 4]))
  }
  plot(log10(LM$kappa), LM$T3, type="l")


LM <- LM[abs(LM$T3) < 1,]
LM <- LM[abs(LM$T4) < 1,]
LM <- LM[! LM$T4 < (5 * LM$T3^2 - 1)/4,]

wnt <- diff(c(9999, LM$T3)) > 0
while(1) {
  plot(LM$kappa, LM$T3, type="l")
  wnt <- diff(c(9999, LM$T3)) > 0
  points(LM$kappa[wnt], LM$T3[wnt], pch=16, col="red")
  LM <- LM[! wnt,]
  if(any(wnt)) next
  break
}
wnt <- diff(c(9999, LM$T3)) > 0
plot(  LM$kappa,      LM$T3, type="l")
points(LM$kappa[wnt], LM$T3[wnt], pch=16, col="red")



wnt <- diff(c(9999, LM$L1)) > 0
plot(  -LM$L1, LM$kappa, log="xy")
points(-LM$L1[wnt], LM$kappa[wnt], pch=16, col="red")
while(1) {
  plot(-LM$L1, LM$kappa, type="l", log="xy")
  wnt <- diff(c(9999, LM$L1)) > 0
  points(-LM$L1[wnt], LM$kappa[wnt], pch=16, col="red")
  LM <- LM[! wnt,]
  if(any(wnt)) next
  break
}
wnt <- diff(c(9999, LM$L1)) > 0
plot(  -LM$L1, LM$kappa, type="l", log="xy")
points(-LM$L1[wnt], LM$kappa[wnt], pch=16, col="red")

row.names(LM) <- NULL
LMo <- LM
#stop()

plot(  -LM$L1, LM$kappa, type="p", log="xy")
plot(   LM$L2, LM$kappa, type="p", log="xy")
plot(   LM$T3, LM$kappa, type="p", log="y")
plot(   LM$T4, LM$kappa, type="p", log="y")
LM$diff_kappa <- c(NA, diff(LM$kappa))
row.names(LM) <- NULL

LM <- LM[-c(3627),]
plot(  -LM$L1, LM$kappa, type="p", log="xy")
LM$diff_kappa <- c(NA, diff(LM$kappa))
row.names(LM) <- NULL
LM <- LM[-c(3567),]
plot(  -LM$L1, LM$kappa, type="p", log="xy")
LM$diff_kappa <- c(NA, diff(LM$kappa))
row.names(LM) <- NULL
LM$diff_kappa <- NULL

  # from Hosking's FORTRAN code in 1990s (C[1-3]; D[1-6])
    A0 <-  0.20466534e1;   A1 <- -0.36544371e+1;
    A2 <-  0.18396733e+1;  A3 <- -0.20360244;
    B1 <- -0.20182173e+1;  B2 <-  0.12420401e+1;  B3 <- -0.21741801
  AoesA <- c(0.20466534e1, -0.36544371e+1,  0.18396733e+1, -0.20360244)
  BoesA <- c(-0.20182173e+1 , 0.12420401e+1, -0.21741801)

  Apara.init <- AoesA
  Bpara.init <- BoesA
  para.init  <- c(Apara.init, Bpara.init)

  "tau3kappa" <-
  function(tau3, A=Apara.init, B=Bpara.init,
           SMALL=.Machine$double.eps) {
     sg <- -sign(  tau3 )
     kk <- sapply(tau3, function(t3) {
         T3 <- -abs(t3)
         #if(T3 <= SMALL) return(0)
         TT <- T3*T3
         K <- -T3*(A[1]+TT*(A[2]+TT*(A[3]+TT*A[4])))/(1+TT*(B[1]+TT*(B[2]+TT*B[3])))
       return(K) })
    return(sg * kk)
  }

  "ofunc" <- function(para, t3=NA, tkappa=NA) {
    gg <- tau3kappa(t3, A=para[1:4], B=para[5:7])
    return(sum((gg - tkappa)^2)) # sum of squared error
  }

  rt <- optim(para.init, fn=ofunc, t3=LM$T3, tkappa=LM$kappa,
              control=list(reltol=.Machine$double.eps^0.5, maxit=2000))
  AoesB <- rt$par[1:4]; BoesB <- rt$par[5:7]

  KappaHosk <- LM$kappa_hosking <- tau3kappa(LM$T3)
  KappaMade <- LM$kappa_lmomco  <- tau3kappa(LM$T3, A=AoesB, B=BoesB)

  plot(LM$kappa, KappaMade, lwd=0.7, pch=16, cex=0.5, log="",
       xlab="Given GNO kappa", ylab="Kappa by new coefficients")
  abline(0,1, col="red")

  plot(LM$T3, abs(KappaHosk-LM$kappa), lwd=0.7, pch=16, cex=0.5, col="salmon",
       xlab="Integrated T3 given kappa", ylab="abs(Kappa difference)")
  points(LM$T3,  abs(KappaMade-LM$kappa), lwd=0.7, pch=16, cex=0.5, col="orchid")

  message("Hosking Aoes: ", paste(AoesA, collapse=", "))
  message("optim() Aoes: ", paste(AoesB, collapse=", "))
  message("Hosking Boes: ", paste(BoesA, collapse=", "))
  message("optim() Boes: ", paste(BoesB, collapse=", "))

  sqrt(mean((KappaHosk-LM$kappa)^2)) # [1] 0.02695981
  sqrt(mean((KappaMade-LM$kappa)^2)) # [1] 0.001721387
  100 * (0.001291257 - 0.02695981)/abs(0.02695981)

  plot(LM$kappa, LM$kappa_lmomco)
  abline(0,1, col="red")

  summary(LM)
# We see in summary(LM) that T3 min is -0.9990979, so let us truncate pargno() at |T3| == 0.999
# instead of Hosking's 0.95.

