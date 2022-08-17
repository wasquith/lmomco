"tlmrpdq4" <-
function(trim=NULL, leftrim=NULL, rightrim=NULL, xi=0, alpha=1,
         kbeg=-5, kend=0, by=.001) {
   ks <- -1*10^(seq(4,-2, by=-.1))
  ks <- c(ks, 10^seq(-5,0, by=.1))
  #ks <- seq(kbeg, kend, by=by)
  n <- length(ks)
  T2 <- T3 <- T4 <- T5 <- T6 <- vector(mode="numeric", length=n)
  i <- 0
  for(k in ks) {
    tmp.para <- vec2par(c(xi, alpha, k), type="pdq4", paracheck=FALSE)
    tmp.lmr  <- theoTLmoms(tmp.para, nmom=6,
                  trim=trim, leftrim=leftrim, rightrim=rightrim)
    i <- i + 1
    T2[i] <- tmp.lmr$ratios[2]
    T3[i] <- tmp.lmr$ratios[3]
    T4[i] <- tmp.lmr$ratios[4]
    T5[i] <- tmp.lmr$ratios[5]
    T6[i] <- tmp.lmr$ratios[6]
  }
  z <- list(ks=ks, tau2=T2, tau3=T3, tau4=T4, tau5=T5, tau6=T6)
  return(z)
}

