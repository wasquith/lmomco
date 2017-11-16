"pars2x" <- function(f, paras, na.rm=FALSE, ...) {
   a <- b <- c <- rep(NA, length(f))
   if(! is.null(paras$lmr)) a <- qlmomco(f, paras$lmr)
   if(! is.null(paras$mle)) b <- qlmomco(f, paras$mle)
   if(! is.null(paras$mps)) c <- qlmomco(f, paras$mps)
   zzz <- data.frame(lmr=a, mle=b, mps=c)
   if(length(a[! is.na(a)]) == 0 &
      length(b[! is.na(b)]) == 0 &
      length(c[! is.na(c)]) == 0) {
      warning("all the parameters methods have NULL")
      return(NULL)
   }
   zz <- t(zzz); n <- 1:length(zz[1,])
   mu <- sapply(n, function(i) {       mean(unique(zz[,i]), na.rm=na.rm)  })
   sd <- sapply(n, function(i) {         sd(unique(zz[,i]), na.rm=na.rm)  })
   nn <- sapply(n, function(i) { g<-zz[,i]; length(unique(g[! is.na(g)])) })
   zzz$f <- f; zzz$lwr <- mu - sd; zzz$mean <- mu; zzz$upr <- mu + sd
   zzz$sd <- sd; zzz$n <- nn; all <- c(zzz$lmr, zzz$mle, zzz$mps)
   tmpA <- c(      f[1],       f, f[length(f)], rev(f))
   tmpB <- c(zzz$upr[1], zzz$lwr, zzz$upr[length(f)], rev(zzz$upr))
   zzz <- zzz[,c(4,1:3,5:9)]; all <- all[! is.na(all)]
   attr(zzz, "f.poly") <- tmpA; attr(zzz, "x.poly") <- tmpB
   attr(zzz, "all.summary") <- summary(all, na.rm=TRUE)
   return(zzz)
}
