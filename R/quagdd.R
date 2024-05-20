"quagdd" <- function(f, para, paracheck=TRUE, ...) {
  if(! check.fs(f)) return()
  if(paracheck) if(! are.pargdd.valid(para)) return()

  A1 <- para$para[1];  B1 <- para$para[2]
  A2 <- para$para[3];  B2 <- para$para[4]
  if(length(para$para) == 5) {
    if(para$para[5] == 1) {
      A2 <- A1
      B2 <- B1
    }
  }

  eps <- .Machine$double.eps^0.5

  f[f <     eps] <-     eps
  f[f > 1 - eps] <- 1 - eps

  "afunc" <- function(x, Fx=NA) {
         theF <- cdfgdd(x, para, paracheck=FALSE)
         #message("theX=", x,"  theF=", theF); message("Fx=", Fx)
         err <- Fx - theF; #message(err)
         return(err)
  }

  mugdd <-       A1/B1 - A2/B2
  sdgdd <- sqrt( A1/B1 + A2/B2 )
  skgdd <- 2*(A1*B2^3 - A2*B1^3) / (A1*B2^2 + A2*B1^2)^(3/2)

  #mx <- mugdd
  #i <- 0; maxit <- 100
  #while(1) {
  #  f1 <- cdfgdd(mx, para)
  #  if(f1 < 1 - eps) {
  #    i <- i + 1
  #    mx <- (mx + sdgdd*i/4); print(mx); print(f1, 16)
  #    if(i > maxit) break
  #    next
  #  }
  #  break
  #}
  #mn <- sdgdd
  #i <- 0
  #while(1) {
  #  f0 <- cdfgdd(mn, para)
  #  if(f0 > eps) {
  #    i <- i + 1
  #    mn <- (mn - sdgdd*i/4); print(mn); print(f0, 16)
  #    if(i > maxit) break
  #    next
  #  }
  #  break
  #}
  #lower <- mn
  #upper <- mx
  lower <- -1000
  upper <- +1000
  x <- sapply(seq_len(length(f)), function(i) {
                 rt <- NULL
                 try( rt <- uniroot(afunc, c(lower, upper), Fx=f[i], ...), silent=TRUE)
                 #print(rt)
                 ifelse(is.null(rt), return(NA), return(rt$root)) })
  return(x)
}
