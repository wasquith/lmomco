"par2qua2" <-
function(f, para1, para2, wfunc=NULL, weight=NULL, as.list=FALSE, ...) {
  Q1 <- par2qua(f, para1, ...)
  Q2 <- par2qua(f, para2, ...)

  if(is.null(wfunc)) {
     if(is.null(weight)) {
        Q <- (1-f)*Q1 + f*Q2
        zz <- data.frame(f=f, qua=Q); tmp <- zz; tmp <- tmp[order(tmp$f),]
        if(any(diff(tmp$qua) < 0)) warning("result is nonmonotonic increasing")
        if(as.list) return(zz)
        return(Q)
     }

     if(length(weight) == 1)  weight <- c((1-weight), weight)

     if(length(weight) == 2) {
        if(sum(weight) != 1) {
           warning("sum of the two weights is not unity, going to rescale as such")
           weight <- weight/sum(weight)
        }
     } else {
        warning("weight can not be a vector longer than 2")
        return()
     }
     Q <- weight[1]*Q1 + weight[2]*Q2
     zz <- data.frame(f=f, qua=Q); tmp <- zz; tmp <- tmp[order(tmp$f),]
     if(any(diff(tmp$qua) < 0)) warning("result is nonmonotonic increasing")
     if(as.list) return(zz)
     return(Q)
  } else {
     if(! is.function(wfunc)) {
        warning("wfunc is not a function")
        return(NA)
     }
     weights <- wfunc(f)
     Q <- Q1*(1-weights) + Q2*weights
     zz <- data.frame(f=f, qua=Q); tmp <- zz; tmp <- tmp[order(tmp$f),]
     if(any(diff(tmp$qua) < 0)) warning("result is nonmonotonic increasing")
     if(as.list) return(zz)
     return(Q)
  }
}
