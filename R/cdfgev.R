"cdfgev" <-
function(x, para, paracheck=TRUE) {
    if(paracheck) {
      if(! are.pargev.valid(para)) return()
    }
    U <- para$para[1]
    A <- para$para[2]
    K <- para$para[3]

    Y <- (x - U)/A
    ZERO <- sqrt(.Machine$double.eps)
    if(abs(K) > ZERO) {
      Y <- suppressWarnings( -log(1-K*Y)/K )
    }
    f <- exp(-exp(-Y))

    if(K < 0) {
       f[!is.finite(f)] <- 0
    } else if(K > 0) {
       f[!is.finite(f)] <- 1
    }
    names(f) <- NULL
    return(f)
}

