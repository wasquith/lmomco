"cdfglo" <-
function(x,para) {
    if(! are.parglo.valid(para)) return()
    XI <- para$para[1]
    A  <- para$para[2]
    K  <- para$para[3]

    f <- vector(mode="numeric", length=length(x))
    Y <- (x-XI)/A
    if(K == 0) {
       f <- 1/(1+exp(-Y))
    } else {
       ARG <- 1-K*Y
       Y <- suppressWarnings( -log(ARG)/K )
       f <- 1/(1+exp(-Y))
    }
    if(K < 0) {
       f[!is.finite(f)] <- 0
    } else if(K > 0) {
       f[!is.finite(f)] <- 1
    }
    names(f) <- NULL
    return(f)
}

