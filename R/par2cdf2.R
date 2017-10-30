"par2cdf2" <-
function(x, para1, para2, weight=NULL, ...) {
    end.min <-   .Machine$double.eps
    end.max <- 1-.Machine$double.eps
    qua.min <- par2qua2(end.min, para1, para2, ...)
    qua.max <- par2qua2(end.max, para1, para2, ...)

    "fn" <- function(ff, quax=NULL) {
       return(quax - par2qua2(ff, para1, para2, weight=weight, ...))
    }

    f <- sapply(1:length(x), function(i) {
                 QUAx <- x[i]
                 if(QUAx <= qua.min) return(end.min)
                 if(QUAx >= qua.max) return(end.max)
                 rt <- NULL
                 try(rt <- uniroot(fn, c(end.min, end.max), quax=QUAx), silent=TRUE)
                 if(is.null(rt)) {
                    warning("uniroot failed for QUAx=", QUAx)
                    return(NA)
                 }
                 return(rt$root) })
    return(f)
}
