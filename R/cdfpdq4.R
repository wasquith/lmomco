"cdfpdq4" <-  function(x, para, paracheck=TRUE) {
  # Check that the parameters are valid one time
  # then use the paracheck switch on quapdq4 for
  # an extreme speed up on this algorithm.
  if(paracheck == TRUE) {
    if(! are.parpdq4.valid(para)) return()
  }

  f <- sapply(1:length(x), function(i) {
            QUAx <- x[i]
            if(QUAx == -Inf) return(0)
            if(QUAx == +Inf) return(1)
            fn <- function(of, target_qua=NA) {
                      qua <- quapdq4(of, para, paracheck=FALSE)
                      val <- target_qua - qua
               return(val)
            }
            root <- uniroot(fn, c(0,1), target_qua=x[i])
            return(root$root)
  })
  names(f) <- NULL
  return(f)
}
