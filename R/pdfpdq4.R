"pdfpdq4" <-
function(x, para, paracheck=TRUE, h=NA, hfactor=0.2) {
  # Check that the parameters are valid one time
  # then use the paracheck switch on quapdq4 for
  # an extreme speed up on this algorithm.
  if(paracheck == TRUE) {
    if(! are.parpdq4.valid(para)) return()
  }
  if(is.na(h)) h <- para$para[2] * hfactor
  # [-f(x+2h)+8f(x+h)-8f(x-h)+f(x-2h)] / 12h
  f1 <-   -cdfpdq4(x+2*h, para, paracheck=FALSE)
  f2 <-  8*cdfpdq4(x + h, para, paracheck=FALSE)
  f3 <- -8*cdfpdq4(x - h, para, paracheck=FALSE)
  f4 <-    cdfpdq4(x-2*h, para, paracheck=FALSE)
  return( (f1+f2+f3+f4) / (12*h) )
  # (f(x+h) - f(x-h)) / 2h
  #f1 <- cdfpdq4(x+h, para, paracheck=FALSE)
  #f2 <- cdfpdq4(x-h, para, paracheck=FALSE)
  #return( (f1 - f2) / (2*h) )
}
