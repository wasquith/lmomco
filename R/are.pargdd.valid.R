"are.pargdd.valid" <-
function(para, nowarn=FALSE) {
  if(! is.gdd(para)) return(FALSE)

  A1 <- para$para[1];  B1 <- para$para[2]
  A2 <- para$para[3];  B2 <- para$para[4]
  if(length(para$para) == 5) {
    if(! is.na(para$para[5]) & para$para[5] == 1) {
      A2 <- A1
      B2 <- B1
    }
  }
  #print(para$para)
  op <- options()
  GO <- TRUE
  if(nowarn == TRUE) options(warn=-1)
  if(is.na(A1)) return(FALSE)
  if(A1 <= 0) {
    warning("Parameter A1 is not > 0, invalid")
    GO <- FALSE
  }
  if(B1 <= 0) {
     warning("Parameter B1 is not > 0, invalid")
     GO <- FALSE
  }
  if(A2 <= 0) {
     warning("Parameter A2 is not > 0, invalid")
     GO <- FALSE
  }
  if(B2 <= 0) {
     warning("Parameter B2 is not > 0, invalid")
     GO <- FALSE
  }
  return(TRUE)
  options(op)
  if(GO) return(TRUE)
  return(FALSE)
}

