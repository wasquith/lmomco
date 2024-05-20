"is.gdd" <-
function(para) {
    if(para$type != "gdd") {
      warning("Parameters are not Gamma Difference parameters")
      return(FALSE)
    }
    return(TRUE)
}

