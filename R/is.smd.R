"is.smd" <-
function(para) {
    if(para$type != "smd") {
      warning("Parameters are not Singh Maddala parameters")
      return(FALSE)
    }
    return(TRUE)
}

