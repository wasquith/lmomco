"is.pdq3" <- function(para) {
  if(para$type != "pdq3") {
    warning("Parameters are not pdq3 parameters")
    return(FALSE)
  }
  return(TRUE)
}
