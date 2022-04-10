"is.pdq4" <- function(para) {
  if(para$type != "pdq4") {
    warning("Parameters are not pdq4 parameters")
    return(FALSE)
  }
  return(TRUE)
}
