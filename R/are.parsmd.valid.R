"are.parsmd.valid" <-
function(para,nowarn=FALSE) {
    if(! is.smd(para)) return(FALSE)
    if(any(is.na(para$para))) return(FALSE)

    U <- para$para[1] # location U >= 0
    A <- para$para[2] # scale A > 0
    B <- para$para[3] # shape B > 0
    Q <- para$para[4] # shape Q > 0
    # 0 <= x <= +Inf

    op <- options()
    GO <- TRUE
    if(nowarn == TRUE) options(warn=-1)
    if(A <= 0) {
      warning("Parameter A is not > 0, invalid")
      GO <- FALSE
    }
    if(B <= 0) {
      warning("Parameter B is not > 0, invalid")
      GO <- FALSE
    }
    if(Q <= 0) {
      warning("Parameter Q is not > 0, invalid")
      GO <- FALSE
    }
    options(op)
    if(GO) return(TRUE)
    return(FALSE)
}

