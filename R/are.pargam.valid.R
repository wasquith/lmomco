"are.pargam.valid" <-
function(para,nowarn=FALSE) {
    if(! is.gam(para)) return(FALSE)
    if(any(is.na(para$para))) return(FALSE)

    op <- options()
    if(nowarn == TRUE) options(warn=-1)
    GO <- TRUE
    if(length(para$para) == 2) {
       ALPHA <- para$para[1]
       BETA  <- para$para[2]
       if(ALPHA <= 0) {
         warning("Parameter ALPHA is not > 0, invalid")
         GO <- FALSE
       }
       if(BETA <= 0) {
         warning("Parameter BETA is not > 0, invalid")
         GO <- FALSE
       }
    } else if(length(para$para) == 3) {
       MU <- para$para[1]; SIGMA <- para$para[2]; NU <- para$para[3]
       if(MU <= 0) {
         warning("Parameter MU is not > 0, invalid")
         GO <- FALSE
       }
       if(SIGMA <= 0) {
         warning("Parameter SIGMA is not > 0, invalid")
         GO <- FALSE
       }
       # -Inf < NU < Inf
    } else {
       stop("should not be here in logic flow")
    }
    options(op)
    if(GO) return(TRUE)
    return(FALSE)
}

