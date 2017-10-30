"mle2par" <-
function(x, type, para.int=NULL, silent=TRUE, null.on.not.converge=TRUE,
                  ptransf=  function(t) return(t),
                  pretransf=function(t) return(t), ...) {

  if(is.null(para.int)) {
     lmr <- lmoms(x)
     if(! are.lmom.valid(lmr)) {
        warning("L-moments of x are not valid for initial parameters, ",
                "try manual initial parameters")
        return(NULL)
     }
     para.int <- lmom2par(lmr, type=type, ...)
     if(is.null(para.int)) {
        warning("L-moments of x are not valid for initial parameters, ",
                "try manual initial parameters")
        return(NULL)
     }
  } else if(is.vector(para.int)) {
     para.int <- vec2par(para.int, type=type)
     if(is.null(para.int)) {
        warning("initial parameters given by vector are not valid for initial parameters, ",
                "try other initial parameters")
        return(NULL)
     }
  }

  if(length(para.int$para) == 1) {
     warning("function is not yet built for single parameter optimization")
     return(NULL)
  }

  "afunc" <- function(para, x=NULL, ...) {
       #print(para)
       lmomco.para <- vec2par(pretransf(para), type=type, paracheck=TRUE)
       if(is.null(lmomco.para)) return(Inf)
       #print(lmomco.para$para)
       pdf <- par2pdf(x,lmomco.para) # pull into local scope, in case of later
       # interception of problems
       # The negative is to accommodate the minimization setup of optim()
       L <- -sum(log(pdf), na.rm=TRUE) # lmomco should fill NAs with zeros by
       # global package design assumptions but just incase some leaked through na.rm=T
       if(! silent) message(" L=",L)
       return(L)
  }

  #   print(ptransf(para.int$para))
  #   raw.afunc.call <- afunc(ptransf(para.int$para), x=x)
  #   print("RAW afunc() CALL WITH INITIAL PARAMETERS")
  #   print(raw.afunc.call)

  # Note for some reason there is an argument name clash(?) that WHA can not get to
  # the bottom of for a case where ... in the main call has say p=3 to trigger the
  # generalized gamma if type="gam" so ... is not wrapped into the optim() call as
  # WHA would instinctively do.

  rt <- NULL
  try(rt <- optim(par=ptransf(para.int$para), fn=afunc, x=x, ...), silent=silent)
  if(is.null(rt)) {
     warning("optim() attempt is NULL")
     return(NULL)
  } else {
     if(null.on.not.converge & rt$convergence != 0) {
        warning("optim() reports convergence error")
        return(NULL)
     }
     lmomco.para <- vec2par(pretransf(rt$par), type=type)
     lmomco.para$AIC <- 2*length(rt$par) - 2*(-1*rt$value) # Note the sign change
     # because the optimize is working in the opposite direciton by default.
     lmomco.para$optim <- rt
     return(lmomco.para)
  }
}
