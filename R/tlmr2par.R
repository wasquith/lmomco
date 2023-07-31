"tlmr2par" <-
function(x, type, init.para=NULL,
                  trim=NULL, leftrim=NULL, rightrim=NULL, ...) {

   if(is.null(type)) {
     warning("must specify distribution type")
     return(NULL)
   }
   if(is.null(init.para)) {
     lmr <- lmoms(x)
     if(!are.lmom.valid(lmr)) {
       warning("L-moments of x are not valid for initial parameters, ",
               "try manual initial parameters")
       return(NULL)
     }
     init.para <- lmom2par(lmr, type=type, ...)
     if(is.null(init.para)) {
       warning("could not estimate initial parameters via L-moments")
       return(NULL)
     }
   } else if(!is.list(init.para) & is.vector(init.para)) {
     init.para <- vec2par(init.para, type = type)
     if(is.null(init.para)) {
       warning("initial parameters given by vector are not valid for initial parameters, ",
               "try other initial parameters")
       return(NULL)
     }
   }

   if(is.null(init.para)) {
     warning(" initial parameters are NULL")
     return(NULL)
   }
   if(init.para$type != type) {
     warning("distribution requested to fit does not match the type of the ",
             "initial parameters")
     return(NULL)
   }
   if(length(init.para$para) == 1) {
     warning("function is not yet built for single parameter optimization")
     return(NULL)
   }

  nmom <- length(init.para$para)

  if(  is.null(trim) & is.null(leftrim) & is.null(rightrim)) trim <- 0
  if(! is.null(trim)) { # this is setting symmetrical trimming as a shortcut
    leftrim  <- trim
    rightrim <- trim
  }
  if(length(unique(x)) == 1) { # L-moments break down if all values are equal
    warning("all values are equal--TLmoments can not be computed")
    return(NULL)
  }

  # compute the sample TL-moments as desired from the sample, this trimming does
  # not throwaway the tails but can be thought of as "deleveraging" the tails
  tlmr <- TLmoms(x, leftrim=leftrim, rightrim=rightrim, nmom=nmom)
  # we set the nmom here only for (1) speed efficiency, don't compute deeper than
  # needed and to (2) avoid a bug of vector recycling by R in sum of squares later

  # the heart of optimization is some type of object function, here we quietly
  # define one on the fly. The function's first argument is a simple vector
  # controlled a level about by the optim() function to come. The rest of the
  # arguments are named arguments to be set when the optim() function is used.
  "afunc" <- function(para, tlmr=NULL, type=type,
                            leftrim=leftrim, rightrim=rightrim) {
     new.para <- vec2par(para, type=type, paracheck=FALSE) # convert to lmomco parameter list
     # there are likely no (or we don't even bother to look for or derive) expressions of
     # distribution parameters in terms of TL-moments, so optimization is the route
     # we have a "guess" at the parameters in new.para and now we know how to convert
     # those to TL-moments
     fit.tlmr <- theoTLmoms(new.para, leftrim=leftrim, rightrim=rightrim, nmom=nmom)
     if(fit.tlmr$lambdas[2] <= 0) return(Inf)
     err <- sum((tlmr$lambdas - fit.tlmr$lambdas)^2) # sum of squares
     return(err) # Sum of square errors, we want to minimize this quantity!!!
  }

  rt <- NULL # standard hack around optim() is to try()
  try(rt <- stats::optim(init.para$para,    afunc, tlmr=tlmr, type=type,
                                 leftrim=leftrim, rightrim=rightrim, ...))
  if(is.null(rt)) {
    warning("failure, so returning NULL, insert further advice to the user")
    return(NULL)
  }
  trim.para <- vec2par(rt$par, type=type) # final the formal lmomco parameter list
  if(is.null(trim.para)) {
    trim.para <- list(para=rep(NA, nmom), text="invalid parameters, see rt, try a different init.para")
  }
  trim.para$source <- "tlmr2par"
  trim.para$rt <- rt # store the results for later use by the user if ever needed
  trim.para$init.para <- init.para
  return(trim.para)
}
