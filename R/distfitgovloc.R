"disfitgovloc" <-
function(x=NULL, loc=NULL, lwr=0, upr=NA, init.para=NULL,
         loctype=c("mean", "median"), objfun=c("rmse", "mad"),
         ptransf=function(p) return(log(p)),
         pretransf=function(p) return(exp(p)),
         silent=TRUE, verbose=FALSE, ...) {

  objfun  <- match.arg(objfun )
  loctype <- match.arg(loctype)

  zz <- list(type="gov", para=rep(NA, 3), source="disfitgovloc",
              supdist=NULL, init.para=init.para, optim=NULL, message="")

  if(is.null(lwr) | is.na(lwr)) {
    zz$message <- "need lower bounds as the lwr argument"
    return(zz)
  }
  if(is.null(upr) | is.na(upr)) {
    zz$message <- "need upper bounds as the upr argument"
    return(zz)
  }
  if(upr <= lwr) {
    zz$message <- "must have the upr argument as greater than the lwr argument"
    return(zz)
  }

  OBJfun <- function(hat, obs) NULL
  if(objfun == "rmse") {
    OBJfun <- function(hat, obs) sqrt(sum(   (obs - hat)^2) / length(obs) )
  } else if(objfun == "mad") {
    OBJfun <- function(hat, obs)      sum(abs(obs - hat))   / length(obs)
  } else {
    stop("should not be here in logic")
  }

  fn <- function(par, lwr, loc, upr) {
    tgov <- vec2par(pretransf(par), type="gov", paracheck=FALSE)
    ploc <- ifelse(loctype == "mean", lmomgov(tgov)$lambdas[1], quagov(0.5, tgov))
    sup  <- quagov(c(0,1), tgov)
    if(verbose) print(c(sup[1], ploc, sup[2], lwr, loc, upr))
    OBJfun(c(sup[1], ploc, sup[2]), c(lwr, loc, upr))
  }

  if(is.null(init.para) & ! is.null(x)) {
    lmr <- lmoms(x, nmom=3, no.stop=TRUE)
    init.para <- pargov(lmr)
    loc <- ifelse(loctype == "mean", lmr$lambdas[1], median(x))
  } else if(is.null(init.para) & ! is.null(loc)) {
    lmr <- lmomtri(vec2par(c(lwr, (lwr+upr)/2, upr), type="tri"))
    init.para <- pargov(lmr)
  } else if(is.null(init.para)) {
    zz$message <- "must have argument x or argument loc or argument init.para set"
    return(zz)
  }

  if(loc < lwr) {
    zz$message <- "location loc is less than given lower bounds, no solution possible"
    return(zz)
  }
  if(loc > upr) {
    zz$message <- "location loc is greater than given upper bounds, no solution possible"
    return(zz)
  }

  rt <- NULL
  try( rt <- optim(ptransf(init.para$para), fn=fn, lwr=lwr, loc=loc, upr=upr, ...), silent=silent )
  if(is.null(rt)) {
    zz$message <- ("optim() attempt is NULL")
    return(zz)
  } else {
    gov <- vec2par(pretransf(rt$par), type="gov")
    gov$source    <- "disfitgovloc"
    gov$supdist   <- quagov(c(0,1), gov)
    gov$init.para <- init.para
    gov$optim     <- rt
    gov$message   <- ""
    if(rt$convergence != 0) gov$message <- "optim() reports nonzero convergence"

    return(gov)
  }
}
