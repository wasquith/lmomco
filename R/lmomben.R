"lmomben" <-
function(para=list(para=c(1, 10)), ...) {
  m <- para$para[1]
  b <- para$para[2]
  if(b != 10) {
    warning("lmomben is only configured for base10 at the moment")
    return(NULL)
  }
  z <- list(L1 = NULL, L2  = NULL, TAU3 = NULL, TAU4 = NULL,
                       LCV = NULL, L3   = NULL, L4   = NULL,
            source = "lmomben")
  attributes(para$para) <- NULL
  if(m == 1) {
    L1 <- 3.43908699617500524
    L2 <- 1.34518434179517077
    T3 <- 0.24794090889493661
    T4 <- 0.01614509742647182
  } else if(m == 2) {
    L1 <- 38.59062918136093145
    L2 <- 13.81767809210059283
    T3 <-  0.22237541787527126
    T4 <-  0.03541037418894027
  } else if(m == 3) {
    L1 <- 390.36783537821605705
    L2 <- 138.21917489739223583
    T3 <-   0.22192482374529940
    T4 <-   0.03571514686148788
  } else {
    warning('specified number of first significant digits is not supported')
    return(NULL)
  }

  z$LCV  <- L2 / L1
  z$L3   <- T3 * L2
  z$L4   <- T4 * L2
  z$L1 <- L1; z$L2 <- L2; z$TAU3 <- T3; z$TAU4 <- T4
  z <- lmomco::lmorph(z)
  return(z)
}
