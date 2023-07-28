lmrloc <- function(x, y=NULL, terse=TRUE) {
  if(is.null(y)) {
    if(ncol(x) != 2) {
      warning("y is NULL, but x does not have exactly two columns")
      return(NULL)
    }
  } else {
    if(length(x) != length(y)) {
      warning("x and y must have same length")
      return(NULL)
    }
    x <- data.frame(x=x, y=y)
  }
  x <- x[complete.cases(x), ]
  y <- x[,2]; x <- x[,1]

  r <- sign( stats::cor(x,y, method="spearman") )
  names(r) <- "Sign Spearman Rho"

  n <- length(x)

  mu_x <- mean(x); names(mu_x) <- "Mean X"
  mu_y <- mean(y); names(mu_y) <- "Mean Y"

  gini_x <- ( 2 / (n * (n-1) ) ) * sum(sort(x) * seq( (1-n), (n-1), by=2) )
  gini_y <- ( 2 / (n * (n-1) ) ) * sum(sort(y) * seq( (1-n), (n-1), by=2) )
  names(gini_x) <- "Gini X"
  names(gini_y) <- "Gini Y"

  lmr_m   <- r * ( gini_y / gini_x )
  lmr_b   <-  mu_y - (lmr_m * mu_x)
  loc_lmr <- c(lmr_b, lmr_m)
  names(loc_lmr) <- c("LMR_Intercept", "LMR_Slope")

  sd_x <- stats::sd(x); names(sd_x) <- "Stdev X"
  sd_y <- stats::sd(y); names(sd_y) <- "Stdev Y"

  pmr_m   <- r * ( sd_y / sd_x )
  pmr_b   <-  mu_y - (pmr_m * mu_x)
  loc_pmr <- c(pmr_b, pmr_m)
  names(loc_pmr) <- c("PMR_Intercept", "PMR_Slope")

  if(terse) {
    return( list(loc_lmr=loc_lmr, loc_pmr=loc_pmr) )
  } else {
    zz <- list(loc_lmr=loc_lmr, loc_pmr=loc_pmr, srho=r,
               mu=c(mu_x, mu_y), gini=c(gini_x, gini_y), sd=c(sd_x, sd_y))
    return(zz)
  }
}
