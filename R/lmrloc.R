lmrloc <- function(x, y=NULL) {
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
  y <- x[,2]
  x <- x[,1]

  r <- sign( stats::cor(x,y, method="spearman") )

  n <- length(x)

  x <- sort(x)
  y <- sort(y)

  gini_x <- ( 2 / (n * (n-1) ) ) * sum(x * seq( (1-n), (n-1), by=2) )
  gini_y <- ( 2 / (n * (n-1) ) ) * sum(y * seq( (1-n), (n-1), by=2) )
  m <- r * gini_y / gini_x
  b <- mean(y) - (m * mean(x))
  loc_lmr <- c(b, m)
  names(loc_lmr) <- c("Intercept", "Slope")

  sd_x <- stats::sd(x)
  sd_y <- stats::sd(y)
  m <- r * sd_y / sd_x
  b <- mean(y) - (m * mean(x))
  loc_pmr <- c(b, m)
  names(loc_pmr) <- c("Intercept", "Slope")

  list(loc_lmr=loc_lmr, loc_pmr=loc_pmr)
}

# x <- rnorm(500)
# y <- -0.4 * x + rnorm(500, sd=0.2)
# y[1] <- 0
# lmrloc(x, y)
# plot(x,y)
# abline(lmrloc(x, y)[[1]], lty=1)
# abline(lmrloc(x, y)[[2]], lty=2)

