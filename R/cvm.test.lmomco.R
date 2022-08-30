"cvm.test.lmomco" <- function(x, para1, ...) {
   x <- sort(x); n <- length(x)
   # next two lines are core mimic to goftest::cvm.test()
   omega.sq <- 1/(12*n) + sum(((2*(1:n)-1)/(2*n) - par2cdf(x, para1, ...))^2)
   pval <- goftest::pCvM(omega.sq, n=n, lower.tail=FALSE)
   names(omega.sq) <- "Omega-Squared"; names(pval) <- "p-value"
   return(list(null.dist=para1, text="Cramer--von Mises test of goodness-of-fit",
               statistic=omega.sq, p.value=pval, source="cvm.test.lmomco"))
}
