\encoding{utf8}
\name{lmompdq3}
\alias{lmompdq3}
\title{L-moments of the Polynomial Density-Quantile3 Distribution}
\description{
This function estimates the L-moments of the Polynomial Density-Quantile3 distribution given the parameters (\eqn{\xi}, \eqn{\alpha}, and \eqn{\kappa}) from \code{\link{parpdq3}}. The L-moments in terms of the parameters are
\deqn{\lambda_1 = \xi + \alpha\bigl[(1+\kappa)\log(1+\kappa) - (1-\kappa)\log(1-\kappa) - \kappa\log(4)\bigr]\mbox{,}}
\deqn{\lambda_2 = \frac{\alpha(1-\kappa^2)}{(1-\kappa\tau_3)}\mbox{,}}
\deqn{\tau_3 = \frac{1}{\kappa} - \frac{1}{\mathrm{arctanh}(\kappa)} \mbox{, and}}
\deqn{\tau_4 = (5\tau_3/\kappa) - 1\mbox{.}}
}
\usage{
lmompdq3(para, paracheck=TRUE)
}
\arguments{
  \item{para}{The parameters of the distribution.}
  \item{paracheck}{A logical switch as to whether the validity of the parameters should be checked. Default is \code{paracheck=TRUE}.}
}
\value{
  An \R \code{list} is returned.

  \item{lambdas}{Vector of the L-moments. First element is
\eqn{\lambda_1}, second element is \eqn{\lambda_2}, and so on.}
  \item{ratios}{Vector of the L-moment ratios. Second element is
\eqn{\tau}, third element is \eqn{\tau_3} and so on.}
  \item{trim}{Level of symmetrical trimming used in the computation, which is \code{0}.}
  \item{leftrim}{Level of left-tail trimming used in the computation, which is \code{NULL}.}
  \item{rightrim}{Level of right-tail trimming used in the computation, which is \code{NULL}.}
  \item{source}{An attribute identifying the computational
source of the L-moments: \dQuote{lmompdq3}.}
}
\note{
Polynomial approximations for the \eqn{\tau_3} and \eqn{\tau_4} are developed here. First, the author's monograph (Asquith, 2011, table 10.1) shows five digits for such approximates for other distributions, so the code below will used the core basis, five digits. Second, an approximation means that \code{\link{lmrdia}} does not have the internal burden of using \code{uniroot()} to solve for the coordinates for the L-moment ratio diagram. The following code represents an exploration towards the definition of a helper function, \code{t4pdq3()}, which is repeated inside the internals of \code{\link{lmrdia}} in order to support the PDQ3. The trajectory of the PDQ3 resides at or above that for the generalized logistic distribution (\code{\link{quaglo}}) that is well known to L-moment theory. In conclusion, the 5-digit approximation provides a maximum absolute \eqn{\tau_4} error of about 0.00055.
\preformatted{
  fn <- function(k, tau3=NA) { t3 <- (1/k - 1/atanh(k))
                               if(is.nan(t3)) t3 <- 0
                               return(t3-tau3) }
  t3s <- seq(-1, 1, by=0.005)
  t4s <- NULL
  for(t3 in t3s) {
    rt  <- uniroot(fn, interval=c(-1,1), tau3=t3)
    t4  <- ((5*t3 / rt$root) - 1) / 4 # Hosking (2007)
    t4s <- c(t4s, t4)
  }
  t4s[is.nan(t4s)] <- 1/6 # by distribution properties

  plotlmrdia(lmrdia())
  points(t3s, t4s, pch=21, cex=0.5, bg=8, col="lightgreen")
  lines( t3s, t4s, col="darkgreen") # above GLO and see Hosking (2007, fig. 1)

  # eight powers as in Hosking and Wallis (1997) coefficient table for
  # many other distributions
  pdq3 <- stats::lm(t4s~I(t3s^1)+I(t3s^2)+I(t3s^3)+I(t3s^4)+
                        I(t3s^5)+I(t3s^6)+I(t3s^7)+I(t3s^8))
  lines(t3s, fitted.values(pdq3), lwd=2, col=grey(0.8))
  pdq3$coefficients # Ah, see the odd coefficients are near zero, so define
  # as such in a repeated linear model but with skips on the odd orders:
  pdq3 <- stats::lm(t4s~I(t3s^2)+I(t3s^4)+I(t3s^6)+I(t3s^8))
  lines(t3s, fitted.values(pdq3), lwd=1, col="red")
  max(abs(t4s - fitted.values(pdq3))) # show the max error in Rs resolution

  # we desire to compare "full resolution" to 5-digit truncation
  print(pdq3$coefficients,       16)  # in the 5 in the next line, c.2022,
  # we can  make new column in Asquith (2011, table 10.1) if ever needed for
  print(round(pdq3$coefficients,  5)) # a second edition

  t4pdq3 <- function(t3, use5digits=TRUE) { # helper to repeat within lmrdia()
    c05 <- c( 0.16688, 0, 0.98951, 0, -0.00526, 0, -0.24074, 0, 0.08906)
    c16 <- c( 0.166875136751297809, 0,  0.989506002306983601, 0,
             -0.005255434641059076, 0, -0.240744479052170501, 0,
              0.089060315246257210)
    ifelse(use5digits, myc <- c05, myc <- c16)
    t4 <- vector(mode="numeric", length(t3))
    for(i in 1:length(t3)) {
      t4[i] <- sum(sapply(2:length(myc), function(k) myc[k]*t3[i]^(k-1)))
    }
    return(t4 + myc[1]) # end with the intercept being added on
  }
  lines(t3s, t4pdq3(t3s), col="darkgreen", lty=2)
  summary(abs(t4s - t4pdq3(t3s, use5digits=TRUE )))
  summary(abs(t4s - t4pdq3(t3s, use5digits=FALSE)))
  max(    abs(t4s - t4pdq3(t3s, use5digits=TRUE )))
  max(    abs(t4s - t4pdq3(t3s, use5digits=FALSE)))
  # further comparisons as needed to understand the aforementioned operations
  plot(  t4s, t4s - t4pdq3(t3s, use5digits=TRUE ), col="red", type="l")
  lines( t4s, t4s - t4pdq3(t3s, use5digits=FALSE), col="blue")
  abline(h=0)
}
}
\references{
Asquith, W.H., 2011, Distributional analysis with L-moment statistics using the R environment for statistical computing: Createspace Independent Publishing Platform, ISBN 978--146350841--8.

Hosking, J.R.M., 2007, Distributions with maximum entropy subject to constraints on their L-moments or expected order statistics: Journal of Statistical Planning and Inference, v. 137, no. 9, pp. 2870--2891, \doi{10.1016/j.jspi.2006.10.010}.

Hosking, J.R.M., and Wallis, J.R., 1997, Regional frequency analysis---An approach based on L-moments: Cambridge University Press.
}
\author{W.H. Asquith}
\seealso{\code{\link{parpdq3}}, \code{\link{cdfpdq3}}, \code{\link{pdfpdq3}}, \code{\link{quapdq3}} }
\examples{
\dontrun{
  para <- list(para=c(20, 1, -0.5), type="pdq3")
  lmoms(quapdq3(runif(100000), para))$lambdas
  lmompdq3(para)$lambdas #}

\dontrun{
  para <- list(para=c(20, 1, +0.5), type="pdq3")
  lmoms(quapdq3(runif(100000), para))$lambdas
  lmompdq3(para)$lambdas #}
}
\keyword{L-moment (distribution)}
\keyword{Distribution: Polynomial Density-Quantile3}
