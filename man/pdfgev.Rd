\name{pdfgev}
\alias{pdfgev}
\title{Probability Density Function of the Generalized Extreme Value Distribution}
\description{
This function computes the probability density
of the Generalized Extreme Value distribution given parameters (\eqn{\xi}, \eqn{\alpha}, and \eqn{\kappa})  computed by \code{\link{pargev}}. The probability density function is
\deqn{f(x) = \alpha^{-1} \exp[-(1-\kappa)Y - \exp(-Y)] \mbox{,}}
where \eqn{Y} is
\deqn{Y = -\kappa^{-1} \log\!\left(1 - \frac{\kappa(x-\xi)}{\alpha}\right)\mbox{,}}
for \eqn{\kappa \ne 0}, and
\deqn{Y = (x-\xi)/\alpha\mbox{,}}
for \eqn{\kappa = 0}, where \eqn{f(x)} is the probability density for quantile \eqn{x}, \eqn{\xi} is a location parameter, \eqn{\alpha} is a scale parameter, and \eqn{\kappa} is a shape parameter. The range of \eqn{x} is \eqn{-\infty < x \le \xi + \alpha/\kappa} if \eqn{k > 0}; \eqn{\xi + \alpha/\kappa \le x < \infty} if \eqn{\kappa \le 0}. Note that the shape parameter \eqn{\kappa} parameterization of the distribution herein follows that in tradition by the greater L-moment community and others use a sign reversal on \eqn{\kappa}. (The \pkg{evd} package is one example.)
}
\usage{
pdfgev(x, para, paracheck=TRUE)
}
\arguments{
  \item{x}{A real value vector.}
  \item{para}{The parameters from \code{\link{pargev}} or \code{\link{vec2par}}.}
  \item{paracheck}{A logical switch as to whether the validity of the parameters should be checked.}
}

\value{
  Probability density (\eqn{f}) for \eqn{x}.
}
\references{
Hosking, J.R.M., 1990, L-moments---Analysis and estimation of distributions using linear combinations of order statistics: Journal of the Royal Statistical Society, Series B, v. 52, pp. 105--124, \doi{10.1111/j.2517-6161.1990.tb01775.x}.

Hosking, J.R.M., 1996, FORTRAN routines for use with the method of L-moments: Version 3, IBM Research Report RC20525, T.J. Watson Research Center, Yorktown Heights, New York.

Hosking, J.R.M. and Wallis, J.R., 1997, Regional frequency analysis---An approach based on L-moments: Cambridge University Press.
}
\author{W.H. Asquith}
\seealso{\code{\link{cdfgev}}, \code{\link{quagev}}, \code{\link{lmomgev}}, \code{\link{pargev}}}
\examples{
  lmr <- lmoms(c(123, 34, 4, 654, 37, 78))
  gev <- pargev(lmr)
  x <- quagev(0.5, gev)
       pdfgev(  x, gev)

\dontrun{
  # We explore using maximum likelihood for GEV estimation on its density function.
  # We check the convergence and check on parameters back estimating the mean.
  small <- .Machine$double.eps
  for(k in c(-2, -1/2, -small, 0, +small, 1/2, 2)) {
    names(k) <- "myKappa"
    gev <- vec2par(c(2, 2, k), type="gev")
    x <- rlmomco(1000, gev)
    mu1 <- mean(x); names(mu1) <- "mean"
    cv1 <-      NA; names(cv1) <- "converge"
    mle <- mle2par(x, type="gev", para.int=pargev(lmoms(x)),
             ptransf=function(t) { c(t[1], log(t[2]), t[3]) },
           pretransf=function(t) { c(t[1], exp(t[2]), t[3]) },
                      null.on.not.converge=FALSE)
    mu2 <- lmomgev(mle)$lambdas[1]; names(mu2) <- "backMean"
    cv2 <- mle$optim$convergence;   names(cv2) <- "converge"
    print(round(c(k, cv1, mu1, gev$para), digits=5))
    print(round(c(k, cv2, mu2, mle$para), digits=5))
  } # }
}
\keyword{distribution}
\keyword{probability density function}
\keyword{Distribution: Generalized Extreme Value}
