\name{pdfgpa}
\alias{pdfgpa}
\title{Probability Density Function of the Generalized Pareto Distribution}
\description{
This function computes the probability density of the Generalized Pareto distribution given parameters (\eqn{\xi}, \eqn{\alpha}, and \eqn{\kappa}) computed by \code{\link{pargpa}}. The probability density function is
\deqn{f(x) = \alpha^{-1} \exp(-(1-\kappa)Y) \mbox{,}}
where \eqn{Y} is
\deqn{Y = -\kappa^{-1} \log\left(1 - \frac{\kappa(x - \xi)}{\alpha}\right)\mbox{,}}
for \eqn{\kappa \ne 0}, and
\deqn{Y = (x - \xi)/\alpha\mbox{,}}
for \eqn{\kappa = 0}, where \eqn{f(x)} is the probability density for quantile \eqn{x}, \eqn{\xi} is a location parameter, \eqn{\alpha} is a scale parameter, and \eqn{\kappa} is a shape parameter. The range of \eqn{x} is \eqn{\xi \le x \le \xi + \alpha/\kappa} if \eqn{k > 0}; \eqn{\xi \le x < \infty} if \eqn{\kappa \le 0}. Note that the shape parameter \eqn{\kappa} parameterization of the distribution herein follows that in tradition by the greater L-moment community and others use a sign reversal on \eqn{\kappa}. (The \pkg{evd} package is one example.)
}
\usage{
pdfgpa(x, para, paracheck=TRUE)
}
\arguments{
  \item{x}{A real value vector.}
  \item{para}{The parameters from \code{\link{pargpa}} or \code{\link{vec2par}}.}
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
\seealso{\code{\link{cdfgpa}}, \code{\link{quagpa}}, \code{\link{lmomgpa}}, \code{\link{pargpa}}}
\examples{
  lmr <- lmoms(c(123, 34, 4, 654, 37, 78))
  gpa <- pargpa(lmr)
  x   <- quagpa(0.5, gpa)
         pdfgpa(  x, gpa)

\dontrun{
  # We explore using maximum likelihood for GPA estimation on its density function
  # with stress testing near the K > -1 lower limit, K near zero, and then large K
  # producing extreme densities. We check the convergence and check on parameters
  # back estimating the mean. The experiment is designed that with repeated
  # operations that convergence "failures" in stats::optim()
  #   1  'indicates that the iteration limit maxit had been reached'
  #   10 'indicates degeneracy of the Nelder-Mead simplex.'
  # With the 10 being a bit more common and 1 but still for many runs convergence
  # at K = 8 is still attainable. Also, note the care in the construction of the
  # ptransf and pretransf for the honoring the GPA parameter space.
  small <- .Machine$double.eps; n <- 1000 # samples
  for(k in c(-1+small, -0.99, -1/2, -small, 0, 1/2, 8)) {
    names(k) <- "myKappa"
    gpa <- vec2par(c(2, 2, k), type="gpa")
    x <- rlmomco(n, gpa)
    mu1 <- mean(x); names(mu1) <- "mean"
    cv1 <-      NA; names(cv1) <- "converge"
    mle <- mle2par(x, type="gpa", para.int=pargpa(lmoms(x)),
             ptransf=function(t) { c(t[1], log(t[2]), log(t[3] +1)) },
           pretransf=function(t) { c(t[1], exp(t[2]), exp(t[3])-1)  },
                      null.on.not.converge=FALSE)
    mu2 <- lmomgpa(mle)$lambdas[1]; names(mu2) <- "backMean"
    cv2 <- mle$optim$convergence;   names(cv2) <- "converge"
    print(round(c(k, cv1, mu1, gpa$para), digits=5))
    print(round(c(k, cv2, mu2, mle$para), digits=5))
  } # }
}
\keyword{distribution}
\keyword{probability density function}
\keyword{Distribution: Generalized Pareto}
