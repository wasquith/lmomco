\name{cdfgev}
\alias{cdfgev}
\title{Cumulative Distribution Function of the Generalized Extreme Value Distribution}
\description{
This function computes the cumulative probability or nonexceedance probability of the Generalized Extreme Value distribution given parameters (\eqn{\xi}, \eqn{\alpha}, and \eqn{\kappa}) computed by \code{\link{pargev}}. The cumulative distribution function is
\deqn{F(x) = \mathrm{exp}(-\mathrm{exp}(-Y)) \mbox{,}}
where \eqn{Y} is
\deqn{Y = -\kappa^{-1} \log\left(1 - \frac{\kappa(x-\xi)}{\alpha}\right)\mbox{,}}
for \eqn{\kappa \ne 0} and
\deqn{Y = (x-\xi)/\alpha\mbox{,}}
for \eqn{\kappa = 0}, where \eqn{F(x)} is the nonexceedance probability for quantile \eqn{x}, \eqn{\xi} is a location parameter, \eqn{\alpha} is a scale parameter, and \eqn{\kappa} is a shape parameter. The range of \eqn{x} is \eqn{-\infty < x \le \xi + \alpha/\kappa} if \eqn{k > 0}; \eqn{\xi + \alpha/\kappa \le x < \infty} if \eqn{\kappa \le 0}. Note that the shape parameter \eqn{\kappa} parameterization of the distribution herein follows that in tradition by the greater L-moment community and others use a sign reversal on \eqn{\kappa}. (The \pkg{evd} package is one example.)
}
\usage{
cdfgev(x, para, paracheck=TRUE)
}
\arguments{
  \item{x}{A real value vector.}
  \item{para}{The parameters from \code{\link{pargev}} or  \code{\link{vec2par}}.}
  \item{paracheck}{A logical switch as to whether the validity of the parameters should be checked.}
}
\value{
  Nonexceedance probability (\eqn{F}) for \eqn{x}.
}
\references{
Hosking, J.R.M., 1990, L-moments---Analysis and estimation of distributions using linear combinations of order statistics: Journal of the Royal Statistical Society, Series B, v. 52, pp. 105--124, \doi{10.1111/j.2517-6161.1990.tb01775.x}.

Hosking, J.R.M., 1996, FORTRAN routines for use with the method of L-moments: Version 3, IBM Research Report RC20525, T.J. Watson Research Center, Yorktown Heights, New York.

Hosking, J.R.M., and Wallis, J.R., 1997, Regional frequency analysis---An approach based on L-moments: Cambridge University Press.
}
\author{W.H. Asquith}
\seealso{\code{\link{pdfgev}}, \code{\link{quagev}}, \code{\link{lmomgev}}, \code{\link{pargev}}}
\examples{
  lmr <- lmoms(c(123, 34, 4, 654, 37, 78))
  cdfgev(50, pargev(lmr))
}
\keyword{distribution}
\keyword{cumulative distribution function}
\keyword{Distribution: Generalized Extreme Value}
