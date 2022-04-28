\name{quagev}
\alias{quagev}
\title{Quantile Function of the Generalized Extreme Value Distribution}
\description{
This function computes the quantiles of the Generalized Extreme Value distribution given parameters (\eqn{\xi}, \eqn{\alpha}, and \eqn{\kappa}) of the distribution computed by \code{\link{pargev}}. The quantile function of the distribution is
\deqn{x(F) = \xi + \frac{\alpha}{\kappa} \left( 1-(-\log(F))^\kappa \right)\mbox{,}}
for \eqn{\kappa \ne 0}, and
\deqn{x(F) = \xi - \alpha \log(-\log(F))\mbox{,}}
for \eqn{\kappa = 0},  where \eqn{x(F)} is the quantile for nonexceedance probability \eqn{F}, \eqn{\xi} is a location parameter, \eqn{\alpha} is a scale parameter, and \eqn{\kappa} is a shape parameter. The range of \eqn{x} is \eqn{-\infty < x \le \xi + \alpha/\kappa} if \eqn{k > 0}; \eqn{\xi + \alpha/\kappa \le x < \infty} if \eqn{\kappa \le 0}. Note that the shape parameter \eqn{\kappa} parameterization of the distribution herein follows that in tradition by the greater L-moment community and others use a sign reversal on \eqn{\kappa}. (The \pkg{evd} package is one example.)
}
\usage{
quagev(f, para, paracheck=TRUE)
}
\arguments{
  \item{f}{Nonexceedance probability (\eqn{0 \le F \le 1}).}
  \item{para}{The parameters from \code{\link{pargev}} or \code{\link{vec2par}}.}
  \item{paracheck}{A logical controlling whether the parameters are checked for validity. Overriding of this check might be extremely important and needed for use of the quantile function in the context of TL-moments with nonzero trimming.}
}
\value{
  Quantile value for nonexceedance probability \eqn{F}.
}
\references{
Hosking, J.R.M., 1990, L-moments---Analysis and estimation of distributions using linear combinations of order statistics: Journal of the Royal Statistical Society, Series B, v. 52, pp. 105--124, \doi{10.1111/j.2517-6161.1990.tb01775.x}.

Hosking, J.R.M., 1996, FORTRAN routines for use with the method of L-moments: Version 3, IBM Research Report RC20525, T.J. Watson Research Center, Yorktown Heights, New York.

Hosking, J.R.M., and Wallis, J.R., 1997, Regional frequency analysis---An approach based on L-moments: Cambridge University Press.
}
\author{W.H. Asquith}
\seealso{\code{\link{cdfgev}}, \code{\link{pdfgev}}, \code{\link{lmomgev}}, \code{\link{pargev}}}
\examples{
  lmr <- lmoms(c(123, 34, 4, 654, 37, 78))
  quagev(0.5, pargev(lmr))
}
\keyword{distribution}
\keyword{quantile function}
\keyword{Distribution: Generalized Extreme Value}
