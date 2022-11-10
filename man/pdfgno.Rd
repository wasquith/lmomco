\name{pdfgno}
\alias{pdfgno}
\title{Probability Density Function of the Generalized Normal Distribution}
\description{
This function computes the probability density of the Generalized Normal distribution given parameters (\eqn{\xi}, \eqn{\alpha}, and \eqn{\kappa})  computed by \code{\link{pargno}}. The probability density function function is
\deqn{f(x) = \frac{\exp(\kappa Y - Y^2/2)}{\alpha \sqrt{2\pi}} \mbox{,} }
where \eqn{Y} is
\deqn{Y = -\kappa^{-1} \log\left(1 - \frac{\kappa(x-\xi)}{\alpha}\right)\mbox{,}}
for \eqn{\kappa \ne 0}, and
\deqn{Y = (x-\xi)/\alpha\mbox{,}}
for \eqn{\kappa = 0}, where \eqn{f(x)} is the probability density for quantile \eqn{x}, \eqn{\xi} is a location parameter, \eqn{\alpha} is a scale parameter, and \eqn{\kappa} is a shape parameter.
}
\usage{
pdfgno(x, para)
}
\arguments{
  \item{x}{A real value vector.}
  \item{para}{The parameters from \code{\link{pargno}} or \code{\link{vec2par}}.}
}
\value{
  Probability density (\eqn{f}) for \eqn{x}.
}
\references{
Hosking, J.R.M., 1990, L-moments---Analysis and estimation of
distributions using linear combinations of order statistics: Journal
of the Royal Statistical Society, Series B, v. 52, pp. 105--124.

Hosking, J.R.M., 1996, FORTRAN routines for use with the method of L-moments:
Version 3, IBM Research Report RC20525, T.J. Watson Research Center,
Yorktown Heights, New York.

Hosking, J.R.M. and Wallis, J.R., 1997, Regional frequency analysis---An
approach based on L-moments: Cambridge University Press.
}
\author{W.H. Asquith}
\seealso{\code{\link{cdfgno}}, \code{\link{quagno}}, \code{\link{lmomgno}}, \code{\link{pargno}}, \code{\link{pdfln3}}
}
\examples{
  lmr <- lmoms(c(123,34,4,654,37,78))
  gno <- pargno(lmr)
  x <- quagno(0.5,gno)
  pdfgno(x,gno)
}
\keyword{distribution}
\keyword{probability density function}
\keyword{Distribution: Generalized Normal}
