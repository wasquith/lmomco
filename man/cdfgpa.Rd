\name{cdfgpa}
\alias{cdfgpa}
\title{Cumulative Distribution Function of the Generalized Pareto Distribution}
\description{
This function computes the cumulative probability or nonexceedance probability of the Generalized Pareto distribution given parameters (\eqn{\xi}, \eqn{\alpha}, and \eqn{\kappa})  computed by \code{\link{pargpa}}. The cumulative distribution function is
\deqn{F(x) = 1 - \mathrm{exp}(-Y) \mbox{,}}
where \eqn{Y} is
\deqn{Y = -\kappa^{-1} \log\left(1 - \frac{\kappa(x-\xi)}{\alpha}\right)\mbox{,}}
for \eqn{\kappa \ne 0} and
\deqn{Y = (x-\xi)/\alpha\mbox{,}}
for \eqn{\kappa = 0}, where \eqn{F(x)} is the nonexceedance probability for quantile \eqn{x}, \eqn{\xi} is a location parameter, \eqn{\alpha} is a scale parameter, and \eqn{\kappa} is a shape parameter.
}
\usage{
cdfgpa(x, para)
}
\arguments{
  \item{x}{A real value vector.}
  \item{para}{The parameters from \code{\link{pargpa}} or  \code{\link{vec2par}}.}
}

\value{
  Nonexceedance probability (\eqn{F}) for \eqn{x}.
}
\references{
Hosking, J.R.M., 1990, L-moments---Analysis and estimation of
distributions using linear combinations of order statistics: Journal
of the Royal Statistical Society, Series B, v. 52, pp. 105--124.

Hosking, J.R.M., 1996, FORTRAN routines for use with the method of L-moments:
Version 3, IBM Research Report RC20525, T.J. Watson Research Center,
Yorktown Heights, New York.

Hosking, J.R.M., and Wallis, J.R., 1997, Regional frequency analysis---An
approach based on L-moments: Cambridge University Press.
}
\author{ W.H. Asquith}
\seealso{\code{\link{pdfgpa}}, \code{\link{quagpa}}, \code{\link{lmomgpa}}, \code{\link{pargpa}}}
\examples{
  lmr <- lmoms(c(123,34,4,654,37,78))
  cdfgpa(50,pargpa(lmr))
}
\keyword{distribution}
\keyword{cumulative distribution function}
\keyword{Distribution: Generalized Pareto}
