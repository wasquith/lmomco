\name{cdfgum}
\alias{cdfgum}
\title{Cumulative Distribution Function of the Gumbel Distribution}
\description{
This function computes the cumulative probability or nonexceedance probability of the Gumbel distribution given parameters (\eqn{\xi} and \eqn{\alpha}) computed by \code{\link{pargum}}. The cumulative distribution function is
\deqn{F(x) = \mathrm{exp}(-\mathrm{exp}(Y)) \mbox{,}}
where
\deqn{Y = -\frac{x - \xi}{\alpha} \mbox{,}}
where \eqn{F(x)} is the nonexceedance probability for quantile \eqn{x}, \eqn{\xi} is a location parameter, and \eqn{\alpha} is a scale parameter.
}
\usage{
cdfgum(x, para)
}
\arguments{
  \item{x}{A real value vector.}
  \item{para}{The parameters from \code{\link{pargum}} or  \code{\link{vec2par}}.}
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
\seealso{\code{\link{pdfgum}}, \code{\link{quagum}}, \code{\link{lmomgum}}, \code{\link{pargum}}}
\examples{
  lmr <- lmoms(c(123,34,4,654,37,78))
  cdfgum(50,pargum(lmr))
}
\keyword{distribution}
\keyword{cumulative distribution function}
\keyword{Distribution: Gumbel}
