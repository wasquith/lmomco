\name{cdfexp}
\alias{cdfexp}
\title{Cumulative Distribution Function of the Exponential Distribution}
\description{
This function computes the cumulative probability or nonexceedance probability of the Exponential distribution given parameters (\eqn{\xi} and \eqn{\alpha}  computed by \code{\link{parexp}}. The cumulative distribution function is
\deqn{F(x) = 1 - \exp(Y)\mbox{,}}
where \eqn{Y} is
\deqn{\frac{-(x - \xi)}{\alpha}\mbox{,}}
where \eqn{F(x)} is the nonexceedance probability for the quantile \eqn{x}, 
\eqn{\xi} is a location parameter, and \eqn{\alpha} is a scale parameter.
}
\usage{
cdfexp(x, para)
}
\arguments{
  \item{x}{A real value vector.}
  \item{para}{The parameters from \code{\link{parexp}} or  \code{\link{vec2par}}.}
}

\value{
  Nonexceedance probability (\eqn{F}) for \eqn{x}.
}
\references{
Hosking, J.R.M., 1990, L-moments---Analysis and estimation of
distributions using linear combinations of order statistics: Journal
of the Royal Statistical Society, Series B, v. 52, p. 105--124.

Hosking, J.R.M., 1996, FORTRAN routines for use with the method of L-moments:
Version 3, IBM Research Report RC20525, T.J. Watson Research Center,
Yorktown Heights, New York.

Hosking, J.R.M., and Wallis, J.R., 1997, Regional frequency analysis---An
approach based on L-moments: Cambridge University Press.
}
\author{ W.H. Asquith}
\seealso{\code{\link{pdfexp}}, \code{\link{quaexp}}, \code{\link{lmomexp}}, \code{\link{parexp}}}
\examples{
  lmr <- lmoms(c(123,34,4,654,37,78))
  cdfexp(50,parexp(lmr))
}
\keyword{distribution}
\keyword{cumulative distribution function}
\keyword{Distribution: Exponential}
