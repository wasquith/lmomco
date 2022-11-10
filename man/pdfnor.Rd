\name{pdfnor}
\alias{pdfnor}
\title{Probability Density Function of the Normal Distribution}
\description{
This function computes the probability density function
of the Normal distribution given parameters  computed
by \code{\link{parnor}}.  The probability density function is
\deqn{f(x) = \frac{1}{\sigma \sqrt{2\pi}} \exp\!\left(\frac{-(x-\mu)^2}{2\sigma^2}\right) \mbox{,}}
where \eqn{f(x)} is the probability density for quantile \eqn{x}, \eqn{\mu} is the arithmetic mean, and \eqn{\sigma} is the standard deviation. The \R function \code{pnorm} is used.
}
\usage{
pdfnor(x, para)
}
\arguments{
  \item{x}{A real value.}
  \item{para}{The parameters from \code{\link{parnor}} or \code{\link{vec2par}}.}
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

Hosking, J.R.M., and Wallis, J.R., 1997, Regional frequency analysis---An
approach based on L-moments: Cambridge University Press.
}
\author{W.H. Asquith}
\seealso{\code{\link{cdfnor}}, \code{\link{quanor}}, \code{\link{lmomnor}}, \code{\link{parnor}}}
\examples{
  lmr <- lmoms(c(123,34,4,654,37,78))
  pdfnor(50,parnor(lmr))
}
\keyword{distribution}
\keyword{probability density function}
\keyword{Distribution: Normal}
