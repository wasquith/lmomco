\name{cdfkap}
\alias{cdfkap}
\title{Cumulative Distribution Function of the Kappa Distribution}
\description{
This function computes the cumulative probability or nonexceedance probability
of the Kappa of the distribution computed by \code{\link{parkap}}. The cumulative distribution function is

\deqn{F(x) = \left(1-h\left(1-\frac{\kappa(x-\xi)}{\alpha}\right)^{1/\kappa}\right)^{1/h} \mbox{,}}

where \eqn{F(x)} is the nonexceedance probability for quantile \eqn{x},
\eqn{\xi} is a location parameter, \eqn{\alpha} is a scale parameter,
\eqn{\kappa} is a shape parameter, and \eqn{h} is another shape parameter.
}
\usage{
cdfkap(x, para)
}
\arguments{
  \item{x}{A real value vector.}
  \item{para}{The parameters from \code{\link{parkap}} or  \code{\link{vec2par}}.}
}

\value{
  Nonexceedance probability (\eqn{F}) for \eqn{x}.
}
\references{
Hosking, J.R.M., 1994, The four-parameter kappa distribution: IBM Journal of Reserach and Development, v. 38, no. 3, pp. 251--258.

Hosking, J.R.M., and Wallis, J.R., 1997, Regional frequency analysis---An approach based on L-moments: Cambridge University Press.
}
\author{W.H. Asquith}
\seealso{\code{\link{pdfkap}}, \code{\link{quakap}}, \code{\link{lmomkap}}, \code{\link{parkap}}}
\examples{
  lmr <- lmoms(c(123,34,4,654,37,78,21,32,231,23))
  cdfkap(50,parkap(lmr))
}
\keyword{distribution}
\keyword{cumulative distribution function}
\keyword{Distribution: Kappa}
