\name{cdfray}
\alias{cdfray}
\title{Cumulative Distribution Function of the Rayleigh Distribution}
\description{
This function computes the cumulative probability or nonexceedance probability of the Rayleigh distribution given parameters (\eqn{\xi} and \eqn{\alpha}) computed by \code{\link{parray}}. The cumulative distribution function is
\deqn{F(x) =  1 - \mathrm{exp}[-(x - \xi)^2/(2\alpha^2)]\mbox{,}} where \eqn{F(x)} is the nonexceedance probability for quantile \eqn{x}, \eqn{\xi} is a location parameter, and \eqn{\alpha} is a scale parameter.
}
\usage{
cdfray(x, para)
}
\arguments{
  \item{x}{A real value vector.}
  \item{para}{The parameters from \code{\link{parray}} or  \code{\link{vec2par}}.}
}

\value{
  Nonexceedance probability (\eqn{F}) for \eqn{x}.
}
\references{
Hosking, J.R.M., 1986, The theory of probability weighted moments:
Research Report RC12210, IBM Research Division, Yorkton Heights, N.Y.
}
\author{ W.H. Asquith}
\seealso{\code{\link{pdfray}}, \code{\link{quaray}}, \code{\link{lmomray}}, \code{\link{parray}}}
\examples{
  lmr <- lmoms(c(123,34,4,654,37,78))
  cdfray(50,parray(lmr))
}
\keyword{distribution}
\keyword{cumulative distribution function}
\keyword{Distribution: Rayleigh}
