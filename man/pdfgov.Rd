\name{pdfgov}
\alias{pdfgov}
\title{Probability Density Function of the Govindarajulu Distribution}
\description{
This function computes the probability density of the Govindarajulu distribution given parameters (\eqn{\xi}, \eqn{\alpha}, and \eqn{\beta})  computed by \code{\link{pargov}}. The probability density function is
\deqn{f(x) = [\alpha\beta(\beta+1)]^{-1} [F(x)]^{1-\beta} [1 - F(x)]^{-1} \mbox{,}}
where \eqn{f(x)} is the probability density for quantile \eqn{x}, \eqn{F(x)} the cumulative distribution function  or nonexceedance probability at \eqn{x}, \eqn{\xi} is a location parameter, \eqn{\alpha} is a scale parameter, and \eqn{\beta} is a shape parameter.
}
\usage{
pdfgov(x, para)
}
\arguments{
  \item{x}{A real value vector.}
  \item{para}{The parameters from \code{\link{pargov}} or \code{\link{vec2par}}.}
}

\value{
  Probability density (\eqn{f}) for \eqn{x}.
}
\references{
Gilchrist, W.G., 2000, Statistical modelling with quantile functions: Chapman and Hall/CRC, Boca Raton.

Nair, N.U., Sankaran, P.G., Balakrishnan, N., 2013, Quantile-based reliability analysis: Springer, New York.

Nair, N.U., Sankaran, P.G., and Vineshkumar, B., 2012, The Govindarajulu distribution---Some Properties and applications: Communications in Statistics, Theory and Methods, v. 41, no. 24, pp. 4391--4406.
}
\author{ W.H. Asquith}
\seealso{\code{\link{cdfgov}}, \code{\link{quagov}}, \code{\link{lmomgov}}, \code{\link{pargov}}}
\examples{
  lmr <- lmoms(c(123,34,4,654,37,78))
  gov <- pargov(lmr)
  x <- quagov(0.5,gov)
  pdfgov(x,gov)
}
\keyword{distribution}
\keyword{probability density function}
\keyword{Distribution: Govindarajulu}
