\name{pdflap}
\alias{pdflap}
\title{Probability Density Function of the Laplace Distribution}
\description{
This function computes the probability density of the Laplace distribution given parameters (\eqn{\xi} and \eqn{\alpha})  computed by \code{\link{parlap}}. The probability density function is
\deqn{f(x) = (2\alpha)^{-1} \exp(Y)\mbox{,}}
where \eqn{Y} is
\deqn{Y = \left(\frac{-|x - \xi|}{\alpha}\right)\mbox{.}}
}
\usage{
pdflap(x, para)
}
\arguments{
  \item{x}{A real value vector.}
  \item{para}{The parameters from \code{\link{parlap}} or \code{\link{vec2par}}.}
}

\value{
  Probability density (\eqn{f}) for \eqn{x}.
}
\references{
Hosking, J.R.M., 1986, The theory of probability weighted moments: IBM Research Report RC12210, T.J. Watson Research Center, Yorktown Heights, New York.
}
\author{W.H. Asquith}
\seealso{\code{\link{cdflap}}, \code{\link{qualap}}, \code{\link{lmomlap}}, \code{\link{parlap}}}
\examples{
  lmr <- lmoms(c(123,34,4,654,37,78))
  lap <- parlap(lmr)
  x <- qualap(0.5,lap)
  pdflap(x,lap)
}
\keyword{distribution}
\keyword{probability density function}
\keyword{Distribution: Laplace}
