\name{quaray}
\alias{quaray}
\title{Quantile Function of the Rayleigh Distribution}
\description{
This function computes the quantiles of the Rayleigh distribution given parameters (\eqn{\xi} and \eqn{\alpha})  computed by \code{\link{parray}}. The quantile function  is
\deqn{x(F) = \xi + \sqrt{-2\alpha^2\log(1-F)} \mbox{,}}
where \eqn{x(F)} is the quantile for nonexceedance probability \eqn{F}, \eqn{\xi} is a location parameter, and \eqn{\alpha} is a scale parameter.
}
\usage{
quaray(f, para, paracheck=TRUE)
}
\arguments{
  \item{f}{Nonexceedance probability (\eqn{0 \le F \le 1}).}
  \item{para}{The parameters from \code{\link{parray}} or \code{\link{vec2par}}.}
  \item{paracheck}{A logical controlling whether the parameters are checked for validity. Overriding of this check might be extremely important and needed for use of the quantile function in the context of TL-moments with nonzero trimming.}
}
\value{
  Quantile value for nonexceedance probability \eqn{F}.
}
\references{
Hosking, J.R.M., 1986, The theory of probability weighted moments:
Research Report RC12210, IBM Research Division, Yorkton Heights, N.Y.
}
\author{ W.H. Asquith}
\seealso{\code{\link{cdfray}}, \code{\link{pdfray}}, \code{\link{lmomray}}, \code{\link{parray}}}
\examples{
  lmr <- lmoms(c(123,34,4,654,37,78))
  quaray(0.5,parray(lmr))
}
\keyword{distribution}
\keyword{quantile function}
\keyword{Distribution: Rayleigh}
