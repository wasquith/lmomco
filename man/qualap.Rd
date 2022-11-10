\name{qualap}
\alias{qualap}
\title{Quantile Function of the Laplace Distribution}
\description{
This function computes the quantiles of the Laplace distribution given parameters (\eqn{\xi} and \eqn{\alpha})  computed by \code{\link{parlap}}. The quantile function  is
\deqn{x(F) = \xi + \alpha\times\log(2F)\mbox{,}}
for \eqn{F \le 0.5}, and
\deqn{x(F) = \xi - \alpha\times\log(2(1-F))\mbox{,}}
for \eqn{F > 0.5}, where \eqn{x(F)} is the quantile for nonexceedance probability \eqn{F}, \eqn{\xi} is a location parameter, and \eqn{\alpha} is a scale parameter.
}
\usage{
qualap(f, para, paracheck=TRUE)
}
\arguments{
  \item{f}{Nonexceedance probability (\eqn{0 \le F \le 1}).}
  \item{para}{The parameters from \code{\link{parlap}} or \code{\link{vec2par}}.}
  \item{paracheck}{A logical controlling whether the parameters are checked for validity. Overriding of this check might be extremely important and needed for use of the quantile function in the context of TL-moments with nonzero trimming.}
}
\value{
  Quantile value for for nonexceedance probability \eqn{F}.
}
\references{
Hosking, J.R.M., 1986, The theory of probability weighted moments: IBM Research Report RC12210, T.J. Watson Research Center, Yorktown Heights, New York.
}
\author{W.H. Asquith}
\seealso{\code{\link{cdflap}}, \code{\link{pdflap}}, \code{\link{lmomlap}}, \code{\link{parlap}}}
\examples{
  lmr <- lmoms(c(123,34,4,654,37,78))
  qualap(0.5,parlap(lmr))
}
\keyword{distribution}
\keyword{quantile function}
\keyword{Distribution: Laplace}
