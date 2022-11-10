\name{quakap}
\alias{quakap}
\title{Quantile Function of the Kappa Distribution}
\description{
This function computes the quantiles of the Kappa distribution given parameters (\eqn{\xi}, \eqn{\alpha}, \eqn{\kappa}, and \eqn{h})  computed by \code{\link{parkap}}. The quantile function  is
\deqn{x(F) = \xi + \frac{\alpha}{\kappa}\left(1-{\left(\frac{1-F^h}{h}\right)}^\kappa\right) \mbox{,}}
where \eqn{x(F)} is the quantile for nonexceedance probability \eqn{F}, \eqn{\xi} is a location parameter, \eqn{\alpha} is a scale parameter, \eqn{\kappa} is a shape parameter, and \eqn{h} is another shape parameter.
}
\usage{
quakap(f, para, paracheck=TRUE)
}
\arguments{
  \item{f}{Nonexceedance probability (\eqn{0 \le F \le 1}).}
  \item{para}{The parameters from \code{\link{parkap}} or \code{\link{vec2par}}.}
  \item{paracheck}{A logical controlling whether the parameters are checked for validity. Overriding of this check might be extremely important and needed for use of the  quantile function in the context of TL-moments with nonzero trimming.}
}
\value{
  Quantile value for nonexceedance probability \eqn{F}.
}
\references{
Hosking, J.R.M., 1994, The four-parameter kappa distribution: IBM Journal of Reserach and Development, v. 38, no. 3, pp. 251--258.

Hosking, J.R.M., and Wallis, J.R., 1997, Regional frequency analysis---An approach based on L-moments: Cambridge University Press.
}
\author{W.H. Asquith}
\seealso{\code{\link{cdfkap}}, \code{\link{pdfkap}}, \code{\link{lmomkap}}, \code{\link{parkap}}}
\examples{
  lmr <- lmoms(c(123,34,4,654,37,78,21,32,231,23))
  quakap(0.5,parkap(lmr))
}
\keyword{distribution}
\keyword{quantile function}
\keyword{Distribution: Kappa}
