\name{quagld}
\alias{quagld}
\title{Quantile Function of the Generalized Lambda Distribution}
\description{
This function computes the quantiles of the Generalized Lambda distribution given parameters (\eqn{\xi}, \eqn{\alpha}, \eqn{\kappa}, and \eqn{h}) of the distribution computed by \code{\link{pargld}}. The quantile function is
\deqn{x(F) = \xi + \alpha(F^{\kappa} - (1-F)^{h}) \mbox{,}}{
      x(F) = x + a*(F^k-(1-F)^h)}
where \eqn{x(F)} is the quantile for nonexceedance probability \eqn{F},
\eqn{\xi} is a location parameter, \eqn{\alpha} is a scale
parameter, and \eqn{\kappa}, and \eqn{h} are shape parameters. Note that in this parameterization, the scale term is shown in the numerator and not the denominator. This is done for \pkg{lmomco} as part of the parallel nature between distributions whose various scale parameters are shown having the same units as the location parameter.
}
\usage{
quagld(f, para, paracheck=TRUE)
}
\arguments{
  \item{f}{Nonexceedance probability (\eqn{0 \le F \le 1}).}
  \item{para}{The parameters from \code{\link{pargld}} or \code{\link{vec2par}}.}
  \item{paracheck}{A logical controlling whether the parameters are checked for validity. Overriding of this check might be extremely important and needed for use quantile function in the context of TL-moments with nonzero trimming.}
}
\value{
  Quantile value for nonexceedance probability \eqn{F}.
}
\references{
Asquith, W.H., 2007, L-moments and TL-moments of the generalized lambda distribution: Computational Statistics and Data Analysis, v. 51, no. 9, pp. 4484--4496.

Karian, Z.A., and Dudewicz, E.J., 2000, Fitting statistical distributions---The generalized lambda distribution and generalized bootstrap methods:
CRC Press, Boca Raton, FL, 438 p.
}
\author{ W.H. Asquith}
\seealso{\code{\link{cdfgld}}, \code{\link{pargld}}, \code{\link{lmomgld}}, \code{\link{lmomTLgld}}, \code{\link{pargld}}, \code{\link{parTLgld}}
}
\examples{
\dontrun{
  para <- vec2par(c(123,34,4,3),type="gld")
  quagld(0.5,para, paracheck=FALSE)
}

}
\keyword{distribution}
\keyword{quantile function}
\keyword{Distribution: Generalized Lambda}
