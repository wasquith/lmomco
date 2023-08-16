\name{pdfsmd}
\alias{pdfsmd}
\title{Probability Density Function of the Singh--Maddala Distribution}
\description{
This function computes the probability density of the Singh--Maddala (Burr Type XII) distribution given parameters (\eqn{a}, \eqn{b}, and \eqn{q})  computed by \code{\link{parsmd}}. The probability density function is

\deqn{f(x) = \frac{bqx^{b-1}}{a^b \biggl(1 + \bigl((x-\xi)/a\bigr)^b \biggr)^{q+1}}\mbox{,}}

where \eqn{f(x)} is the probability density for quantile \eqn{x} with \eqn{0 \le x \le \infty}, \eqn{\xi} is a location parameter, \eqn{a} is a scale parameter (\eqn{a > 0}), \eqn{b} is a shape parameter (\eqn{b > 0}), and \eqn{q} is another shape parameter (\eqn{q > 0}).
}
\usage{
pdfsmd(x, para)
}
\arguments{
  \item{x}{A real value vector.}
  \item{para}{The parameters from \code{\link{parsmd}} or \code{\link{vec2par}}.}
}
\value{
  Probability density (\eqn{f}) for \eqn{x}.
}
\references{
Shahzad, M.N., and Zahid, A., 2013, Parameter estimation of Singh Maddala distribution by moments: International Journal of Advanced Statistics and Probability, v. 1, no. 3, pp. 121--131, \doi{10.14419/ijasp.v1i3.1206}.
}
\author{W.H. Asquith}
\seealso{\code{\link{cdfsmd}}, \code{\link{quasmd}}, \code{\link{lmomsmd}}, \code{\link{parsmd}}}
\examples{

}
\keyword{distribution}
\keyword{probability density function}
\keyword{Distribution: Singh-Maddala}
\keyword{Distribution: Burr Type XII}
