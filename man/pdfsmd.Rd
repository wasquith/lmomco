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
# The SMD approximating the normal and use x=0
tau4_of_normal <- 30 * pi^-1 * atan(sqrt(2)) - 9 # from theory
pdfsmd(0, parsmd( vec2lmom( c( -pi, pi, 0,tau4_of_normal ) ) ) ) # 0.061953
dnorm( 0, mean=-pi, sd=pi*sqrt(pi))                              # 0.061103
}
\keyword{distribution}
\keyword{probability density function}
\keyword{Distribution: Singh-Maddala}
\keyword{Distribution: Burr Type XII}