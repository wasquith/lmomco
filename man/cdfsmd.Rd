\name{cdfsmd}
\alias{cdfsmd}
\title{Cumulative Distribution Function of the Singh--Maddala Distribution}
\description{
This function computes the cumulative probability or nonexceedance probability of the Singh--Maddala (Burr Type XII) distribution given parameters (\eqn{a}, \eqn{b}, and \eqn{q}) of the distribution computed by \code{\link{parsmd}}. The cumulative distribution function is
\deqn{F(x) = 1 - \biggl(1 + ((x-\xi)/a)^b \biggl)^{-q}\mbox{,}}
where \eqn{F(x)} is the nonexceedance probability for quantile \eqn{x} with \eqn{0 \le x \le \infty}, \eqn{\xi} is a location parameter, \eqn{a} is a scale parameter (\eqn{a > 0}), \eqn{b} is a shape parameter (\eqn{b > 0}), and \eqn{q} is another shape parameter (\eqn{q > 0}).
}
\usage{
cdfsmd(x, para)
}
\arguments{
  \item{x}{A real value vector.}
  \item{para}{The parameters from \code{\link{parsmd}} or  \code{\link{vec2par}}.}
}

\value{
  Nonexceedance probability (\eqn{F}) for \eqn{x}.
}
\references{
Shahzad, M.N., and Zahid, A., 2013, Parameter estimation of Singh Maddala distribution by moments: International Journal of Advanced Statistics and Probability, v. 1, no. 3, pp. 121--131, \doi{10.14419/ijasp.v1i3.1206}.
}
\author{W.H. Asquith}
\seealso{\code{\link{pdfsmd}}, \code{\link{quasmd}}, \code{\link{lmomsmd}}, \code{\link{parsmd}}}
\examples{

}
\keyword{distribution}
\keyword{cumulative distribution function}
\keyword{Distribution: Singh-Maddala}
\keyword{Distribution: Burr Type XII}
