\name{pdfsla}
\alias{pdfsla}
\title{Probability Density Function of the Slash Distribution}
\description{
This function computes the probability density of the Slash distribution given parameters (\eqn{\xi} and \eqn{\alpha})  provided by \code{\link{parsla}}.  The probability density function is
\deqn{f(x) = \frac{\phi(0) - \phi(y)}{y^2} \mbox{,}}
where \eqn{f(x)} is the probability density for quantile \eqn{x}, \eqn{y = (x - \xi)/\alpha}, \eqn{\xi} is a location parameter, and \eqn{\alpha} is a scale parameter. The function \eqn{\phi(y)} is the probability density function of the Standard Normal distribution.
}
\usage{
pdfsla(x, para)
}
\arguments{
  \item{x}{A real value vector.}
  \item{para}{The parameters from \code{\link{parsla}} or \code{\link{vec2par}}.}
}

\value{
  Probability density (\eqn{f}) for \eqn{x}.
}
\references{
Rogers, W.H., and Tukey, J.W., 1972, Understanding some long-tailed symmetrical distributions: Statistica Neerlandica, v. 26, no. 3, pp. 211--226.
}
\author{W.H. Asquith}
\seealso{\code{\link{cdfsla}}, \code{\link{quasla}}, \code{\link{lmomsla}}, \code{\link{parsla}}}
\examples{
  sla <- vec2par(c(12,1.2),type='sla')
  x <- quasla(0.5,sla)
  pdfsla(x,sla)
}
\keyword{distribution}
\keyword{probability density function}
\keyword{Distribution: Slash}
