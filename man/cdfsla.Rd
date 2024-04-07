\name{cdfsla}
\alias{cdfsla}
\title{Cumulative Distribution Function of the Slash Distribution}
\description{
This function computes the cumulative probability or nonexceedance probability of the Slash distribution given parameters (\eqn{\xi} and \eqn{\alpha}) of the distribution provided by \code{\link{parsla}} or \code{\link{vec2par}}.  The cumulative distribution function is
\deqn{F(x) = \Phi(Y) - [\phi(0) - \phi(Y)]/Y \mbox{,}}
for \eqn{Y \ne 0} and
\deqn{F(x) = 1/2 \mbox{,}}
for \eqn{Y = 0}, where \eqn{F(x)} is the nonexceedance probability for quantile \eqn{x}, \eqn{Y = (x - \xi)/\alpha}, \eqn{\xi} is a location parameter, and \eqn{\alpha} is a scale parameter. The function \eqn{\Phi(Y)} is the cumulative distribution function of the Standard Normal distribution, and \eqn{\phi(Y)} is the probability density function of the Standard Normal distribution.
}
\usage{
cdfsla(x, para)
}
\arguments{
  \item{x}{A real value vector.}
  \item{para}{The parameters from \code{\link{parsla}} or \code{\link{vec2par}}.}
}

\value{
  Nonexceedance probability (\eqn{F}) for \eqn{x}.
}
\references{
Rogers, W.H., and Tukey, J.W., 1972, Understanding some long-tailed symmetrical distributions: Statistica Neerlandica, v. 26, no. 3, pp. 211--226.
}
\author{W.H. Asquith}
\seealso{\code{\link{pdfsla}}, \code{\link{quasla}}, \code{\link{lmomsla}}, \code{\link{parsla}}}
\examples{
  para <- c(12, 1.2)
  cdfsla(50, vec2par(para, type="sla"))
}
\keyword{distribution}
\keyword{cumulative distribution function}
\keyword{Distribution: Slash}
