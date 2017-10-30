\name{cdftri}
\alias{cdftri}
\title{Cumulative Distribution Function of the Asymmetric Triangular Distribution}
\description{
This function computes the cumulative probability or nonexceedance probability of the Asymmetric Triangular distribution given parameters (\eqn{\nu}, \eqn{\omega}, and \eqn{\psi}) computed by \code{\link{partri}}. The cumulative distribution function is
\deqn{F(x) = \frac{(x - \nu)^2}{(\omega-\nu)(\psi-\nu)}\mbox{,}}
for \eqn{x < \omega},
\deqn{F(x) = 1 - \frac{(\psi - x)^2}{(\psi - \omega)(\psi - \nu)}\mbox{,}}
for \eqn{x >  \omega}, and
\deqn{F(x) = \frac{(\omega - \nu)}{(\psi - \nu)}\mbox{,}}
for \eqn{x = \omega}
where \eqn{x(F)} is the quantile for nonexceedance probability \eqn{F}, \eqn{\nu} is the minimum, \eqn{\psi} is the maximum, and \eqn{\omega} is the mode of the distribution.
}
\usage{
cdftri(x, para)
}
\arguments{
  \item{x}{A real value vector.}
  \item{para}{The parameters from \code{\link{partri}} or  \code{\link{vec2par}}.}
}
\value{
  Nonexceedance probability (\eqn{F}) for \eqn{x}.
}
\author{ W.H. Asquith}
\seealso{\code{\link{pdftri}}, \code{\link{quatri}}, \code{\link{lmomtri}}, \code{\link{partri}}}
\examples{
  lmr <- lmoms(c(46, 70, 59, 36, 71, 48, 46, 63, 35, 52))
  cdftri(50,partri(lmr))
}
\keyword{distribution}
\keyword{cumulative distribution function}
\keyword{Distribution: Asymmetric Triangular}
\keyword{Distribution: Triangular}
