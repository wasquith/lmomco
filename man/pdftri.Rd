\name{pdftri}
\alias{pdftri}
\title{Probability Density Function of the Asymmetric Triangular Distribution}
\description{
This function computes the probability density of the Asymmetric Triangular distribution given parameters (\eqn{\nu}, \eqn{\omega}, and \eqn{\psi})  computed by \code{\link{partri}}. The probability density function is
\deqn{f(x) = \frac{2(x-\nu)}{(\omega - \nu)(\psi - \nu)}\mbox{,}}
for \eqn{x < \omega},
\deqn{f(x) = \frac{2(\psi-x)}{(\psi - \omega)(\psi - \nu)}\mbox{,}}
for \eqn{x >  \omega}, and
\deqn{f(x) = \frac{2}{(\psi - \nu)}\mbox{,}}
for \eqn{x = \omega}
where \eqn{x(F)} is the quantile for nonexceedance probability \eqn{F}, \eqn{\nu} is the minimum, \eqn{\psi} is the maximum, and \eqn{\omega} is the mode of the distribution.
}
\usage{
pdftri(x, para)
}
\arguments{
  \item{x}{A real value vector.}
  \item{para}{The parameters from \code{\link{partri}} or \code{\link{vec2par}}.}
}

\value{
  Probability density (\eqn{f}) for \eqn{x}.
}
\author{ W.H. Asquith}
\seealso{\code{\link{pdftri}}, \code{\link{quatri}}, \code{\link{lmomtri}}, \code{\link{partri}}}
\examples{
  tri <- vec2par(c(-120, 102, 320), type="tri")
  x <- quatri(nonexceeds(),tri)
  pdftri(x,tri)
}
\keyword{distribution}
\keyword{probability density function}
\keyword{Distribution: Asymmetric Triangular}
\keyword{Distribution: Triangular}
