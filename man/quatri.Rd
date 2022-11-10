\name{quatri}
\alias{quatri}
\title{Quantile Function of the Asymmetric Triangular Distribution}
\description{
This function computes the quantiles of the Asymmetric Triangular distribution given parameters (\eqn{\nu}, \eqn{\omega}, and \eqn{\psi}) of the distribution computed by \code{\link{partri}}. The quantile function of the distribution is
\deqn{x(F) = \nu + \sqrt{(\psi - \nu)(\omega - \nu)F}\mbox{,}}
for \eqn{F < P},
\deqn{x(F) = \psi - \sqrt{(\psi - \nu)(\psi - \omega)(1-F)}\mbox{,}}
for \eqn{F > P}, and
\deqn{x(F) = \omega\mbox{,}}
for \eqn{F = P}
where \eqn{x(F)} is the quantile for nonexceedance probability \eqn{F}, \eqn{\nu} is the minimum, \eqn{\psi} is the maximum, and \eqn{\omega} is the mode of the distribution and
\deqn{P = \frac{(\omega - \nu)}{(\psi - \nu)}\mbox{.}}
}
\usage{
quatri(f, para, paracheck=TRUE)
}
\arguments{
  \item{f}{Nonexceedance probability (\eqn{0 \le F \le 1}).}
  \item{para}{The parameters from \code{\link{partri}} or \code{\link{vec2par}}.}
  \item{paracheck}{A logical controlling whether the parameters are checked for validity. Overriding of this check might be extremely important and needed for use of the quantile function in the context of TL-moments with nonzero trimming.}
}
\value{
  Quantile value for nonexceedance probability \eqn{F}.
}
\author{W.H. Asquith}
\seealso{\code{\link{cdftri}}, \code{\link{pdftri}}, \code{\link{lmomtri}}, \code{\link{partri}}}
\examples{
  lmr <- lmoms(c(46, 70, 59, 36, 71, 48, 46, 63, 35, 52))
  quatri(0.5,partri(lmr))
}
\keyword{distribution}
\keyword{quantile function}
\keyword{Distribution: Asymmetric Triangular}
\keyword{Distribution: Triangular}

