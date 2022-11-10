\name{pdfgld}
\alias{pdfgld}
\title{Probability Density Function of the Generalized Lambda Distribution}
\description{
This function computes the probability density function of the Generalized Lambda distribution given parameters (\eqn{\xi}, \eqn{\alpha}, \eqn{\kappa}, and \eqn{h})  computed by \code{\link{pargld}} or similar. The probability density function is
\deqn{f(x) = {[(\kappa[F(x)^{\kappa-1}] + h[1-F(x)])^{h-1})\times\alpha]}^{-1} \mbox{,}}
where \eqn{f(x)} is the probability density function at \eqn{x}, \eqn{F(x)} is the cumulative distribution function at \eqn{x}.
}
\usage{
pdfgld(x, para, paracheck)
}
\arguments{
  \item{x}{A real value vector.}
  \item{para}{The parameters from \code{\link{pargld}} or \code{\link{vec2par}}.}
  \item{paracheck}{A logical switch as to whether the validity of the parameters should be checked. Default is \code{paracheck=TRUE}. This switch is made so that the root solution needed for \code{\link{cdfgld}} exhibits an extreme speed increase because of the repeated calls to \code{\link{quagld}}.}
}
\value{
  Probability density (\eqn{f}) for \eqn{x}.
}
\references{
Asquith, W.H., 2007, L-moments and TL-moments of the generalized lambda distribution: Computational Statistics and Data Analysis, v. 51, no. 9, pp. 4484--4496.

Karian, Z.A., and Dudewicz, E.J., 2000, Fitting statistical distributions---The generalized lambda distribution and generalized bootstrap methods: CRC Press, Boca Raton, FL, 438 p.
}
\author{W.H. Asquith}
\seealso{\code{\link{cdfgld}}, \code{\link{quagld}}, \code{\link{lmomgld}}, \code{\link{pargld}}}
\examples{
\dontrun{
# Using Karian and Dudewicz, 2000, p. 10
gld <- vec2par(c(0.0305,1/1.3673,0.004581,0.01020),type='gld')
quagld(0.25,gld) # which equals about 0.028013 as reported by K&D
pdfgld(0.028013,gld) # which equals about 43.04 as reported by K&D
F <- seq(.001,.999,by=.001)
x <- quagld(F,gld)
plot(x, pdfgld(x,gld), type='l', xlim=c(0,.1))
}
}
\keyword{distribution}
\keyword{probability density function}
\keyword{Distribution: Generalized Lambda}
