\name{pdfst3}
\alias{pdfst3}
\title{Probability Density Function of the 3-Parameter Student t Distribution}
\description{
This function computes the probability density of the 3-parameter Student t distribution given parameters (\eqn{\xi}, \eqn{\alpha}, \eqn{\nu})  computed by \code{\link{parst3}}.  The probability density function is
\deqn{
f(x) = \frac{\Gamma(\frac{1}{2} + \frac{1}{2}\nu)}{\alpha\nu^{1/2}\,\Gamma(\frac{1}{2})\Gamma(\frac{1}{2}\nu)}(1+t^2/\nu)^{-(\nu+1)/2}\mbox{,}
}
where \eqn{f(x)} is the probability density for quantile \eqn{x}, \eqn{t} is defined as \eqn{t = (x - \xi)/\alpha}, \eqn{\xi} is a location parameter, \eqn{\alpha} is a scale parameter, and \eqn{\nu} is a shape parameter in terms of the degrees of freedom as for the more familiar Student t distribution in \R.

 For value \code{X}, the built-in \R functions can be used. For \eqn{\nu \ge 1000}, one can use \code{dnorm(X, mean=U,} \code{sd=A)} and for \code{U} = \eqn{\xi} and \code{A}=\eqn{\alpha} for \eqn{1.000001 \le \nu \le 1000}, one can use \code{dt((X-U)/A, N)/A} for \code{N}=\eqn{\nu}. The \R function \code{dnorm} is used for the Normal distribution and the \R function \code{dt} is used for the 1-parameter Student t distribution.
}
\usage{
pdfst3(x, para, paracheck=TRUE)
}
\arguments{
  \item{x}{A real value vector.}
  \item{para}{The parameters from \code{\link{parst3}} or \code{\link{vec2par}}.}
  \item{paracheck}{A logical on whether the parameter should be check for validity.}
}
\value{
  Probability density (\eqn{f}) for \eqn{x}.
}
\author{ W.H. Asquith}
\references{
Asquith, W.H., 2011, Distributional analysis with L-moment statistics using the R environment for statistical computing: Createspace Independent Publishing Platform, ISBN 978--146350841--8.
}
\seealso{\code{\link{cdfst3}}, \code{\link{quast3}}, \code{\link{lmomst3}},  \code{\link{parst3}}}
\examples{
\dontrun{
xs <- -200:200
  para <- vec2par(c(37, 25,  114), type="st3")
plot(xs, pdfst3(xs, para), type="l")
  para <- vec2par(c(11, 36, 1000), type="st3")
lines(xs, pdfst3(xs, para), lty=2)
  para <- vec2par(c(-7, 60,   40), type="st3")
lines(xs, pdfst3(xs, para), lty=3)
}
}
\keyword{distribution}
\keyword{probability density function}
\keyword{Distribution: Student t (3-parameter)}
