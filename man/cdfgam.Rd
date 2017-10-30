\name{cdfgam}
\alias{cdfgam}
\title{Cumulative Distribution Function of the Gamma Distribution}
\description{
This function computes the cumulative probability or nonexceedance probability of the Gamma distribution given parameters (\eqn{\alpha} and \eqn{\beta}) computed by \code{\link{pargam}}.  The cumulative distribution function has no explicit form but is expressed as an integral:
\deqn{F(x) = \frac{\beta^{-\alpha}}{\Gamma(\alpha)}\int_0^x t^{\alpha - 1}
\exp(-t/\beta)\; \mbox{d}t \mbox{,}}
where \eqn{F(x)} is the nonexceedance probability for the quantile \eqn{x}, \eqn{\alpha} is a shape parameter, and \eqn{\beta} is a scale parameter.

Alternatively, a three-parameter version is available following the parameterization of the Generalized Gamma distribution used in the \pkg{gamlss.dist} package and is

\deqn{F(x) =\frac{\theta^\theta\, |\nu|}{\Gamma(\theta)}\int_0^x \frac{z^\theta}{x}\,\mathrm{exp}(-z\theta)\; \mbox{d}x \mbox{,}}

where \eqn{z =(x/\mu)^\nu}, \eqn{\theta = 1/(\sigma^2\,|\nu|^2)} for \eqn{x > 0}, location parameter \eqn{\mu > 0}, scale parameter \eqn{\sigma > 0}, and shape parameter \eqn{-\infty < \nu < \infty}. The three parameter version is automatically triggered if the length of the \code{para} element is three and not two.
}
\usage{
cdfgam(x, para)
}
\arguments{
  \item{x}{A real value vector.}
  \item{para}{The parameters from \code{\link{pargam}} or  \code{\link{vec2par}}.}
}
\value{
  Nonexceedance probability (\eqn{F}) for \eqn{x}.
}
\references{
Hosking, J.R.M., 1990, L-moments---Analysis and estimation of
distributions using linear combinations of order statistics: Journal
of the Royal Statistical Society, Series B, v. 52, pp. 105--124.

Hosking, J.R.M., and Wallis, J.R., 1997, Regional frequency analysis---An
approach based on L-moments: Cambridge University Press.
}
\author{ W.H. Asquith}
\seealso{\code{\link{pdfgam}}, \code{\link{quagam}}, \code{\link{lmomgam}}, \code{\link{pargam}}}
\examples{
  lmr <- lmoms(c(123,34,4,654,37,78))
  cdfgam(50,pargam(lmr))

  # A manual demonstration of a gamma parent
  G  <- vec2par(c(0.6333,1.579),type='gam') # the parent
  F1 <- 0.25         # nonexceedance probability
  x  <- quagam(F1,G) # the lower quartile (F=0.25)
  a  <- 0.6333       # gamma parameter
  b  <- 1.579        # gamma parameter
  # compute the integral
  xf <- function(t,A,B) { t^(A-1)*exp(-t/B) }
  Q  <- integrate(xf,0,x,A=a,B=b)
  # finish the math
  F2 <- Q$val*b^(-a)/gamma(a)
  # check the result
  if(abs(F1-F2) < 1e-8) print("yes")

\dontrun{
# 3-p Generalized Gamma Distribution and gamlss.dist package parameterization
gg <- vec2par(c(7.4, 0.2, 14), type="gam"); X <- seq(0.04,9, by=.01)
GGa <- gamlss.dist::pGG(X, mu=7.4, sigma=0.2, nu=14)
GGb <- cdfgam(X, gg) # lets compare the two cumulative probabilities
plot( X, GGa, type="l", xlab="X", ylab="PROBABILITY", col=3, lwd=6)
lines(X, GGb, col=2, lwd=2) #}

\dontrun{
# 3-p Generalized Gamma Distribution and gamlss.dist package parameterization
gg <- vec2par(c(4, 1.5, -.6), type="gam"); X <- seq(0,1000, by=1)
GGa <- 1-gamlss.dist::pGG(X, mu=4, sigma=1.5, nu=-.6) # Note 1-... (pGG bug?)
GGb <- cdfgam(X, gg) # lets compare the two cumulative probabilities
plot( X, GGa, type="l", xlab="X", ylab="PROBABILITY", col=3, lwd=6)
lines(X, GGb, col=2, lwd=2) #}
}
\keyword{distribution}
\keyword{cumulative distribution function}
\keyword{Distribution: Gamma}
\keyword{Distribution: Generalized Gamma}
