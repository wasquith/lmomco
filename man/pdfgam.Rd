\name{pdfgam}
\alias{pdfgam}
\title{Probability Density Function of the Gamma Distribution}
\description{
This function computes the probability density function
of the Gamma distribution given parameters (\eqn{\alpha}, shape, and \eqn{\beta}, scale)  computed by \code{\link{pargam}}.  The probability density function has no explicit form, but is expressed as an integral

\deqn{f(x|\alpha, \beta)^{\mathrm{lmomco}} = \frac{1}{\beta^\alpha\,\Gamma(\alpha)}\, x^{\alpha - 1}\, \mathrm{exp}(-x/\beta) \mbox{,}}

where \eqn{f(x)} is the probability density for the quantile \eqn{x}, \eqn{\alpha} is a shape parameter, and \eqn{\beta} is a scale parameter.

Alternatively, a three-parameter version is available for this package following the parameterization of the Generalized Gamma distribution used in the \pkg{gamlss.dist} package and is

\deqn{f(x|\mu,\sigma,\nu)_{\mathrm{gamlss.dist}}^{\mathrm{lmomco}}=\frac{\theta^\theta\, |\nu|}{\Gamma(\theta)}\,\frac{z^\theta}{x}\,\mathrm{exp}(-z\theta)\mbox{,}}

where \eqn{z =(x/\mu)^\nu}, \eqn{\theta = 1/(\sigma^2\,|\nu|^2)} for \eqn{x > 0}, location parameter \eqn{\mu > 0}, scale parameter \eqn{\sigma > 0}, and shape parameter \eqn{-\infty < \nu < \infty}. Note that for \eqn{\nu = 0} the distribution is log-Normal. The three parameter version is automatically triggered if the length of the \code{para} element is three and not two.
}
\usage{
pdfgam(x, para)
}
\arguments{
  \item{x}{A real value vector.}
  \item{para}{The parameters from \code{\link{pargam}} or \code{\link{vec2par}}.}
}
\value{
  Probability density (\eqn{f}) for \eqn{x}.
}
\note{
\bold{Two Parameter \eqn{\equiv} Three Parameter}\cr
For \eqn{\nu = 1}, the parameter conversion between the two gamma forms is \eqn{\alpha = \sigma^{-2}} and \eqn{\beta = \mu\sigma^2} and this can be readily verified:
\preformatted{
  mu <- 5; sig <- 0.7; nu <- 0
  x <- exp(seq(-3,3,by=.1))
  para2 <- vec2par(c(1/sig^2, (mu*sig^2)  ), type="gam")
  para3 <- vec2par(c(      mu,    sig,   1), type="gam")
  plot(x, pdfgam(x, para2), ylab="Gamma Density"); lines(x, pdfgam(x, para3))
}

\bold{Package flexsurv Generalized Gamma}\cr
The \pkg{flexsurv} package provides an \dQuote{original} (\code{GenGamma.orig}) and \dQuote{preferred} parameterization (\code{GenGamma}) of the Generalized Gamma distribution and discusses parameter conversion between the two. Here the parameterization of the preferred form is compared to that in \pkg{lmomco}. The probability density function of \code{dgengamma()} from \pkg{flexsurv} is

\deqn{f(x|\mu_2, \sigma_2, Q)_{\mathrm{flexsurv}} = \frac{\eta^\eta|Q|}{\sigma_2\Gamma(\eta)}\frac{1}{x}\, \mathrm{exp}\bigr\{\eta\times[wQ - \mathrm{exp}(wQ)]\bigr\}\mbox{,}}

where \eqn{\eta = Q^{-2}}, \eqn{w = \log(g/\eta)/Q} for \eqn{g \sim \mathrm{Gamma}(\eta, 1)} where \eqn{\mathrm{Gamma}} is the cumulative distribution function (presumably, need to verify this) of the Gamma distribution, and

\deqn{x \sim \mathrm{exp}(\mu_2 + w\sigma_2)\mbox{,}}

where \eqn{\mu_2 > 0}, \eqn{\sigma_2 > 0}, and \eqn{-\infty < Q < \infty}, and the log-Normal distribution results for \eqn{Q=0}.  These definitions for \pkg{flexsurv} seem incomplete to this author and further auditing is needed.

\bold{Additional Generalized Gamma Comparison}\cr
The default \pkg{gamlss.dist} package version uses so-called \emph{log.links} for \eqn{\mu} and \eqn{\sigma}, and so-called \emph{identity.link} for \eqn{\nu} and these links are implicit for \pkg{lmomco}. The parameters can be converted to \pkg{flexsurv} package equivalents by \eqn{\mu_2 = \log(\mu)}, \eqn{\sigma_2 = \sigma}, and \eqn{Q = \sigma\nu}, which is readily verified by
\preformatted{
  mu <- 2; sig <- 0.8; nu <- 0.2; x <- exp(seq(-3,1,by=0.1))
  para <- vec2par(c(mu,sig,nu), type="gam")
  dGG <- gamlss.dist::dGG(x, mu=mu, sigma=sig, nu=nu)
  plot( x, dGG, ylab="density", lwd=0.8, cex=2)
  lines(x, flexsurv::dgengamma(x, log(mu), sig, Q=sig*nu), col=8, lwd=5)
  lines(x, pdfgam(x, para), col=2)
}
What complicates the discussion further is that seemingly only the \emph{log.link} concept is manifested in the use of \eqn{\log(mu)} to provide the \eqn{\mu_2} for \code{flexsurv::}\code{dgengamma}.

\bold{On the Log-Normal via Generalized Gamma}\cr
The \pkg{gamlss.dist} package uses an \eqn{|\nu| < 1\mathrm{e{-}6}} trigger for the log-Normal calls. Further testing and the initial independent origin of \pkg{lmomco} code suggests that a primary trigger though can be based on the finiteness of the \code{lgamma(theta)} for \eqn{\theta}. This is used in \code{pdfgam} as well as \code{\link{cdfgam}} and \code{\link{quagam}}.
}
\references{
Hosking, J.R.M., 1990, L-moments---Analysis and estimation of
distributions using linear combinations of order statistics: Journal
of the Royal Statistical Society, Series B, v. 52, pp. 105--124.

Hosking, J.R.M. and Wallis, J.R., 1997, Regional frequency analysis---An
approach based on L-moments: Cambridge University Press.
}
\author{ W.H. Asquith}
\seealso{\code{\link{cdfgam}}, \code{\link{quagam}}, \code{\link{lmomgam}}, \code{\link{pargam}}}
\examples{
  lmr <- lmoms(c(123,34,4,654,37,78))
  gam <- pargam(lmr)
  x <- quagam(0.5,gam)
  pdfgam(x,gam)

\dontrun{
# 3-p Generalized Gamma Distribution and gamlss.dist package parameterization
gg <- vec2par(c(7.4, 0.2, 14), type="gam"); X <- seq(0.04,9, by=.01)
GGa <- gamlss.dist::dGG(X, mu=7.4, sigma=0.2, nu=14)
GGb <- pdfgam(X, gg) # We now compare the two densities.
plot( X, GGa, type="l", xlab="X", ylab="PROBABILITY DENSITY", col=3, lwd=6)
lines(X, GGb, col=2, lwd=2) #}

\dontrun{
# 3-p Generalized Gamma Distribution and gamlss.dist package parameterization
gg <- vec2par(c(1.7, 3, -4), type="gam"); X <- seq(0.04,9, by=.01)
GGa <- gamlss.dist::dGG(X, mu=1.7, sigma=3, nu=-4)
GGb <- pdfgam(X, gg) # We now compare the two densities.
plot( X, GGa, type="l", xlab="X", ylab="PROBABILITY DENSITY", col=3, lwd=6)
lines(X, GGb, col=2, lwd=2) #}
}
\keyword{distribution}
\keyword{probability density function}
\keyword{Distribution: Gamma}
\keyword{Distribution: Generalized Gamma}
