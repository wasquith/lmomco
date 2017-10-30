\name{pdfln3}
\alias{pdfln3}
\title{Probability Density Function of the 3-Parameter Log-Normal Distribution}
\description{
This function computes the probability density
of the Log-Normal3 distribution given parameters (\eqn{\zeta}, lower bounds; \eqn{\mu_{\mathrm{log}}}, location; and \eqn{\sigma_{\mathrm{log}}}, scale)  computed by \code{\link{parln3}}. The probability density function function (same as Generalized Normal distribution, \code{\link{pdfgno}}) is
\deqn{f(x) = \frac{\exp(\kappa Y - Y^2/2)}{\alpha \sqrt{2\pi}} \mbox{,} }
where \eqn{Y} is
\deqn{
Y = \frac{\log(x - \zeta) - \mu_{\mathrm{log}}}{\sigma_{\mathrm{log}}}\mbox{,}
}
where \eqn{\zeta} is the lower bounds (real space) for which \eqn{\zeta < \lambda_1 - \lambda_2} (checked in \code{\link{are.parln3.valid}}), \eqn{\mu_{\mathrm{log}}} be the mean in natural logarithmic space, and \eqn{\sigma_{\mathrm{log}}} be the standard deviation in natural logarithm space for which \eqn{\sigma_{\mathrm{log}} > 0} (checked in \code{\link{are.parln3.valid}}) is obvious because this parameter has an analogy to the second product moment. Letting \eqn{\eta = \exp(\mu_{\mathrm{log}})}, the parameters of the Generalized Normal are \eqn{\zeta + \eta}, \eqn{\alpha = \eta\sigma_{\mathrm{log}}}, and \eqn{\kappa = -\sigma_{\mathrm{log}}}. At this point, the algorithms (\code{\link{pdfgno}}) for the Generalized Normal provide the functional core.
}
\usage{
pdfln3(x, para)
}
\arguments{
  \item{x}{A real value vector.}
  \item{para}{The parameters from \code{\link{parln3}} or \code{\link{vec2par}}.}
}
\value{
  Probability density (\eqn{f}) for \eqn{x}.
}
\note{
The parameterization of the Log-Normal3 results in ready support for either a known or unknown lower bounds. Details regarding the parameter fitting and control of the \eqn{\zeta} parameter can be seen under the Details section in \code{\link{parln3}}.
}
\author{ W.H. Asquith}
\references{
Asquith, W.H., 2011, Distributional analysis with L-moment statistics using the R environment for statistical computing: Createspace Independent Publishing Platform, ISBN 978--146350841--8.
}
\seealso{\code{\link{cdfln3}}, \code{\link{qualn3}}, \code{\link{lmomln3}}, \code{\link{parln3}}, \code{\link{pdfgno}}}
\examples{
  lmr <- lmoms(c(123,34,4,654,37,78))
  ln3 <- parln3(lmr); gno <- pargno(lmr)
  x <- qualn3(0.5,ln3)
  pdfln3(x,ln3) # 0.008053616
  pdfgno(x,gno) # 0.008053616 (the distributions are the same, but see Note)
}
\keyword{distribution}
\keyword{probability density function}
\keyword{Distribution: Log-Normal (3-parameter)}
