\name{qualn3}
\alias{qualn3}
\title{Quantile Function of the 3-Parameter Log-Normal Distribution}
\description{
This function computes the quantiles of the Log-Normal3 distribution given parameters (\eqn{\zeta}, lower bounds; \eqn{\mu_{\mathrm{log}}}, location; and \eqn{\sigma_{\mathrm{log}}}, scale) of the distribution computed by \code{\link{parln3}}.  The quantile function (same as Generalized Normal distribution, \code{\link{quagno}}) is
\deqn{x = \Phi^{(-1)}(Y) \mbox{,} }
where \eqn{\Phi^{(-1)}} is the quantile function of the Standard Normal distribution and \eqn{Y} is
\deqn{
Y = \frac{\log(x - \zeta) - \mu_{\mathrm{log}}}{\sigma_{\mathrm{log}}}\mbox{,}
}
where  \eqn{\zeta} is the lower bounds (real space) for which \eqn{\zeta < \lambda_1 - \lambda_2} (checked in \code{\link{are.parln3.valid}}), \eqn{\mu_{\mathrm{log}}} be the mean in natural logarithmic space, and \eqn{\sigma_{\mathrm{log}}} be the standard deviation in natural logarithm space for which \eqn{\sigma_{\mathrm{log}} > 0} (checked in \code{\link{are.parln3.valid}}) is obvious because this parameter has an analogy to the second product moment. Letting \eqn{\eta = \exp(\mu_{\mathrm{log}})}, the parameters of the Generalized Normal are \eqn{\zeta + \eta}, \eqn{\alpha = \eta\sigma_{\mathrm{log}}}, and \eqn{\kappa = -\sigma_{\mathrm{log}}}. At this point, the algorithms (\code{\link{quagno}}) for the Generalized Normal provide the functional core.
}
\usage{
qualn3(f, para, paracheck=TRUE)
}
\arguments{
  \item{f}{Nonexceedance probability (\eqn{0 \le F \le 1}).}
  \item{para}{The parameters from \code{\link{parln3}} or  \code{\link{vec2par}}.}
  \item{paracheck}{A logical controlling whether the parameters are checked for validity. Overriding of this check might be extremely important and needed for use of the distribution quantile function in the context of TL-moments with nonzero trimming.}
}
\value{
  Quantile value for nonexceedance probability \eqn{F}.
}
\note{
The parameterization of the Log-Normal3 results in ready support for either a known or unknown lower bounds. More information regarding the parameter fitting and control of the \eqn{\zeta} parameter can be seen in the Details section under \code{\link{parln3}}.
}
\author{W.H. Asquith}
\references{
Asquith, W.H., 2011, Distributional analysis with L-moment statistics using the R environment for statistical computing: Createspace Independent Publishing Platform, ISBN 978--146350841--8.
}
\seealso{\code{\link{cdfln3}}, \code{\link{pdfln3}}, \code{\link{lmomln3}}, \code{\link{parln3}}, \code{\link{quagno}}}
\examples{
  lmr <- lmoms(c(123,34,4,654,37,78))
  qualn3(0.5,parln3(lmr))
}
\keyword{distribution}
\keyword{quantile function}
\keyword{Distribution: Log-Normal (3-parameter)}
