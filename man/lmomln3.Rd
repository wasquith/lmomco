\name{lmomln3}
\alias{lmomln3}
\title{L-moments of the 3-Parameter Log-Normal Distribution}
\description{
This function estimates the L-moments of the Log-Normal3 distribution given the parameters (\eqn{\zeta}, lower bounds; \eqn{\mu_{\mathrm{log}}}, location; and \eqn{\sigma_{\mathrm{log}}}, scale) from \code{\link{parln3}}.  The distribution is the same as the Generalized Normal with algebraic manipulation of the parameters, and \pkg{lmomco} does not have truly separate algorithms for the Log-Normal3 but uses those of the Generalized Normal. The discussion begins with the later distribution.

The two L-moments in terms of the Geneneralized Normal  distribution parameters (\code{\link{lmomgno}}) are
\deqn{\lambda_1 = \xi + \frac{\alpha}{\kappa}[1-\mathrm{exp}(\kappa^2/2)] \mbox{, and}}
\deqn{\lambda_2 = \frac{\alpha}{\kappa}(\mathrm{exp}(\kappa^2/2)(1-2\Phi(-\kappa/\sqrt{2})) \mbox{,}}
where \eqn{\Phi} is the cumulative distribution of the Standard Normal distribution. There are no simple expressions for \eqn{\tau_3}, \eqn{\tau_4}, and \eqn{\tau_5}, and numerical methods are used.

Let \eqn{\zeta} be the lower bounds (real space) for which \eqn{\zeta < \lambda_1 - \lambda_2} (checked in \code{\link{are.parln3.valid}}), \eqn{\mu_{\mathrm{log}}} be the mean in natural logarithmic space, and \eqn{\sigma_{\mathrm{log}}} be the standard deviation in natural logarithm space for which \eqn{\sigma_{\mathrm{log}} > 0} (checked in \code{\link{are.parln3.valid}}) is obvious because this parameter has an analogy to the second product moment. Letting \eqn{\eta = \exp(\mu_{\mathrm{log}})}, the parameters of the Generalized Normal are \eqn{\zeta + \eta}, \eqn{\alpha = \eta\sigma_{\mathrm{log}}}, and \eqn{\kappa = -\sigma_{\mathrm{log}}}. At this point the L-moments can be solved for using algorithms for the Generalized Normal.
}
\usage{
lmomln3(para)
}
\arguments{
  \item{para}{The parameters of the distribution.}
}
\value{
  An \R \code{list} is returned.

  \item{lambdas}{Vector of the L-moments. First element is
\eqn{\lambda_1}, second element is \eqn{\lambda_2}, and so on.}
  \item{ratios}{Vector of the L-moment ratios. Second element is
\eqn{\tau}, third element is \eqn{\tau_3} and so on.}
  \item{trim}{Level of symmetrical trimming used in the computation, which is \code{0}.}
  \item{leftrim}{Level of left-tail trimming used in the computation, which is \code{NULL}.}
  \item{rightrim}{Level of right-tail trimming used in the computation, which is \code{NULL}.}
  \item{source}{An attribute identifying the computational
source of the L-moments: \dQuote{lmomln3}.}
}
\author{ W.H. Asquith}
\references{
Asquith, W.H., 2011, Distributional analysis with L-moment statistics using the R environment for statistical computing: Createspace Independent Publishing Platform, ISBN 978--146350841--8.
}
\seealso{\code{\link{parln3}}, \code{\link{cdfln3}}, \code{\link{pdfln3}}, \code{\link{qualn3}}, \code{\link{lmomgno}} }
\examples{
X <- exp(rnorm(10))
pargno(lmoms(X))$para
parln3(lmoms(X))$para
}
\keyword{L-moment (distribution)}
\keyword{Distribution: Log-Normal (3-parameter)}
