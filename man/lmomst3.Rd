\name{lmomst3}
\alias{lmomst3}
\title{L-moments of the 3-Parameter Student t Distribution}
\description{
This function estimates the first six L-moments of the 3-parameter Student t distribution given the parameters (\eqn{\xi}, \eqn{\alpha}, \eqn{\nu}) from \code{\link{parst3}}. The L-moments in terms of the parameters are
\deqn{\lambda_1 = \xi\mbox{,}}
\deqn{\lambda_2 = 2^{6-4\nu}\pi\alpha\nu^{1/2}\,\Gamma(2\nu-2)/[\Gamma(\frac{1}{2}\nu)]^4\mbox{\, and}}
\deqn{\tau_4 = \frac{15}{2} \frac{\Gamma(\nu)}{\Gamma(\frac{1}{2})\Gamma(\nu - \frac{1}{2})} \int_0^1 \! \frac{(1-x)^{\nu - 3/2}[I_x(\frac{1}{2},\frac{1}{2}\nu)]^2}{\sqrt{x}}\; \mathrm{d} x - \frac{3}{2}\mbox{,}}
where \eqn{I_x(\frac{1}{2}, \frac{1}{2}\nu)} is the cumulative distribution function of the Beta distribution.  The distribution is symmetrical so that \eqn{\tau_r = 0} for odd values of \eqn{r: r \ge 3}.

Numerical integration of is made to estimate \eqn{\tau_4}. The other two parameters are readily solved for when \eqn{\nu} is available. A polynomial approximation is used to estimate the \eqn{\tau_6} as a function of \eqn{\tau_4}; the polynomial was based on the \code{\link{theoLmoms}} estimating \eqn{\tau_4} and \eqn{\tau_6}. The \eqn{\tau_6} polynomial has nine coefficients with a maximum absolute residual value of 2.065e-06 for 4,000 degrees of freedom (see \code{inst/doc/t4t6/studyST3.R}).
}
\usage{
lmomst3(para, ...)
}
\arguments{
  \item{para}{The parameters of the distribution.}
  \item{...}{Additional arguments to pass.}
}
\value{
  An \R \code{list} is returned.

  \item{lambdas}{Vector of the L-moments. First element is \eqn{\lambda_1}, second element is \eqn{\lambda_2}, and so on.}
  \item{ratios}{Vector of the L-moment ratios. Second element is \eqn{\tau}, third element is \eqn{\tau_3} and so on.}
  \item{trim}{Level of symmetrical trimming used in the computation, which is \code{0}.}
  \item{leftrim}{Level of left-tail trimming used in the computation, which is \code{NULL}.}
  \item{rightrim}{Level of right-tail trimming used in the computation, which is \code{NULL}.}
  \item{source}{An attribute identifying the computational source of the L-moments: \dQuote{lmomst3}.}
}
\author{W.H. Asquith with A.R. Biessen}
\references{
Asquith, W.H., 2011, Distributional analysis with L-moment statistics using the R environment for statistical computing: Createspace Independent Publishing Platform, ISBN 978--146350841--8.
}
\seealso{\code{\link{parst3}}, \code{\link{cdfst3}}, \code{\link{pdfst3}}, \code{\link{quast3}} }
\examples{
lmomst3(vec2par(c(1124, 12.123, 10), type="st3"))
}
\keyword{L-moment (distribution)}
\keyword{Distribution: Student t (3-parameter)}
