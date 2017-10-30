\name{lmomkap}
\alias{lmomkap}
\title{L-moments of the Kappa Distribution}
\description{
This function estimates the L-moments of the Kappa distribution given the parameters (\eqn{\xi}, \eqn{\alpha}, \eqn{\kappa}, and \eqn{h}) from \code{\link{parkap}}. The L-moments in terms of the parameters are complicated and are solved numerically. If the parameter \eqn{k = 0} (is small or near zero) then let
\deqn{d_r = \gamma + \log(-h) + \mathrm{digamma}(-r/h)\ \mbox{for}\ h < 0}
\deqn{d_r = \gamma + \log(r)\ \mbox{for}\ h = 0\ \mbox{(is small)}}
\deqn{d_r = \gamma + \log(h) + \mathrm{digamma}(1+r/h)\ \mbox{for}\ h > 0}
or if \eqn{k > -1} (nonzero) then let
\deqn{g_r = \frac{\Gamma(1+k)\Gamma(-r/h-k)}{-h^k\,\Gamma(-r/h)}\ \mbox{for}\ h < 0}
\deqn{g_r = \frac{\Gamma(1+k)}{r^k} \times (1-0.5hk(1+k)/r)\ \mbox{for}\ h = 0\ \mbox{(is small)}}
\deqn{g_r = \frac{\Gamma(1+k)\Gamma(1+r/h)}{h^g\,\Gamma(1+k+r/h)}\ \mbox{for}\ h > 0}
where \eqn{r} is L-moment order, \eqn{\gamma} is Euler's constant, and for \eqn{h = 0} the term to the right of the multiplication is not in Hosking (1994) or Hosking and Wallis (1997) for exists within Hosking's FORTRAN code base.

The probability-weighted moments (\eqn{\beta_r}; \code{\link{pwm2lmom}}) for \eqn{k = 0} (is small or near zero) are
\deqn{r\beta_{r-1} = \xi + (\alpha/\kappa)[1 - d_r]}
or if \eqn{k > -1} (nonzero) then 
\deqn{r\beta_{r-1} = \xi + (\alpha/\kappa)[1 - g_r]}

}
\usage{
lmomkap(para, nmom=5)
}
\arguments{
  \item{para}{The parameters of the distribution.}
  \item{nmom}{The number of moments to compute. Default is 5.}
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
source of the L-moments: \dQuote{lmomkap}.}
}
\references{
Hosking, J.R.M., 1994, The four-parameter kappa distribution: IBM Journal of Reserach and Development, v. 38, no. 3, pp. 251--258.

Hosking, J.R.M., and Wallis, J.R., 1997, Regional frequency analysis---An approach based on L-moments: Cambridge University Press.
}
\author{ W.H. Asquith}
\seealso{\code{\link{parkap}}, \code{\link{cdfkap}}, \code{\link{pdfkap}}, \code{\link{quakap}} }
\examples{
lmr <- lmoms(c(123,34,4,654,37,78))
lmomkap(parkap(lmr))
}
\keyword{L-moment (distribution)}
\keyword{Distribution: Kappa}
