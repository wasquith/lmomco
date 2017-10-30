\name{lmomaep4}
\alias{lmomaep4}
\title{L-moments of the 4-Parameter Asymmetric Exponential Power Distribution}
\description{
This function computes the L-moments of the 4-parameter Asymmetric Exponential Power distribution given the parameters (\eqn{\xi}, \eqn{\alpha}, \eqn{\kappa}, and \eqn{h}) from \code{\link{paraep4}}. The first four L-moments are complex. The mean \eqn{\lambda_1} is
\deqn{\lambda_1 = \xi + \alpha(1/\kappa - \kappa)\frac{\Gamma(2/h)}{\Gamma(1/h)}\mbox{,}}
where \eqn{\Gamma(x)} is the complete gamma function or \code{gamma()} in \R.

The L-scale \eqn{\lambda_2} is
\deqn{\lambda_2 = -\frac{\alpha\kappa(1/\kappa - \kappa)^2\Gamma(2/h)}
                        {(1+\kappa^2)\Gamma(1/h)}
                + 2\frac{\alpha\kappa^2(1/\kappa^3 + \kappa^3)\Gamma(2/h)I_{1/2}(1/h,2/h)}
                        {(1+\kappa^2)^2\Gamma(1/h)}\mbox{,}}
where \eqn{I_{1/2}(1/h,2/h)} is the cumulative distribution function of the Beta distribution (\eqn{I_x(a,b)}) or \code{pbeta(1/2,} \code{shape1=1/h,} \code{shape2=2/h)} in \R. This function is also referred to as the normalized incomplete beta function (Delicado and Goria, 2008) and defined as
\deqn{I_x(a,b) =  \frac{\int_0^x t^{a-1} (1-t)^{b-1}\; \mathrm{d}t}{\beta(a,b)}\mbox{,}}
where \eqn{\beta(1/h, 2/h)} is the complete beta function or \code{beta(1/h, 2/h)} in \R.

The third L-moment \eqn{\lambda_3} is
\deqn{\lambda_3 = A_1 + A_2 + A_3\mbox{,}}
where the \eqn{A_i} are
\deqn{A_1 = \frac{\alpha(1/\kappa - \kappa)(\kappa^4 - 4\kappa^2 + 1)\Gamma(2/h)}
                 {(1+\kappa^2)^2\Gamma(1/h)}\mbox{,}}
\deqn{A_2 = -6\frac{\alpha\kappa^3(1/\kappa - \kappa)(1/\kappa^3 + \kappa^3)\Gamma(2/h)I_{1/2}(1/h,2/h)}
                 {(1+\kappa^2)^3\Gamma(1/h)}\mbox{,}}
\deqn{A_3 = 6\frac{\alpha(1+\kappa^4)(1/\kappa - \kappa)\Gamma(2/h)\Delta}
                 {(1+\kappa^2)^2\Gamma(1/h)}\mbox{,}}
and where \eqn{\Delta} is
\deqn{\Delta = \frac{1}{\beta(1/h, 2/h)}\int_0^{1/2} t^{1/h - 1} (1-t)^{2/h - 1} I_{(1-t)/(2-t)}(1/h, 3/h) \; \mathrm{d}t\mbox{.}}


The fourth L-moment \eqn{\lambda_4} is
\deqn{\lambda_4 = B_1 + B_2 + B_3 + B_4\mbox{,}}
where the \eqn{B_i} are
\deqn{B_1 = -\frac{\alpha\kappa(1/\kappa - \kappa)^2(\kappa^4 - 8\kappa^2 + 1)\Gamma(2/h)}
                  {(1+\kappa^2)^3\Gamma(1/h)}\mbox{,}}
\deqn{B_2 = 12\frac{\alpha\kappa^2(\kappa^3 + 1/\kappa^3)(\kappa^4 - 3\kappa^2 + 1)\Gamma(2/h)I_{1/2}(1/h,2/h)}
                   {(1+\kappa^2)^4\Gamma(1/h)}\mbox{,}}
\deqn{B_3 = -30\frac{\alpha\kappa^3(1/\kappa - \kappa)^2(1/\kappa^2 + \kappa^2)\Gamma(2/h)\Delta}
                    {(1+\kappa^2)^3\Gamma(1/h)}\mbox{,}}
\deqn{B_4 = 20\frac{\alpha\kappa^4(1/\kappa^5 + \kappa^5)\Gamma(2/h)\Delta_1}
                   {(1+\kappa^2)^4\Gamma(1/h)}\mbox{,}}
and where \eqn{\Delta_1} is
\deqn{\Delta_1 =  \frac{\int_0^{1/2}  \int_0^{(1-y)/(2-y)}   y^{1/h - 1} (1-y)^{2/h - 1}
                       z^{1/h - 1} (1-z)^{3/h - 1}
        \;I'\; \mathrm{d}z\,\mathrm{d}y}{\beta(1/h, 2/h)\beta(1/h, 3/h)}\mbox{,}}
for which \eqn{I' = I_{(1-z)(1-y)/(1+(1-z)(1-y))}(1/h, 2/h)} is the cumulative distribution function of the beta distribution (\eqn{I_x(a,b)}) or \code{pbeta((1-z)(1-y)/(1+(1-z)(1-y)), shape1=1/h, shape2=2/h)} in \R.
}
\usage{
lmomaep4(para, paracheck=TRUE, t3t4only=FALSE)
}
\arguments{
  \item{para}{The parameters of the distribution.}
  \item{paracheck}{Should the parameters be checked for validity by the \code{\link{are.paraep4.valid}} function.}
  \item{t3t4only}{Return only the \eqn{\tau_3} and \eqn{\tau_4} for the parameters \eqn{\kappa} and \eqn{h}. The \eqn{\lambda_1} and \eqn{\lambda_2} are not explicitly used although numerical values for these two L-moments are required only to avoid computational errors. Care is made so that the \eqn{\alpha} parameter that is in numerator of \eqn{\lambda_{2,3,4}} is not used in the computation of \eqn{\tau_3} and \eqn{\tau_4}. Hence, this option permits the computation of \eqn{\tau_3} and \eqn{\tau_4} when \eqn{\alpha} is unknown. This features is largely available for research and development purposes. Mostly this feature was used for the \eqn{\{\tau_3, \tau_4\}} trajectory for \code{\link{lmrdia}}}.
}
\value{
  An \R \code{list} is returned.

  \item{lambdas}{Vector of the L-moments. First element is
\eqn{\lambda_1}, second element is \eqn{\lambda_2}, and so on.}
  \item{ratios}{Vector of the L-moment ratios. Second element is
\eqn{\tau}, third element is \eqn{\tau_3} and so on. }
  \item{trim}{Level of symmetrical trimming used in the computation, which is \code{0}.}
  \item{leftrim}{Level of left-tail trimming used in the computation, which is \code{NULL}.}
  \item{rightrim}{Level of right-tail trimming used in the computation, which is \code{NULL}.}
  \item{source}{An attribute identifying the computational
source of the L-moments: \dQuote{lmomaep4}.}

   or an alternative \R list is returned if \code{t3t4only=TRUE}

  \item{T3}{L-skew, \eqn{\tau_3}.}
  \item{T4}{L-kurtosis, \eqn{\tau_4}.}
}
\author{ W.H. Asquith}
\references{
Asquith, W.H., 2014, Parameter estimation for the 4-parameter asymmetric exponential power distribution by the method of L-moments using R: Computational Statistics and Data Analysis, v. 71, pp. 955--970.

Delicado, P., and Goria, M.N., 2008, A small sample comparison of maximum likelihood,
moments and L-moments methods for the asymmetric exponential power distribution:
Computational Statistics and Data Analysis, v. 52, no. 3, pp. 1661--1673.
}
\seealso{\code{\link{paraep4}}, \code{\link{cdfaep4}}, \code{\link{pdfaep4}}, \code{\link{quaaep4}} }
\examples{
\dontrun{
para <- vec2par(c(0, 1, 0.5, 4), type="aep4")
lmomaep4(para)
}
}
\keyword{L-moment (distribution)}
\keyword{Distribution: Asymmetric Exponential Power}
