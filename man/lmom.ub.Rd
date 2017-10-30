\name{lmom.ub}
\alias{lmom.ub}
\title{ Unbiased Sample L-moments by Direct Sample Estimators }
\description{
Unbiased sample L-moments are computed for a vector using the direct sample estimation method as opposed to the use of sample probability-weighted moments. The L-moments are the ordinary L-moments and not the trimmed L-moments (see \code{\link{TLmoms}}). The mean, L-scale, coefficient of L-variation (\eqn{\tau}, LCV, L-scale/mean), L-skew (\eqn{\tau_3}, TAU3, L3/L2), L-kurtosis (\eqn{\tau_4}, TAU4, L4/L2), and \eqn{\tau_5} (TAU5, L5/L2) are computed. In conventional nomenclature, the L-moments are
\deqn{ \hat{\lambda}_1 = \mbox{L1} = \mbox{mean, }}
\deqn{ \hat{\lambda}_2 = \mbox{L2} = \mbox{L-scale, }}
\deqn{ \hat{\lambda}_3 = \mbox{L3} = \mbox{third L-moment, }}
\deqn{ \hat{\lambda}_4 = \mbox{L4} = \mbox{fourth L-moment, and }}
\deqn{ \hat{\lambda}_5 = \mbox{L5} = \mbox{fifth L-moment. }}

The L-moment ratios are
\deqn{ \hat{\tau} = \mbox{LCV} = \lambda_2/\lambda_1 = \mbox{coefficient of L-variation, }}
\deqn{ \hat{\tau}_3 = \mbox{TAU3} = \lambda_3/\lambda_2 = \mbox{L-skew, }}
\deqn{ \hat{\tau}_4 = \mbox{TAU4} = \lambda_4/\lambda_2 = \mbox{L-kurtosis, and}}
\deqn{ \hat{\tau}_5 = \mbox{TAU5} = \lambda_5/\lambda_2 = \mbox{not named.}}

It is common amongst practitioners to lump the L-moment ratios into the general term \dQuote{L-moments} and remain inclusive of the L-moment ratios. For example, L-skew then is referred to as the 3rd L-moment when it technically is the 3rd L-moment ratio.  The first L-moment ratio has no definition; the \code{\link{lmoms}} function uses the \code{NA} of \R in its vector representation of the ratios.

The mathematical expression for sample L-moment computation is shown under \code{\link{TLmoms}}. The formula jointly handles sample L-moment computation and sample TL-moment computation.
}
\usage{
lmom.ub(x)
}
\arguments{
  \item{x}{A vector of data values.}
}
\details{
The L-moment ratios (\eqn{\tau}, \eqn{\tau_3}, \eqn{\tau_4}, and \eqn{\tau_5}) are the primary higher L-moments for application, such as for distribution parameter estimation. However, the actual L-moments (\eqn{\lambda_3}, \eqn{\lambda_4}, and \eqn{\lambda_5}) are also reported. The implementation of \code{\link{lmom.ub}} requires a minimum of five data points. If more or fewer L-moments are needed then use the function \code{\link{lmoms}}.
}
\value{
  An \R \code{list} is returned.

  \item{L1}{Arithmetic mean.}
  \item{L2}{L-scale---analogous to standard deviation (see also \code{\link{gini.mean.diff}}.}
  \item{LCV}{coefficient of L-variation---analogous to coe. of variation.}
  \item{TAU3}{The third L-moment ratio or L-skew---analogous to skew.}
  \item{TAU4}{The fourth L-moment ratio or L-kurtosis---analogous to kurtosis.}
  \item{TAU5}{The fifth L-moment ratio.}
  \item{L3}{The third L-moment.}
  \item{L4}{The fourth L-moment.}
  \item{L5}{The fifth L-moment.}
  \item{source}{An attribute identifying the computational source of the L-moments: \dQuote{lmom.ub}.}
}
\note{
The \code{\link{lmom.ub}} function was among the first functions written for \pkg{lmomco} and actually written before \pkg{lmomco} was initiated. The \code{ub} was to be contrasted with plotting-position-based estimation methods: \code{\link{pwm.pp}} \eqn{\rightarrow} \code{\link{pwm2lmom}}. Further, at the time of development the radical expansion of \pkg{lmomco} beyond the Hosking (1996) FORTRAN libraries was not anticipated.  The author now exclusively uses \code{\link{lmoms}} but the numerical results should be identical. The direct sample estimator algorithm by Wang (1996) is used in \code{\link{lmom.ub}} and a more generalized algorithm is associated with \code{\link{lmoms}}.
}
\source{
The Perl code base of W.H. Asquith
}
\references{
Asquith, W.H., 2011, Distributional analysis with L-moment statistics using the R environment for statistical computing: Createspace Independent Publishing Platform, ISBN 978--146350841--8.

Hosking, J.R.M., 1990, L-moments---Analysis and estimation of
distributions using linear combinations of order statistics: Journal of the Royal Statistical Society, Series B, v. 52, pp. 105--124.

Wang, Q.J., 1996, Direct sample estimators of L-moments: Water Resources Research, v. 32, no. 12., pp. 3617--3619.
}
\author{ W.H. Asquith}
\seealso{\code{\link{lmom2pwm}}, \code{\link{pwm.ub}}, \code{\link{pwm2lmom}}, \code{\link{lmoms}}, \code{\link{lmorph}}}
\examples{
lmr <- lmom.ub(c(123,34,4,654,37,78))
lmorph(lmr)
lmom.ub(rnorm(100))
}
\keyword{L-moment (sample)}

