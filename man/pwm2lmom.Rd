\name{pwm2lmom}
\alias{pwm2lmom}
\title{Probability-Weighted Moments to L-moments}
\description{
Converts the probability-weighted moments (PWM) to the L-moments. The conversion is linear so procedures based on PWMs are identical to those based on L-moments through a system of linear equations
\deqn{\lambda_1 = \beta_0 \mbox{,}}
\deqn{\lambda_2 = 2\beta_1 - \beta_0 \mbox{,}}
\deqn{\lambda_3 = 6\beta_2 - 6\beta_1 + \beta_0 \mbox{,}}
\deqn{\lambda_4 = 20\beta_3 - 30\beta_2 + 12\beta_1 - \beta_0 \mbox{,}}
\deqn{\lambda_5 = 70\beta_4 - 140\beta_3 + 90\beta_2 - 20\beta_1 + \beta_0 \mbox{,}}
\deqn{\tau = \lambda_2/\lambda_1 \mbox{,}}
\deqn{\tau_3 = \lambda_3/\lambda_2 \mbox{,}}
\deqn{\tau_4 = \lambda_4/\lambda_2 \mbox{, and}}
\deqn{\tau_5 = \lambda_5/\lambda_2 \mbox{.}}

The general expression and the expression used for computation if the argument is a vector of PWMs is
\deqn{\lambda_{r+1} = \sum^r_{k=0} (-1)^{r-k}{r \choose k}{r+k \choose k} \beta_{k+1}\mbox{.}}
}
\usage{
pwm2lmom(pwm)
}
\arguments{
  \item{pwm}{ A PWM object created by \code{\link{pwm.ub}} or similar. }
}
\details{
The probability-weighted moments (PWMs) are linear combinations of the L-moments and therefore contain the same statistical information of the data as the L-moments. However, the PWMs are harder to interpret as measures of probability distributions. The linearity between L-moments and PWMs means that procedures base on one are equivalent to the other.

The function can take a variety of PWM argument types in \code{\link{pwm}}. The function checks whether the argument is an \R \code{list} and if so attempts to extract the \eqn{\beta_r}'s from \code{list} names such as \code{BETA0}, \code{BETA1}, and so on. If the extraction is successful, then a list of L-moments similar to \code{\link{lmom.ub}} is returned. If the extraction was not successful, then an \R \code{list} name \code{betas} is checked; if \code{betas} is found, then this vector of PWMs is used to compute the L-moments. If \code{pwm} is a \code{list} but can not be routed in the function, a \code{warning} is made and \code{NULL} is returned. If the \code{pwm} argument is a \code{vector}, then this vector of PWMs is used. to compute the L-moments are returned.
}
\value{
  One of two \R \code{list}s are returned. Version I is

  \item{L1}{Arithmetic mean.}
  \item{L2}{L-scale---analogous to standard deviation.}
  \item{LCV}{coefficient of L-variation---analogous to coe. of variation.}
  \item{TAU3}{The third L-moment ratio or L-skew---analogous to skew.}
  \item{TAU4}{The fourth L-moment ratio or L-kurtosis---analogous to kurtosis.}
  \item{TAU5}{The fifth L-moment ratio.}
  \item{L3}{The third L-moment.}
  \item{L4}{The fourth L-moment.}
  \item{L5}{The fifth L-moment.}


Version II is 

  \item{lambdas}{The L-moments.}
  \item{ratios}{The L-moment ratios.}
  \item{source}{Source of the L-moments \dQuote{pwm2lmom}.}
}
\references{
Greenwood, J.A., Landwehr, J.M., Matalas, N.C., and Wallis, J.R., 1979, Probability weighted moments---Definition and relation to parameters of several distributions expressable in inverse form: Water Resources Research, v. 15, pp. 1,049--1,054.

Hosking, J.R.M., 1990, L-moments--Analysis and estimation of distributions using linear combinations of order statistics: Journal of the Royal Statistical Society, Series B, v. 52, pp. 105--124.
}
\author{ W.H. Asquith}
\seealso{\code{\link{lmom.ub}}, \code{\link{pwm.ub}}, \code{\link{pwm}}, \code{\link{lmom2pwm}}}
\examples{
D <- c(123,34,4,654,37,78)
pwm2lmom(pwm.ub(D))
pwm2lmom(pwm(D))
pwm2lmom(pwm(rnorm(100)))
}
\keyword{L-moment (theoretical)}
\keyword{probability-weighted moment (theoretical)}
