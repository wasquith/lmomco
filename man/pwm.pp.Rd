\name{pwm.pp}
\alias{pwm.pp}
\title{ Plotting-Position Sample Probability-Weighted Moments }
\description{
The sample probability-weighted moments (PWMs) are computed from the plotting positions of the data. The first five \eqn{\beta_r}'s are computed by default. The plotting-position formula for a sample size of \eqn{n} is
\deqn{pp_i = \frac{i+A}{n+B} \mbox{,}}
where \eqn{pp_i} is the nonexceedance probability \eqn{F} of the \eqn{i}th ascending data values. An alternative form of the plotting position equation is
\deqn{pp_i = \frac{i + a}{n + 1 - 2a}\mbox{,}}
where \eqn{a} is a single plotting position coefficient. Having \eqn{a} provides \eqn{A} and \eqn{B}, therefore the parameters \eqn{A} and \eqn{B} together specify the plotting-position type. The PWMs are computed by
\deqn{\beta_r = n^{-1}\sum_{i=1}^{n}pp_i^r \times x_{j:n} \mbox{,}}
where \eqn{x_{j:n}} is the \eqn{j}th order statistic \eqn{x_{1:n} \le x_{2:n} \le x_{j:n} \dots \le x_{n:n}} of random variable X, and \eqn{r} is \eqn{0, 1, 2, \dots} for the PWM order.
}
\usage{
pwm.pp(x, pp=NULL, A=NULL, B=NULL, a=0, nmom=5, sort=TRUE)
}
\arguments{
  \item{x}{A vector of data values.}
  \item{pp}{An optional vector of nonexceedance probabilities. If present then \code{A} and \code{B} or \code{a} are ignored.}
  \item{A}{A value for the plotting-position formula. If \code{A} and \code{B} are both zero then the unbiased PWMs are computed through \code{\link{pwm.ub}}.}
  \item{B}{Another value for the plotting-position formula. If \code{A} and \code{B} are both zero then the unbiased PWMs are computed through \code{\link{pwm.ub}}.}
  \item{a}{A single plotting position coefficient from which, if not \code{NULL}, \eqn{A} and \eqn{B} will be internally computed;}
  \item{nmom}{Number of PWMs to return.}
  \item{sort}{Do the data need sorting? The computations require sorted data. This option is provided to optimize processing speed if presorted data already exists.}
}
\value{
  An \R \code{list} is returned.

  \item{betas}{The PWMs. Note that convention is the have a \eqn{\beta_0}, but this is placed in the first index \code{i=1} of the \code{betas} vector.}
  \item{source}{Source of the PWMs: \dQuote{pwm.pp}.}
}
\references{
Greenwood, J.A., Landwehr, J.M., Matalas, N.C., and Wallis, J.R., 1979, Probability weighted moments---Definition and relation to parameters of several distributions expressable in inverse form: Water Resources Research, v. 15, pp. 1,049--1,054.

Hosking, J.R.M., 1990, L-moments---Analysis and estimation of
distributions using linear combinations of order statistics: Journal
of the Royal Statistical Society, Series B, v. 52, pp. 105--124.
}
\author{ W.H. Asquith}
\seealso{ \code{\link{pwm.ub}}, \code{\link{pwm.gev}}, \code{\link{pwm2lmom}} }
\examples{
pwm <- pwm.pp(rnorm(20), A=-0.35, B=0)

X <- rnorm(20)
pwm <- pwm.pp(X, pp=pp(X)) # weibull plotting positions
}
\keyword{probability-weighted moment (sample)}
