\name{pwm.gev}
\alias{pwm.gev}
\title{ Generalized Extreme Value Plotting-Position Probability-Weighted Moments }
\description{
Generalized Extreme Value plotting-position probability-weighted moments (PWMs) are computed from a sample.  The first five \eqn{\beta_r}'s are computed by default. The plotting-position formula for the Generalized Extreme Value distribution is
\deqn{pp_i = \frac{i-0.35}{n} \mbox{,}}
where \eqn{pp_i} is the nonexceedance probability \eqn{F} of the \eqn{i}th ascending values of the sample of size \eqn{n}. The PWMs are computed by
\deqn{\beta_r = n^{-1}\sum_{i=1}^{n}pp_i^r \times x_{j:n} \mbox{,}}
where \eqn{x_{j:n}} is the \eqn{j}th order statistic
\eqn{x_{1:n} \le x_{2:n} \le x_{j:n} \dots \le x_{n:n}} of random variable X, and \eqn{r} is \eqn{0, 1, 2, \dots}. Finally, \code{\link{pwm.gev}} dispatches to \code{pwm.pp(data,A=-0.35,B=0)} and does not have its own logic.
}
\usage{
pwm.gev(x, nmom=5, sort=TRUE)
}
\arguments{
  \item{x}{A vector of data values.}
  \item{nmom}{Number of PWMs to return.}
  \item{sort}{Do the data need sorting? The computations require sorted data. This option is provided to optimize processing speed if presorted data already exists.}
}
\value{
  An \R \code{list} is returned.

  \item{betas}{The PWMs. Note that convention is the have a \eqn{\beta_0}, but this is placed in the first index \code{i=1} of the \code{betas} vector.}
  \item{source}{Source of the PWMs: \dQuote{pwm.gev}.}
}
\references{
Greenwood, J.A., Landwehr, J.M., Matalas, N.C., and Wallis, J.R., 1979, Probability weighted moments---Definition and relation to parameters of several distributions expressable in inverse form: Water Resources Research, v. 15, pp. 1,049--1,054.

Hosking, J.R.M., 1990, L-moments---Analysis and estimation of
distributions using linear combinations of order statistics: Journal
of the Royal Statistical Society, Series B, v. 52, pp. 105--124.
}
\author{ W.H. Asquith}
\seealso{ \code{\link{pwm.ub}}, \code{\link{pwm.pp}}, \code{\link{pwm2lmom}} }
\examples{
pwm <- pwm.gev(rnorm(20))
}
\keyword{probability-weighted moment (sample)}
\keyword{Distribution: Generalized Extreme Value}
