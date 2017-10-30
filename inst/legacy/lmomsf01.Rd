\name{lmomsf01}
\alias{lmomsf01}
\title{ The Sample L-moments and L-moment Ratios from Nonexceedance Probabilities on the Bounds of 0 to 1 }
\description{
Compute the sample L-moments for a vector given another vector of nonexceedance probabilities \eqn{F} according to the following

\deqn{\lambda_r = \int_0^1 x(F) P^\star_{r-1}(F) \; \mathrm{d}F}

where \eqn{P^\star_{r}(F)} is

\deqn{P^\star_{r}(F)  = \sum_{k=0}^r (-1)^{r-k} {r \choose k} {r+k \choose k} F^k \mathrm{.}}

The quantile function \eqn{x(F)} is approximated by linear interpolation using the \code{approx} function with flat line (constant) extrapolation to the left and right of the end points of the input vector \code{x}. The extrapolation however can be controlled by functional arguments.
}
\usage{
lmomsf01(x, f=NULL, a=0, nmom=5, subdivisionscale=50,
       sort=TRUE, constantails=TRUE, efs=1e-3)
}
\arguments{
  \item{x}{A vector of data values.}
  \item{f}{A vector of nonexceedance probabilities for each value in \code{x}.}
  \item{a}{A plotting position coefficient as also used by the \code{\link{pp}} function for estimation of nonexceedance probabilities in \code{f} if \code{f=NULL}. The default is for the Weibull plotting positions.}
  \item{nmom}{The number of moments to compute. Default is 5.}
  \item{subdivisionscale}{A multiplicative factor on the sample size (\code{length(x)}) for the \code{subdivision} argument of the \code{integrate} function. The default value is not optimal but determined by some limited trial and error by the author.}
  \item{sort}{Does the data need sorting? The computations require sorted data. This option is provided to optimize processing speed if presorted data already exists.}
  \item{constantails}{During the linear interpolation, is \code{rule=2} of the \code{approx} function to be used if \code{constantails=TRUE} or will tail extrapolation using the last left and right tail slopes of the \eqn{x} and \eqn{F} values be used.}
  \item{efs}{A near zero value for a difference between the pair of left and right edge \eqn{F} values that are used for computation of the demoninator of a slope computation. The default value is relatively ad hoc determined but provides protection from slopes becoming too large. This argument only has effect if \code{constantails=FALSE}.}
}
\value{
  An \R \code{list} is returned.

  \item{lambdas}{Vector of the L-moments. First element is
\eqn{\hat{\lambda}^{(0,0)}_1}, second element is \eqn{\hat{\lambda}^{(0,0)}_2}, and so on.}
  \item{ratios}{Vector of the L-moment ratios. Second element is
\eqn{\hat{\tau}^{(0,0)}}, third element is \eqn{\hat{\tau}^{(0,0)}_3} and so on.}
  \item{trim}{Level of symmetrical trimming used in the computation, which will equal \code{NULL} if asymmetrical trimming was used.}
  \item{leftrim}{Level of left-tail trimming used in the computation.}
  \item{rightrim}{Level of right-tail trimming used in the computation.}
  \item{source}{An attribute identifying the computational source of the L-moments: \dQuote{lmomsf}.}
}
\note{
Trimming of the L-moments is not yet supported.
}
\references{
Asquith, W.H., 2011, Distributional analysis with L-moment statistics using the R environment for statistical computing: CreateSpace, ISBN 978--1--46350841--8, 344 p. [eqs. 6.2 and 6.3, p. 91]
}
\author{ W.H. Asquith }
\seealso{ \code{\link{lmoms}} }
\examples{

\dontrun{
   set.seed(11)
   x <- rnorm(35)
   lmomsf01(x)$lambdas # linear approx via weibull
   lmomsf01(sort(x), f=pp(x))$lambdas # same
   lmomsf01(x, f=pp(x, sort=FALSE))$lambdas # same
   pwm2lmom(pwm.pp(x, pp=pp(x, sort=FALSE)))$lambdas # directly by weibull
   pwm2lmom(pwm.pp(x, B=1))$lambdas # same
   lmoms(x)$lambdas # unbiased estimates

  lmr0 <- lmoms(x)
  lmr1 <- lmomsf01(sort(x), f=pp(x))
  lmr2 <- pwm2lmom(pwm.pp(x, B=1))
  Fs <- seq(0.001, 0.999, by=0.001)
  F <- pp(x, sort=FALSE)
  xofF <- approx(x=F, y=x, xout=Fs, rule=2)$y
  plot(qnorm(Fs), xofF, type="l",
       xlab="STANDARD NORMAL VARIATE",
       ylab="VALUE")
  points(qnorm(F), x)
  lines(qnorm(Fs), qlmomco(Fs, parnor(lmr0)), col=2)
  lines(qnorm(Fs), qlmomco(Fs, parnor(lmr1)))
  lines(qnorm(Fs), qlmomco(Fs, parnor(lmr2)), col=3)
}
}
\keyword{L-moment (sample)}
