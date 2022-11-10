\name{theopwms}
\alias{theopwms}
\title{The Theoretical Probability-Weighted Moments using Integration of the Quantile Function }
\description{
Compute the theoretrical probability-weighted moments (PWMs) for a distribution. A theoretrical PWM in integral form is
\deqn{ \beta_r = \int^1_0 x(F)\,F^r\,\mathrm{d}F \mbox{,}}
where \eqn{x(F)} is the quantile function of the random variable \eqn{X} for nonexceedance probability \eqn{F} and \eqn{r} represents the order of the PWM. This function loops across the above equation for each \code{nmom} set in the argument list. The function \eqn{x(F)} is computed through the \code{\link{par2qua}} function. The distribution type is determined using the \code{type} attribute of the \code{para} argument, which is a parameter object of \pkg{lmomco} (see \code{\link{vec2par}}).
}
\usage{
theopwms(para, nmom=5, verbose=FALSE)
}
\arguments{
  \item{para}{A distribution parameter object such as that by \code{\link{lmom2par}} or \code{\link{vec2par}}.}
  \item{nmom}{The number of moments to compute. Default is 5.}
  \item{verbose}{Toggle verbose output. Because the \R function \code{integrate} is used to perform the numerical integration, it might be useful to see selected messages regarding the numerical integration.}
}
\value{
  An \R \code{list} is returned.

  \item{betas}{The PWMs. Note that convention is the have a \eqn{\beta_0}, but this is placed in the first index \code{i=1} of the \code{betas} vector.}
  \item{source}{An attribute identifying the computational source of the probability-weighted moments: \dQuote{theopwms}.}
}
\references{
Hosking, J.R.M., 1990, L-moments--Analysis and estimation of distributions using linear combinations of order statistics: Journal of the Royal Statistical Society, Series B, v. 52, p. 105--124.
}
\author{W.H. Asquith}

\seealso{\code{\link{theoLmoms}}, \code{\link{pwm}}, \code{\link{pwm2lmom}}}
\examples{
para     <- vec2par(c(0,1),type='nor') # standard normal
the.pwms <- theopwms(para) # compute PWMs
str(the.pwms)
}
\keyword{probability-weighted moment (theoretical)}

