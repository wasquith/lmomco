\name{vec2pwm}
\alias{vec2pwm}
\title{Convert a Vector of Probability-Weighted Moments to a Probability-Weighted Moments Object}
\description{
This function converts a vector of probability-weighted moments (PWM) to a PWM
object of \pkg{lmomco}. The object is an \R \code{list}. This function is intended
to facilitate the use of PWM that the user might have from other sources. The
first five PWMs are supported (\eqn{\beta_0}, \eqn{\beta_1}, \eqn{\beta_2},
\eqn{\beta_3}, \eqn{\beta_4}) if \code{as.list=FALSE} otherwise the \eqn{\beta_r} are unlimited.
}
\usage{
vec2pwm(vec, as.list=FALSE)
}
\arguments{
  \item{vec}{A vector of PWM values in (\eqn{\beta_0}, \eqn{\beta_1}, \eqn{\beta_2}, \eqn{\beta_3}, \eqn{\beta_4}) order.}
  \item{as.list}{A logical controlling the returned data structure.}
}
\value{
  An \R \code{list} is returned if \code{as.list=TRUE}.

  \item{BETA0}{The first PWM, which is equal to the arithmetic mean.}
  \item{BETA1}{The second PWM.}
  \item{BETA2}{The third PWM.}
  \item{BETA3}{The fourth PWM.}
  \item{BETA4}{The fifth PWM.}
  \item{source}{Source of the PWMs: \dQuote{vec2pwm}.}

  Another \R \code{list} is returned if \code{as.list=FALSE}.

  \item{betas}{The PWMs.}
  \item{source}{Source of the PWMs: \dQuote{vec2pwm}.}
}
\author{W.H. Asquith}
\seealso{\code{\link{vec2lmom}}, \code{\link{lmom2pwm}}, \code{\link{pwm2lmom}}}
\examples{
pwm <- vec2pwm(c(12,123,12,12,54))
}
\keyword{utility (vector)}
