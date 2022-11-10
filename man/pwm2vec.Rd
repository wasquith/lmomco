\name{pwm2vec}
\alias{pwm2vec}
\title{Convert Probability-Weighted Moment object to a Vector}
\description{
This function converts a probability-weighted moment object in the structure used by \pkg{lmomco} into a simple vector of \eqn{\beta_0}, \eqn{\beta_1}, \eqn{\beta_2}, \eqn{\beta_3}, \eqn{\beta_4}, \dots, \eqn{\beta_{r-1}}.
}
\usage{
pwm2vec(pwm, ...)
}
\arguments{
  \item{pwm}{Probability-weighted moment object such as from \code{\link{pwm}} and \code{\link{vec2pwm}}.}
  \item{...}{Not presently used.}
}
\value{
A vector of the first five probability-weighted moments if available. The \code{$betas} field of the \code{pwm} argument is simply returned by this function.
}
\author{W.H. Asquith}
\seealso{\code{\link{pwm}}, \code{\link{vec2pwm}}, \code{\link{lmom2vec}}
}
\examples{
  pmr <- pwm(rnorm(40));             pwm2vec(pmr)
  pmr <- vec2pwm(c(140,150,45,21));  pwm2vec(pmr)
}
\keyword{utility (vector)}
