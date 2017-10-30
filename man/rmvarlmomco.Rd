\name{rmvarlmomco}
\alias{rmvarlmomco}
\title{Variance Residual Quantile Function of the Distributions}
\description{
This function computes the Variance Residual Quantile Function for a quantile function \eqn{x(F)} (\code{\link{par2qua}}, \code{\link{qlmomco}}).  The variance is defined by Nair et al. (2013, p. 55) as
\deqn{V(u) = \frac{1}{1-u} \int_u^1 M(u)^2\; \mathrm{d}p\mbox{,}}
where \eqn{V(u)} is the variance of \eqn{M(u)} (the residual mean quantile function, \code{\link{rmlmomco}}) for nonexceedance probability \eqn{u}.
}
\usage{
rmvarlmomco(f, para)
}
\arguments{
  \item{f}{Nonexceedance probability (\eqn{0 \le F \le 1}).}
  \item{para}{The parameters from \code{\link{lmom2par}} or \code{\link{vec2par}}.}
}
\value{
  Residual variance value for \eqn{F}.
}
\author{ W.H. Asquith}
\references{
Nair, N.U., Sankaran, P.G., and Balakrishnan, N., 2013, Quantile-based reliability analysis: Springer, New York.
}
\seealso{\code{\link{qlmomco}}, \code{\link{rmlmomco}}
}
\examples{
# It is easiest to think about residual life as starting at the origin, units in days.
A <- vec2par(c(0.0, 2649, 2.11), type="gov") # so set lower bounds = 0.0
qlmomco(0.5, A)  # The median lifetime = 1261 days
rmlmomco(0.5, A) # The average remaining life given survival to the median = 861 days
rmvarlmomco(0.5, A) # and the variance of that value.
\dontrun{
A <- lmom2par(vec2lmom(c(2000, 450, 0.14, 0.1)), type="kap")
F <- nonexceeds(f01=TRUE)
plot(F, qlmomco(F,A), type="l", ylim=c(100,6000),
     xlab="NONEXCEEDANCE PROBABILITY", ylab="LIFETIME OR SQRT(VAR LIFE), IN DAYS")
lines(F, sqrt( rmvarlmomco(F, A)), col=4, lwd=4) # thick blue, residual mean life
lines(F, sqrt(rrmvarlmomco(F, A)), col=2, lwd=4) # thick red, reversed resid. mean life
lines(F,   rmlmomco(F,A), col=4, lty=2); lines(F, rrmlmomco(F,A), col=2, lty=2)
lines(F,  tttlmomco(F,A), col=3, lty=2); lines(F,  cmlmomco(F,A), col=3)
}
}
\keyword{quantile function}
\keyword{variance residual quantile function}
\keyword{lifetime/reliability analysis}
\concept{lifetime}
\concept{reliability}
\keyword{The lmomco functions}
