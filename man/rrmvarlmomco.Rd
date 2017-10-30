\name{rrmvarlmomco}
\alias{rrmvarlmomco}
\title{Reversed Variance Residual Quantile Function of the Distributions}
\description{
This function computes the Reversed Variance Residual Quantile Function for a quantile function \eqn{x{F}} (\code{\link{par2qua}}, \code{\link{qlmomco}}). The variance is defined by Nair et al. (2013, p. 58) as
\deqn{D(u) = \frac{1}{u} \int_0^u R(u)^2\; \mathrm{d}p\mbox{,}}
where \eqn{D(u)} is the variance of \eqn{R(u)}  (the reversed mean residual quantile function, \code{\link{rrmlmomco}}) for nonexceedance probability \eqn{u}.  The variance of \eqn{M(u)} is provided in \code{\link{rmvarlmomco}}.
}
\usage{
rrmvarlmomco(f, para)
}
\arguments{
  \item{f}{Nonexceedance probability (\eqn{0 \le F \le 1}).}
  \item{para}{The parameters from \code{\link{lmom2par}} or \code{\link{vec2par}}.}
}
\value{
  Reversed residual variance value for \eqn{F}.
}
\author{ W.H. Asquith}
\references{
Nair, N.U., Sankaran, P.G., and Balakrishnan, N., 2013, Quantile-based reliability analysis: Springer, New York.
}
\seealso{\code{\link{qlmomco}}, \code{\link{rrmlmomco}}
}
\examples{
# It is easiest to think about residual life as starting at the origin, units in days.
A <- vec2par(c(0.0, 264, 1.6), type="gov") # so set lower bounds = 0.0
rrmvarlmomco(0.5, A) # variance at the median reversed mean residual life
\dontrun{
A <- vec2par(c(-100, 264, 1.6), type="gov")
F <- nonexceeds(f01=TRUE)
plot(F, rmvarlmomco(F,A), type="l")
lines(F, rrmvarlmomco(F,A), col=2)
}
}
\keyword{quantile function}
\keyword{reversed variance residual quantile function}
\keyword{lifetime/reliability analysis}
\concept{lifetime}
\concept{lifetime}
\concept{reliability}
\keyword{The lmomco functions}
