\name{ralmomco}
\alias{ralmomco}
\title{Alpha-Percentile Residual Quantile Function of the Distributions}
\description{
This function computes the \eqn{\alpha}-Percentile Residual Quantile Function for quantile function \eqn{x(F)} (\code{\link{par2qua}}, \code{\link{qlmomco}}). The function is defined by Nair and Vineshkumar (2011, p. 85) and Nair et al. (2013, p. 56) as
\deqn{P_\alpha(u) = x(1 - [1-\alpha][1-u]) - x(u)\mbox{,}}
where \eqn{P_\alpha(u)} is the  \eqn{\alpha}-percentile residual quantile for nonexceedance probability \eqn{u} and percentile \eqn{\alpha} and \eqn{x(u)} is a constant for \eqn{x(F = u)}. The reversed \eqn{\alpha}-percentile residual quantile is available under \code{\link{rralmomco}}.
}
\usage{
ralmomco(f, para, alpha=0)
}
\arguments{
  \item{f}{Nonexceedance probability (\eqn{0 \le F \le 1}).}
  \item{para}{The parameters from \code{\link{lmom2par}} or \code{\link{vec2par}}.}
  \item{alpha}{The \eqn{\alpha} percentile, which is divided by \eqn{100} inside the function ahead of calling the quantile function of the distribution.}
}
\value{
   \eqn{\alpha}-percentile residual quantile value for \eqn{F}.
}
\author{W.H. Asquith}
\references{
Nair, N.U., and Vineshkumar, B., 2011, Reversed percentile residual life and related concepts: Journal of the Korean Statistical Society, v. 40, no. 1, pp. 85--92.

Nair, N.U., Sankaran, P.G., and Balakrishnan, N., 2013, Quantile-based reliability analysis: Springer, New York.
}
\seealso{\code{\link{qlmomco}}, \code{\link{rmlmomco}}, \code{\link{rralmomco}}
}
\examples{
# It is easiest to think about residual life as starting at the origin, units in days.
A <- vec2par(c(0.0, 2649, 2.11), type="gov") # so set lower bounds = 0.0
maximum.lifetime <- quagov(1,A) # 2649 days
ralmomco(0,A,alpha=0)   #    0 days
ralmomco(0,A,alpha=100) # 2649 days
ralmomco(1,A,alpha=0)   #    0 days (death certain)
ralmomco(1,A,alpha=100) #    0 days (death certain)
\dontrun{
F <- nonexceeds(f01=TRUE)
plot(F, qlmomco(F,A), type="l",
     xlab="NONEXCEEDANCE PROBABILITY", ylab="LIFETIME, IN DAYS")
lines(F, rmlmomco(F, A), col=4, lwd=4) # thick blue, residual mean life
lines(F, ralmomco(F, A, alpha=50), col=2) # solid red, median residual life
lines(F, ralmomco(F, A, alpha=10), col=2, lty=2) # lower dashed line,
                                              # the 10th percentile of residual life
lines(F, ralmomco(F, A, alpha=90), col=2, lty=2) # upper dashed line,
                                              # 10th percentile of residual life}
}
\keyword{quantile function}
\keyword{alpha-percentile residual quantile function}
\keyword{lifetime/reliability analysis}
\concept{lifetime}
\concept{lifetime}
\concept{reliability}
\keyword{The lmomco functions}
