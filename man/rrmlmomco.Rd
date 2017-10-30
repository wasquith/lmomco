\name{rrmlmomco}
\alias{rrmlmomco}
\title{Reversed Mean Residual Quantile Function of the Distributions}
\description{
This function computes the Reversed Mean Residual Quantile Function for quantile function \eqn{x(F)} (\code{\link{par2qua}}, \code{\link{qlmomco}}). The function is defined by Nair et al. (2013, p.57) as
\deqn{R(u) = x(u) - \frac{1}{u}\int_0^u x(p)\; \mathrm{d}p\mbox{,}}
where \eqn{R(u)} is the reversed mean residual for nonexceedance probability \eqn{u} and \eqn{x(u)} is a constant for \eqn{x(F = u)}.
}
\usage{
rrmlmomco(f, para)
}
\arguments{
  \item{f}{Nonexceedance probability (\eqn{0 \le F \le 1}).}
  \item{para}{The parameters from \code{\link{lmom2par}} or \code{\link{vec2par}}.}
}
\value{
  Reversed mean residual value for \eqn{F}.
}
\author{ W.H. Asquith}
\references{
Nair, N.U., Sankaran, P.G., and Balakrishnan, N., 2013, Quantile-based reliability analysis: Springer, New York.
}
\seealso{\code{\link{qlmomco}}, \code{\link{rrmvarlmomco}}
}
\examples{
# It is easiest to think about residual life as starting at the origin, units in days.
A <- vec2par(c(0.0, 2649, 2.6), type="gov") # so set lower bounds = 0.0
qlmomco(0.5, A)  # The median lifetime = 1005 days
rrmlmomco(0.5, A) # The reversed mean remaining life given median survival = 691 days

\dontrun{
F <- nonexceeds(f01=TRUE)
plot(F, qlmomco(F,A), type="l", # life
     xlab="NONEXCEEDANCE PROBABILITY", ylab="LIFETIME, IN DAYS")
lines(F,  rmlmomco(F, A), col=4, lwd=4) # thick blue, mean residual life
lines(F, rrmlmomco(F, A), col=2, lwd=4) # thick red, reversed mean residual life
}
}
\keyword{quantile function}
\keyword{reversed mean residual quantile function}
\keyword{lifetime/reliability analysis}
\concept{lifetime}
\concept{lifetime}
\concept{reliability}
\keyword{The lmomco functions}
