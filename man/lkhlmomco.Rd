\name{lkhlmomco}
\alias{lkhlmomco}
\title{Leimkuhler Curve of the Distributions}
\description{
This function computes the Leimkuhler Curve  for quantile function \eqn{x(F)} (\code{\link{par2qua}}, \code{\link{qlmomco}}). The function is defined by Nair et al. (2013, p. 181) as
\deqn{K(u) = 1 - \frac{1}{\mu}\int_0^{1-u} x(p)\; \mathrm{d}p\mbox{,}}
where \eqn{K(u)} is Leimkuhler curve for nonexceedance probability \eqn{u}. The Leimkuhler curve is related to the Lorenz curve (\eqn{L(u)}, \code{\link{lrzlmomco}}) by
\deqn{K(u) = 1-L(1-u)\mbox{,}}
and related to the reversed residual mean quantile function (\eqn{R(u)}, \code{\link{rrmlmomco}}) and conditional mean (\eqn{\mu}, \code{\link{cmlmomco}}) for \eqn{u=0} by
\deqn{K(u) = \frac{1}{\mu} [\mu - (1-u)(x(1-u) - R(1-u))] \mbox{.}}
}
\usage{
lkhlmomco(f, para)
}
\arguments{
  \item{f}{Nonexceedance probability (\eqn{0 \le F \le 1}).}
  \item{para}{The parameters from \code{\link{lmom2par}} or \code{\link{vec2par}}.}
}
\value{
  Leimkuhler curve value for \eqn{F}.
}
\author{W.H. Asquith}
\references{
Nair, N.U., Sankaran, P.G., and Balakrishnan, N., 2013, Quantile-based reliability analysis: Springer, New York.
}
\seealso{\code{\link{qlmomco}}, \code{\link{lrzlmomco}}
}
\examples{
# It is easiest to think about residual life as starting at the origin, units in days.
A <- vec2par(c(0.0, 2649, 2.11), type="gov") # so set lower bounds = 0.0

"afunc" <- function(u) { return(par2qua(u,A,paracheck=FALSE)) }
f <- 0.35 # All three computations report: Ku = 0.6413727
Ku1 <- 1 - 1/cmlmomco(f=0,A) * integrate(afunc,0,1-f)$value
Ku2 <- (cmlmomco(0,A) - (1-f)*(quagov(1-f,A) - rrmlmomco(1-f,A)))/cmlmomco(0,A)
Ku3 <- lkhlmomco(f, A)
}
\keyword{quantile function}
\keyword{Leimkuhler curve}
\keyword{lifetime/reliability analysis}
\concept{lifetime}
\concept{lifetime}
\concept{reliability}
\keyword{The lmomco functions}
