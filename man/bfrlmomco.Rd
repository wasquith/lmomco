\name{bfrlmomco}
\alias{bfrlmomco}
\title{Bonferroni Curve of the Distributions}
\description{
This function computes the Bonferroni Curve for quantile function \eqn{x(F)} (\code{\link{par2qua}}, \code{\link{qlmomco}}). The function is defined by Nair et al. (2013, p. 179) as
\deqn{B(u) = \frac{1}{\mu u}\int_0^u x(p)\; \mathrm{d}p\mbox{,}}
where \eqn{B(u)} is Bonferroni curve for quantile function \eqn{x(F)} and \eqn{\mu} is the conditional mean for quantile \eqn{u=0} (\code{\link{cmlmomco}}). The Bonferroni curve is related to the Lorenz curve (\eqn{L(u)}, \code{\link{lrzlmomco}}) by
\deqn{B(u) = \frac{L(u)}{u}\mbox{.}}
}
\usage{
bfrlmomco(f, para)
}
\arguments{
  \item{f}{Nonexceedance probability (\eqn{0 \le F \le 1}).}
  \item{para}{The parameters from \code{\link{lmom2par}} or \code{\link{vec2par}}.}
}
\value{
  Bonferroni curve value for \eqn{F}.
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
f <- 0.65 # Both computations report: 0.5517342
Bu1 <- 1/(cmlmomco(f=0,A)*f) * integrate(afunc, 0, f)$value
Bu2 <- bfrlmomco(f, A)
}
\keyword{quantile function}
\keyword{Bonferroni curve}
\keyword{lifetime/reliability analysis}
\concept{lifetime}
\concept{lifetime}
\concept{reliability}
\keyword{The lmomco functions}
