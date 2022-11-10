\name{lrzlmomco}
\alias{lrzlmomco}
\title{Lorenz Curve  of the Distributions}
\description{
This function computes the Lorenz Curve for quantile function \eqn{x(F)} (\code{\link{par2qua}}, \code{\link{qlmomco}}). The function is defined by Nair et al. (2013, p. 174) as
\deqn{L(u) = \frac{1}{\mu}\int_0^u x(p)\; \mathrm{d}p\mbox{,}}
where \eqn{L(u)} is the Lorenz curve for nonexceedance probability \eqn{u}. The Lorenz curve is related to the Bonferroni curve (\eqn{B(u)}, \code{\link{bfrlmomco}}) by
\deqn{L(u) = \mu B(u)\mbox{.}}
}
\usage{
lrzlmomco(f, para)
}
\arguments{
  \item{f}{Nonexceedance probability (\eqn{0 \le F \le 1}).}
  \item{para}{The parameters from \code{\link{lmom2par}} or \code{\link{vec2par}}.}
}
\value{
  Lorzen curve value for \eqn{F}.
}
\author{W.H. Asquith}
\references{
Nair, N.U., Sankaran, P.G., and Balakrishnan, N., 2013, Quantile-based reliability analysis: Springer, New York.
}
\seealso{\code{\link{qlmomco}}, \code{\link{bfrlmomco}}
}
\examples{
# It is easiest to think about residual life as starting at the origin, units in days.
A <- vec2par(c(0.0, 2649, 2.11), type="gov") # so set lower bounds = 0.0
f <- c(0.25, 0.75) # Both computations report: 0.02402977 and 0.51653731
Lu1 <-   lrzlmomco(f, A)
Lu2 <- f*bfrlmomco(f, A)

# The Lorenz curve is related to the Gini index (G), which is L-CV:
"afunc" <- function(u) { return(lrzlmomco(f=u, A)) }
L <- integrate(afunc, lower=0, upper=1)$value
G <- 1 - 2*L                                                    # 0.4129159
G <- 1 - expect.min.ostat(2,para=A,qua=quagov)*cmlmomco(f=0,A)  # 0.4129159
LCV <- lmomgov(A)$ratios[2]                                     # 0.41291585
}
\keyword{quantile function}
\keyword{Lorenz curve}
\keyword{lifetime/reliability analysis}
\concept{lifetime}
\concept{lifetime}
\concept{reliability}
\keyword{The lmomco functions}
