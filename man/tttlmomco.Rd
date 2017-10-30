\name{tttlmomco}
\alias{tttlmomco}
\title{Total Time on Test Transform of Distributions}
\description{
This function computes the Total Time on Test Transform Quantile Function for a quantile function \eqn{x(F)} (\code{\link{par2qua}}, \code{\link{qlmomco}}).  The TTT is defined by Nair et al. (2013, p. 171--172, 176) has several expressions
\deqn{T(u) = \mu - (1 - u) M(u)\mbox{,}}
\deqn{T(u) = x(u) - u R(u)\mbox{,}}
\deqn{T(u) = (1-u) x(u) + \mu L(u)\mbox{,}}
where \eqn{T(u)} is the total time on test for nonexceedance probability \eqn{u},  \eqn{M(u)} is the residual mean quantile function (\code{\link{rmlmomco}}), \eqn{x(u)} is a constant for \eqn{x(F = u)}, \eqn{R(u)} is the reversed mean residual quantile function (\code{\link{rrmlmomco}}), \eqn{L(u)} is the Lorenz curve (\code{\link{lrzlmomco}}), and \eqn{\mu} as the following definitions
\deqn{\mu \equiv \lambda_1(u=0)\mbox{\ first L-moment of residual life for\ }u=0\mbox{,}}
\deqn{\mu \equiv \lambda_1(x(F))\mbox{\ first L-moment of the quantile function}\mbox{,}}
\deqn{\mu \equiv \mu(0)\mbox{\ conditional mean for\ }u=0\mbox{.}}
The definitions imply that within numerical tolerances that \eqn{\mu(0)} (\code{\link{cmlmomco}}) should be equal to \eqn{T(1)}, which means that the conditional mean that the 0th percentile in life has been reached equals that total time on test for the 100th percentile. The later can be interpreted as meaning that each of realization of the lifetime distribution for the respective sample size lived to its expected ordered lifetimes.
}
\usage{
tttlmomco(f, para)
}
\arguments{
  \item{f}{Nonexceedance probability (\eqn{0 \le F \le 1}).}
  \item{para}{The parameters from \code{\link{lmom2par}} or \code{\link{vec2par}}.}
}
\value{
  Total time on test value for \eqn{F}.
}
\note{
The second definition for \eqn{\mu} is used and in \pkg{lmomco} code the implementation for nonexceedance probability \code{f} and parameter object \code{para} is
\preformatted{
Tu <- par2qua(f, para) - f*rrmlmomco(f, para) # 2nd def.
}
but other possible implementations for the first and third definitions respectively are
\preformatted{
Tu <- cmlmomco(f=0, para) - (1-f)*rmlmomco(f, para) # 1st def.
Tu <- (1-f)*par2qua(f, para) + cmlmomco(f=0, para)*lrzlmomco(f, para) # 3rd def.
}
}
\author{ W.H. Asquith}
\references{
Nair, N.U., Sankaran, P.G., and Balakrishnan, N., 2013, Quantile-based reliability analysis: Springer, New York.
}
\seealso{\code{\link{qlmomco}}, \code{\link{rmlmomco}}, \code{\link{rrmlmomco}}, \code{\link{lrzlmomco}}
}
\examples{
# It is easiest to think about residual life as starting at the origin, units in days.
A <- vec2par(c(0.0, 2649, 2.11), type="gov") # so set lower bounds = 0.0
tttlmomco(0.5, A)  # The median lifetime = 859 days

f <- c(0.25,0.75) # All three computations report: 306.2951 and 1217.1360 days.
Tu1 <- cmlmomco(f=0, A) - (1-f)* rmlmomco(f, A)
Tu2 <-    par2qua(f, A) -    f * rrmlmomco(f, A)
Tu3 <- (1-f)*par2qua(f, A) + cmlmomco(f=0, A)*lrzlmomco(f, A)

if(abs(cmlmomco(0,A) - tttlmomco(1,A)) < 1E-4) {
   print("These two quantities should be nearly identical.\n")
}
}
\keyword{quantile function}
\keyword{total time on test transform}
\keyword{lifetime/reliability analysis}
\concept{lifetime}
\concept{lifetime}
\concept{reliability}
\concept{ttt transform}
\keyword{The lmomco functions}
