\name{cmlmomco}
\alias{cmlmomco}
\title{Conditional Mean Residual Quantile Function of the Distributions}
\description{
This function computes the Conditional Mean Residual Quantile Function for quantile function \eqn{x(F)} (\code{\link{par2qua}}, \code{\link{qlmomco}}). The function is defined by Nair et al. (2013, p. 68) as
\deqn{\mu(u) = \frac{1}{1-u}\int_u^1 x(p)\; \mathrm{d}p\mbox{,}}
where \eqn{\mu(u)} is the conditional mean for nonexceedance probability \eqn{u}. The \eqn{\mu(u)} is the expectation \eqn{\mathrm{E}[X | X > x]}.
The \eqn{\mu(u)} also is known as the \emph{vitality function}. Details can be found in Nair et al. (2013, p. 68) and Kupka and Loo (1989). Mathematically, the vitality function simply is
\deqn{\mu(u) = M(u) + x(u)\mbox{,}}
where \eqn{M(u)} is the mean residual quantile function (\code{\link{rmlmomco}}), \eqn{x(u)} is a constant for \eqn{x(F = u)}.
}
\usage{
cmlmomco(f, para)
}
\arguments{
  \item{f}{Nonexceedance probability (\eqn{0 \le F \le 1}).}
  \item{para}{The parameters from \code{\link{lmom2par}} or \code{\link{vec2par}}.}
}
\value{
  Conditional mean residual value for \eqn{F} or conditional mean life for \eqn{F}.
}
\author{W.H. Asquith}
\references{
Kupka, J., and Loo, S., 1989, The hazard and vitality measures of ageing: Journal of Applied Probability, v. 26, pp. 532--542.

Nair, N.U., Sankaran, P.G., and Balakrishnan, N., 2013, Quantile-based reliability analysis: Springer, New York.
}
\seealso{\code{\link{qlmomco}}, \code{\link{rmlmomco}}
}
\examples{
# It is easiest to think about residual life as starting at the origin, units in days.
A <- vec2par(c(0.0, 2649, 2.11), type="gov") # so set lower bounds = 0.0
qlmomco(0.5, A)  # The median lifetime = 1261 days
rmlmomco(0.5, A) # The average remaining life given survival to the median = 861 days
cmlmomco(0.5, A) # The average total life given survival to the median = 2122 days

# Now create with a nonzero origin
A <- vec2par(c(100, 2649, 2.11), type="gov") # so set lower bounds = 0.0
qlmomco(0.5, A)  # The median lifetime = 1361 days
rmlmomco(0.5, A) # The average remaining life given survival to the median = 861 days
cmlmomco(0.5, A) # The average total life given survival to the median = 2222 days

# Mean life (mu), which shows up in several expressions listed under rmlmomco.
mu1 <- cmlmomco(0,A)
mu2 <- par2lmom(A)$lambdas[1]
mu3 <- reslife.lmoms(0,A)$lambdas[1]
# Each mu is 1289.051 days.
}
\keyword{quantile function}
\keyword{conditional mean quantile function}
\keyword{lifetime/reliability analysis}
\concept{lifetime}
\concept{lifetime}
\concept{reliability}
\keyword{The lmomco functions}
