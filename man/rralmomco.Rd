\name{rralmomco}
\alias{rralmomco}
\title{Reversed Alpha-Percentile Residual Quantile Function of the Distributions}
\description{
This function computes the Reversed \eqn{\alpha}-Percentile Residual Quantile Function for quantile function \eqn{x(F)} (\code{\link{par2qua}}, \code{\link{qlmomco}}). The function is defined by Nair and Vineshkumar (2011, p. 87) and Midhu et al. (2013, p. 13) as
\deqn{R_\alpha(u) = x(u) - x(u[1-\alpha])\mbox{,}}
where \eqn{R_\alpha(u)} is the reversed \eqn{\alpha}-percentile residual quantile for nonexceedance probability \eqn{u} and percentile \eqn{\alpha} and \eqn{x(u[1-\alpha])} is a constant for \eqn{x(F = u[1-\alpha])}. The nonreversed \eqn{\alpha}-percentile residual quantile is available under \code{\link{ralmomco}}.
}
\usage{
rralmomco(f, para, alpha=0)
}
\arguments{
  \item{f}{Nonexceedance probability (\eqn{0 \le F \le 1}).}
  \item{para}{The parameters from \code{\link{lmom2par}} or \code{\link{vec2par}}.}
  \item{alpha}{The \eqn{\alpha} percentile, which is divided by \eqn{100} inside the function ahead of calling the quantile function of the distribution.}
}
\value{
   Reversed \eqn{\alpha}-percentile residual quantile value for \eqn{F}.
}
\note{
Technically it seems that Nair et al. (2013) do not explictly define the reversed \eqn{\alpha}-percentile residual quantile but their index points to pp. 69--70 for a derivation involving the Generalized Lambda distribution (GLD) but that derivation (top of p. 70) has incorrect algebra.  A possibilty is that Nair et al. (2013) forgot to include \eqn{R_\alpha(u)} as an explicit definition in juxtaposition to \eqn{P_\alpha(u)} (\code{\link{ralmomco}}) and then apparently made an easy-to-see algebra error in trying to collect terms for the GLD.
}
\author{W.H. Asquith}
\references{
Nair, N.U., and Vineshkumar, B., 2011, Reversed percentile residual life and related concepts: Journal of the Korean Statistical Society, v. 40, no. 1, pp. 85--92.

Midhu, N.N., Sankaran, P.G., and Nair, N.U., 2013, A class of distributions with linear mean residual quantile function and it's generalizations: Statistical Methodology, v. 15, pp. 1--24.

Nair, N.U., Sankaran, P.G., and Balakrishnan, N., 2013, Quantile-based reliability analysis: Springer, New York.
}
\seealso{\code{\link{qlmomco}}, \code{\link{ralmomco}}
}
\examples{
# It is easiest to think about residual life as starting at the origin, units in days.
A <- vec2par(c(145, 2649, 2.11), type="gov") # so set lower bounds = 0.0
rralmomco(0.78, A, alpha=50)
\dontrun{
F <- nonexceeds(f01=TRUE); r <- range(rralmomco(F,A, alpha=50), ralmomco(F,A, alpha=50))
plot(F, rralmomco(F,A, alpha=50), type="l", xlab="NONEXCEEDANCE PROBABILITY",
                  ylim=r, ylab="MEDIAN RESIDUAL OR REVERSED LIFETIME, IN DAYS")
lines(F, ralmomco(F, A, alpha=50), col=2) # notice the lack of symmetry
}
}
\keyword{quantile function}
\keyword{alpha-percentile residual quantile function}
\keyword{lifetime/reliability analysis}
\concept{lifetime}
\concept{lifetime}
\concept{reliability}
\keyword{The lmomco functions}
