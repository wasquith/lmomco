\name{rmlmomco}
\alias{rmlmomco}
\title{Mean Residual Quantile Function of the Distributions}
\description{
This function computes the Mean Residual Quantile Function for quantile function \eqn{x(F)} (\code{\link{par2qua}}, \code{\link{qlmomco}}). The function is defined by Nair et al. (2013, p. 51) as
\deqn{M(u) = \frac{1}{1-u}\int_u^1 [x(p) - x(u)]\; \mathrm{d}p\mbox{,}}
where \eqn{M(u)} is the mean residual quantile for nonexceedance probability \eqn{u} and \eqn{x(u)} is a constant for \eqn{x(F = u)}. The variance of \eqn{M(u)} is provided in \code{\link{rmvarlmomco}}.

The integration instead of from \eqn{0 \rightarrow 1} for the usual quantile function is \eqn{u \rightarrow 1}. Note that \eqn{x(u)} is a constant, so
\deqn{M(u) = \frac{1}{1-u}\int_u^1 x(p)\; \mathrm{d}p - x(u)\mbox{,}}
is equivalent and the basis for the implementation in \code{\link{rmlmomco}}. Assuming that \eqn{x(F)} is a life distribution, the \eqn{M(u)} is interpreted (see Nair et al. [2013, p. 51]) as the average remaining life beyond the \eqn{100(1-F)\%} of the distribution. Alternatively, \eqn{M(u)} is the mean residual life conditioned that survival to lifetime \eqn{x(F)} has occurred.

If \eqn{u = 0}, then \eqn{M(0)} is the expectation of the life distribution or in otherwords \eqn{M(0) = \lambda_1} of the parent quantile function.  If \eqn{u = 1}, then \eqn{M(u) = 0} (death has occurred)---there is zero residual life remaining. The implementation intercepts an intermediate \eqn{\infty} and returns 0 for \eqn{u = 1}.

The \eqn{M(u)} is referred to as a quantile function but this quantity is not to be interpreted as a type of probability distribution. The second example produces a \eqn{M(u)} that is not monotonic increasing with \eqn{u} and therefore it is immediately apparent that \eqn{M(u)} is not the quantile function of some probability distribution by itself. Nair et al. (2013) provide extensive details on quantile-based reliability analysis.
}
\usage{
rmlmomco(f, para)
}
\arguments{
  \item{f}{Nonexceedance probability (\eqn{0 \le F \le 1}).}
  \item{para}{The parameters from \code{\link{lmom2par}} or \code{\link{vec2par}}.}
}
\value{
  Mean residual value for \eqn{F}.
}
\note{
The Mean Residual Quantile Function is the first of many other functions and \dQuote{curves} associated with lifetime/reliability analysis operations that at their root use the quantile function (QF, \eqn{x(F)}) of a distribution. Nair et al. (2013) (NSB) is the authoritative text on which the following functions in \pkg{lmomco} were based
\tabular{lclr}{
\cr
Residual mean QF \tab \eqn{M(u)} \tab \code{\link{rmlmomco}} \tab NSB[p.51] \cr
Variance residual QF   \tab \eqn{V(u)} \tab \code{\link{rmvarlmomco}} \tab NSB[p.54] \cr
\eqn{\alpha}-percentile residual QF   \tab \eqn{P_\alpha(u)} \tab \code{\link{ralmomco}} \tab NSB[p.56] \cr
Reversed \eqn{\alpha}-percentile residual QF   \tab \eqn{R_\alpha(u)} \tab \code{\link{rralmomco}} \tab NSB[p.69--70] \cr
Reversed residual mean QF \tab  \eqn{R(u)} \tab \code{\link{rrmlmomco}} \tab NSB[p.57] \cr
Reversed variance residual  QF \tab  \eqn{D(u)} \tab \code{\link{rrmvarlmomco}} \tab NSB[p.58] \cr
Conditional mean QF \tab \eqn{\mu(u)} \tab \code{\link{cmlmomco}}  \tab NSB[p.68] \cr
Vitality function (see conditional mean)\cr
Total time on test transform QF \tab  \eqn{T(u)} \tab \code{\link{tttlmomco}} \tab NSB[p.171--172, 176] \cr
Scaled total time on test transform QF \tab  \eqn{\phi(u)} \tab \code{\link{stttlmomco}} \tab NSB[p.173] \cr
Lorenz curve \tab  \eqn{L(u)} \tab \code{\link{lrzlmomco}} \tab NSB[p.174] \cr
Bonferroni curve \tab  \eqn{B(u)} \tab \code{\link{bfrlmomco}} \tab NSB[p.179] \cr
Leimkuhler curve \tab  \eqn{K(u)} \tab \code{\link{lkhlmomco}} \tab NSB[p.181] \cr
Income gap ratio curve  \tab  \eqn{G(u)} \tab \code{\link{riglmomco}} \tab NSB[p.230] \cr
Mean life: \eqn{\mu \equiv \mu(0) \equiv \lambda_1(u=0)  \equiv \lambda_1}\cr
L-moments of residual life \tab \eqn{\lambda_r(u)} \tab \code{\link{reslife.lmoms}} \tab NSB[p.202]\cr
L-moments of reversed residual life \tab \eqn{{}_\mathrm{r}\lambda_r(u)} \tab \code{\link{rreslife.lmoms}} \tab NSB[p.211]\cr
}
}
\author{ W.H. Asquith}
\references{
Kupka, J., and Loo, S., 1989, The hazard and vitality measures of ageing: Journal of Applied Probability, v. 26, pp. 532--542.

Nair, N.U., Sankaran, P.G., and Balakrishnan, N., 2013, Quantile-based reliability analysis: Springer, New York.
}
\seealso{\code{\link{qlmomco}}, \code{\link{cmlmomco}}, \code{\link{rmvarlmomco}}
}
\examples{
# It is easiest to think about residual life as starting at the origin, units in days.
A <- vec2par(c(0.0, 2649, 2.11), type="gov") # so set lower bounds = 0.0
qlmomco(0.5, A)  # The median lifetime = 1261 days
rmlmomco(0.5, A) # The average remaining life given survival to the median = 861 days

# 2nd example with discussion points
F <- nonexceeds(f01=TRUE)
plot(F, qlmomco(F, A), type="l", # usual quantile plot as seen throughout lmomco
     xlab="NONEXCEEDANCE PROBABILITY", ylab="LIFETIME, IN DAYS")
lines(F, rmlmomco(F, A), col=2, lwd=3)           # mean residual life
L1 <- lmomgov(A)$lambdas[1]                      # mean lifetime at start/birth
lines(c(0,1), c(L1,L1), lty=2)                   # line "ML" (mean life)
# Notice how ML intersects M(F|F=0) and again later in "time" (about F = 1/4)  showing
# that this Govindarajulu as a peak mean resiudal life that is **greater** than the
# expected lifetime at start. The M(F) then tapers off to zero at infinity time (F=1).
# M(F) is non-monotonic for this example---not a proper probability dist.
}
\keyword{quantile function}
\keyword{mean residual quantile function}
\keyword{lifetime/reliability analysis}
\concept{lifetime}
\concept{lifetime}
\concept{reliability}
\keyword{The lmomco functions}
