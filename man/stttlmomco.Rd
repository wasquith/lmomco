\name{stttlmomco}
\alias{stttlmomco}
\title{Scaled Total Time on Test Transform of Distributions}
\description{
This function computes the Scaled Total Time on Test Transform Quantile Function for a quantile function \eqn{x(F)} (\code{\link{par2qua}}, \code{\link{qlmomco}}).  The TTT is defined by Nair et al. (2013, p. 173) as
\deqn{\phi(u) = \frac{1}{\mu}\left[(1-u)x(u) + \int_0^u x(p)\; \mathrm{d}p \right]\mbox{,}}
where \eqn{\phi(u)} is the scaled  total time on test for nonexceedance probability \eqn{u}, and \eqn{x(u)} is a constant  for \eqn{x(F = u)}. The \eqn{\phi(u)} is also expressible in terms of total time on test transform quantile function (\eqn{T(u)}, \code{\link{tttlmomco}}) as
\deqn{\phi(u) = \frac{T(u)}{\mu}\mbox{,}}
where \eqn{\mu} is the conditional mean (\code{\link{cmlmomco}}) at \eqn{u = 0} and the later definition is the basis for implementation in \pkg{lmomco}. The integral in the first definition is closely related to the structure of the reversed residual mean quantile function (\eqn{R(u)}, \code{\link{rrmlmomco}}).
}
\usage{
stttlmomco(f, para)
}
\arguments{
  \item{f}{Nonexceedance probability (\eqn{0 \le F \le 1}).}
  \item{para}{The parameters from \code{\link{lmom2par}} or \code{\link{vec2par}}.}
}
\value{
  Scaled total time on test value for \eqn{F}.
}
\author{W.H. Asquith}
\references{
Nair, N.U., Sankaran, P.G., and Balakrishnan, N., 2013, Quantile-based reliability analysis: Springer, New York.
}
\seealso{\code{\link{qlmomco}}, \code{\link{tttlmomco}}
}
\examples{
# It is easiest to think about residual life as starting at the origin,
# but for this example, let us set the lower limit at 100 days.
A <- vec2par(c(100, 2649, 2.11), type="gov")
f <- 0.47  # Both computations of Phi show 0.6455061
"afunc" <- function(p) { return(par2qua(p,A,paracheck=FALSE)) }
tmpa <- 1/cmlmomco(f=0, A); tmpb <- (1-f)*par2qua(f,A,paracheck=FALSE)
Phiu1 <- tmpa * ( tmpb + integrate(afunc,0,f)$value )
Phiu2 <- stttlmomco(f, A)
\dontrun{
# The TTT-plot (see Nair et al. (2013, p. 173))
n <- 30; X <- sort(rlmomco(n, A)); lmr <- lmoms(X)  # simulated lives and their L-moments
# recognize here that the "fit" is to the lifetime data themselves and not to special
# curves or projections of the data to other scales
"Phir" <- function(r, X, sort=TRUE) {
   n <- length(X); if(sort) X <- sort(X)
   if(r == 0) return(0) # can use 2:r as X_{0:n} is zero
   Tau.rOFn <- sapply(1:r, function(j) { Xlo <- ifelse((j-1) == 0, 0, X[(j-1)]);
                                         return((n-j+1)*(X[j] - Xlo)) })
   return(sum(Tau.rOFn))
}
Xbar <- mean(X); rOFn <- (1:n)/n # Nair et al. (2013) are clear r/n used in the Phi(u)
Phi <- sapply(1:n, function(r) { return(Phir(r,X, sort=FALSE)) }) / (n*Xbar)
layout(matrix(1:3, ncol=1))
plot(rOFn, Phi, type="b",
     xlab="NONEXCEEDANCE PROBABILITY", ylab="SCALED TOTAL TIME ON TEST")
lines(rOFn, stttlmomco(rOFn, A), lwd=2, col=8) # solid grey, the parent distribution
par1 <- pargov(lmr); par2 <- pargov(lmr, xi=min(X)) # notice attempt to "fit at minimum"
lines(pp(X), stttlmomco(rOFn, par1)) # now Weibull (i/(n+1)) being used for F via pp()
lines(pp(X), stttlmomco(rOFn, par2), lty=2) # perhaps better, but could miss short lives
F <- nonexceeds(f01=TRUE)
plot(pp(X), sort(X), xlab="NONEXCEEDANCE PROBABILITY", ylab="TOTAL TIME ON TEST (DAYS)")
lines(F, qlmomco(F, A), lwd=2, col=8) # the parent again
lines(F, qlmomco(F, par1), lty=1); lines(F, qlmomco(F, par2), lty=2) # two estimated fits
plot(F,  lrzlmomco(F, par2), col=2, type="l")  # Lorenz curve from L-moment fit (red)
lines(F, bfrlmomco(F, par2), col=3, lty=2) # Bonferroni curve from L-moment fit (green)
lines(F, lkhlmomco(F, par2), col=4, lty=4) # Leimkuhler curve from L-moment fit (blue)
lines(rOFn, Phi) # Scaled Total Time on Test
}
}
\keyword{quantile function}
\keyword{scaled total time on test transform}
\keyword{total time on test transform}
\keyword{lifetime/reliability analysis}
\concept{lifetime}
\concept{lifetime}
\concept{reliability}
\concept{ttt transform}
\keyword{The lmomco functions}
