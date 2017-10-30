\encoding{utf8}
\name{dat2bernquaf}
\alias{dat2bernquaf}
\title{ Equivalent Nonexceedance Probability for a Given Value through Observed Data to Empirical Quantiles through Bernstein or Kantorovich Polynomials }
\description{
This function computes an equivalent nonexceedance probability \eqn{F} of a single value \eqn{x} for the sample data set (\eqn{\hat{X}}) through inversion of the empricial quantile function as computable through Bernstein or Kantorovich Polynomials by the \code{\link{dat2bernqua}} function.
}
\usage{
dat2bernquaf(x, data, interval=NA, ...)
}
\arguments{
  \item{x}{A scalar value for which the equivalent nonexceedance probability \eqn{F} through the function \code{\link{dat2bernqua}} is to be computed.}
  \item{data}{A vector of data values that directly correspond to the argument \code{x} in the function \code{\link{dat2bernqua}}.}
  \item{interval}{The search interval. If \code{NA}, then \eqn{[1/(n+1), 1 - 1/(n+1)]} is used. If \code{interval} is a single value \eqn{a}, then the interval is computed as \eqn{[a, 1 - a]}.}
  \item{...}{Additional arguments passed to \code{\link{dat2bernqua}} through the \code{uniroot()} function in \R.}
}
\details{
The basic logic is thus. The \eqn{\hat{X}} in conjunction with the settings for the polynomials provides the empirical quantile function (EQF). The \code{dat2bernquaf} function then takes the EQF (through dynamic recomputation) and seeks a root for the single value also given.

The critical piece likely is the search interval, which can be modified by the \code{interval} argument if the internal defaults are not sufficient. The default interval is determined as the first and last Weibull plotting positions of \eqn{\hat{X}} having a sample size \eqn{n}: \eqn{[1/(n+1), 1 - 1/(n+1)]}.  Because the \code{\link{dat2bernqua}} function has a substantial set of options that control how the empirical curve is (might be) extrapolated beyond the range of \eqn{\hat{X}}, it is difficult to determine an always suitable interval for the rooting.  However, it should be considered obvious that the result is more of an interpolation if \eqn{F(x)} is within \eqn{F \in  [1/(n+1), 1 - 1/(n+1)]} and increasingly becomes an accurate interpolation as \eqn{F(x) \rightarrow 1/2} (the median).

If the value \eqn{x} is too far beyond the data or if the search interval is not sufficient then the following error will be triggered:
\preformatted{
Error in uniroot(afunc, interval, ...) :
  f() values at end points not of opposite sign
}
The Examples section explores this aspect.
}
\value{
  An \R \code{list} is returned.
  \item{x}{An echoing of the \eqn{x} value via the \code{x} argument.}
  \item{f}{The equivalent nonexceedance probability \eqn{F(x{\mid}\hat{X})}.}
  \item{interval}{The search interval of \eqn{F} used.}
  \item{afunc.root}{Corresponds to the \code{f.root} element returned by the \code{uniroot()} function.}
  \item{iter}{Corresponds to the \code{iter} element returned by the \code{uniroot()} function.}
  \item{estim.prec}{Corresponds to the \code{estim.prec} element returned by the \code{uniroot()} function.}
  \item{source}{An attribute identifying the computational source: \dQuote{dat2bernquaf}}
}
\author{ W.H. Asquith}
\seealso{ \code{\link{dat2bernqua}} }
\examples{
dat2bernquaf(6, c(2,10)) # median 1/2 of 2 and 10 is 6 (trivial and fast)
\dontrun{
set.seed(5135)
lmr <- vec2lmom(c(1000, 400, 0.2, 0.3, 0.045))
par <- lmom2par(lmr, type="wak")
Q <- rlmomco(83, par) # n = 83 and extremely non-Normal data
lgQ <- max(Q) # 5551.052 by theory
dat2bernquaf(median(Q), Q)$f  # returns 0.5100523 (nearly 1/2)
dat2bernquaf(lgQ, Q)$f                   # unable to root
dat2bernquaf(lgQ, Q, bound.type="sd")$f  # unable to root
itf <- c(0.5, 0.99999)
f <- dat2bernquaf(lgQ, Q, interval=itf, bound.type="sd")$f
print(f) # F=0.9961118
qlmomco(f, par) # 5045.784 for the estimate F=0.9961118
# If we were not using the maximum and something more near the center of the
# distribution then that estimate would be closer to qlmomco(f, par).
# You might consider lqQ <- qlmomco(0.99, Q) # theoretical 99th percentile and
# let the random seed wander and see the various results.
}
}
\keyword{smoothing (Bernstein-Kantorovich)}
\keyword{smoothing (Bernstein polynomial)}
\keyword{smoothing (Kantorovich polynomial)}
\keyword{smoothing (Parzen weighting)}
\keyword{quantile function}

