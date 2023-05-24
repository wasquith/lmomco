\name{pp}
\alias{pp}
\title{Plotting-Position Formula}
\description{
The plotting positions of a data vector (\code{x}) are returned in ascending order. The plotting-position formula is
\deqn{pp_i = \frac{i-a}{n+1-2a} \mbox{,}}
where \eqn{pp_i} is the nonexceedance probability \eqn{F} of the \eqn{i}th ascending data value. The parameter \eqn{a} specifies the plotting-position type, and \eqn{n} is the sample size (\code{length(x)}). Alternatively, the plotting positions can be computed by
\deqn{pp_i = \frac{i+A}{n+B} \mbox{,}}
where \eqn{A} and \eqn{B} can obviously be expressed in terms of \eqn{a}.
}
\usage{
pp(x, A=NULL, B=NULL, a=0, sort=TRUE, ties.method="first", ...)
}
\arguments{
  \item{x}{A vector of data values. The vector is used to get sample size through \code{length}.}
  \item{A}{A value for the plotting-position coefficient \eqn{A}.}
  \item{B}{A value for the plotting-position coefficient \eqn{B}.}
  \item{a}{A value for the plotting-position formula from which \eqn{A} and \eqn{B} are computed, default is \code{a=0}, which returns the Weibull plotting positions.}
  \item{sort}{A logical whether the ranks of the data are sorted prior to \eqn{F} computation. It was a design mistake years ago to default this function to a sort, but it is now far too late to risk changing the logic now. The function originally lacked the \code{sort} argument for many years.}
  \item{ties.method}{This is the argument of the same name passed to \code{rank}.}
  \item{...}{Additional arguments to pass.}
}
\value{
  An \R \code{vector} is returned.
}
\note{
Various plotting positions have been suggested in the literature. Stedinger and others (1992, p.18.25) comment that \dQuote{all plotting positions give crude estimates of the unknown [non]exceedance probabilities associated with the largest (and smallest) events.} The various plotting positions are summarized in the follow table.

\describe{
 \item{Weibull}{\eqn{a=0}, Unbiased exceedance probability for all distributions (see discussion in \code{\link{pp.f}}).}
 \item{Median}{\eqn{a=0.3175}, Median exceedance probabilities for all distributions (if so, see \code{\link{pp.median}}).}
 \item{APL}{\eqn{\approx 0.35}, Often used with probability-weighted moments.}
 \item{Blom}{\eqn{a=0.375}, Nearly unbiased quantiles for normal distribution.}
 \item{Cunnane}{\eqn{a=0.40}, Approximately quantile unbiased.}
 \item{Gringorten}{\eqn{a=0.44}, Optimized for Gumbel distribution.}
 \item{Hazen}{\eqn{a=0.50}, A traditional choice.}
}

The function uses the \R \code{rank} function, which has specific settings to handle tied data. For implementation here, the \code{ties.method="first"} method to \code{rank} is used. The user has flexibility in changing this to their own custom purposes.
}
\references{
Stedinger, J.R., Vogel, R.M., and Foufoula-Georgiou, E., 1992, Frequency analysis of extreme events, in Handbook of Hydrology, chapter 18, editor-in-chief D. A. Maidment: McGraw-Hill, New York.
}
\author{W.H. Asquith}
\seealso{\code{\link{nonexceeds}}, \code{\link{pwm.pp}}, \code{\link{pp.f}}, \code{\link{pp.median}}, \code{\link{headrick.sheng.lalpha}} }
\examples{
Q  <- rnorm(20)
PP <- pp(Q)
plot(PP, sort(Q))

Q <- rweibull(30, 1.4, scale=400)
WEI <- parwei(lmoms(Q))
PP <- pp(Q)
plot( PP, sort(Q))
lines(PP, quawei(PP, WEI))

# This plot looks similar, but when connecting lines are added
# the nature of the sorting is obvious.
plot( pp(Q, sort=FALSE), Q)
lines(pp(Q, sort=FALSE), Q, col=2)
}
\keyword{plotting position}
