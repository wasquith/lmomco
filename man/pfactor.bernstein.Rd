\encoding{utf8}
\name{pfactor.bernstein}
\alias{pfactor.bernstein}
\title{Estimation of Optimal p-factor of Distributional Support Estimation for Smoothed Quantiles from the Bernstein or Kantorovich Polynomials }
\description{
Compute the optimal p-factor through numerical integration of the smoothed empirical quantile function to estimate the L-moments of the distribution. This function attempts to report an optimal \dQuote{p-factor}  (author's term) for the given parent distribution in \code{para} based on estimating the crossing of the origin of an error between the given L-moment ratio \eqn{\tau_r} for 3, 4, and 5 that will come from either the distribution parameter object or given as an argument in \code{lmr.dist}.  The estimated support of the distribution is that shown by  Turnbull and Ghosh (2014) and is computed as follows
\deqn{\biggl(x_{0:n},\: x_{n+1:n}\biggr) = \biggl(x_{1:n} - \frac{(x_{2:n} - x_{1:n})}{(1 - p)^{-2} - 1},\: x_{n:n} + \frac{(x_{n:n} - x_{n-1:n})}{(1 - p)^{-2} - 1}\biggr)\mbox{,}}
where \eqn{p} is the p-factor. The support will honor natural bounds if given by either \code{fix.lower} or \code{fix.upper}. The polynomial type for smooth is provided in \code{poly.type}. These three arguments are the same as those for \code{\link{dat2bernqua}} and \code{\link{lmoms.bernstein}}. The statistic type used to measure central tendency of the errors for the \code{nsim} simulations per \eqn{p}. The function has its own hardwired p-factors to compute but these can be superseded by the \code{pfactors} argument. The \code{p.lo} and \code{p.hi} are the lower and upper bounds to truncate on immediately after the p-factors to use are assembled. These are made for three purposes: (1) protection against numerical problems for mathematical upper limits (unity), (2) to potentially provide for much faster execution if the user already knows the approximate optimal value for the p-factor, and (3) to potentially use this function in a direct optimization framework using the \R functions \code{optim} or \code{uniroot}. It is strongly suggested to keep \code{plot.em} set so the user can inspect the computations.
}
\usage{
pfactor.bernstein(para, x=NULL, n=NULL,
                        bern.control=NULL,
                        poly.type=c("Bernstein", "Kantorovich"),
                        stat.type=c("Mean", "Median"),
                        fix.lower=NULL, fix.upper=NULL,
                        lmr.dist=NULL, lmr.n=c("3", "4", "5"),
                        nsim=500, plot.em=TRUE, pfactors=NULL,
                        p.lo=.Machine$double.eps, p.hi=1)
}
\arguments{
  \item{para}{A mandatory \dQuote{parent} distribution defined by a usual \pkg{lmomco} distribution parameter object for a distribution. The simulations are based on this distribution, although optimization for \eqn{p} can be set to a different L-moment value by \code{lmr.dist}.}
  \item{x}{An optional vector of data values.}
  \item{n}{An optional sample size to run the simulations on. This value is computed by \code{length(x)} if \code{x} is provided. If set by argument, then that size supersedes the length of the optional observed sample.}
  \item{bern.control}{A \code{list} that holds \code{poly.type}, \code{stat.type}, \code{fix.lower}, and \code{fix.upper}. And this list will supersede the respective
values provided as separate arguments. There is an implicit \code{bound.type} of \code{"Carv"}.}
  \item{poly.type}{Same argument as for \code{\link{dat2bernqua}}.}
  \item{stat.type}{The central estimation statistic for each p-factor evaluated.}
  \item{fix.lower}{Same argument as for \code{\link{dat2bernqua}}.}
  \item{fix.upper}{Same argument as for \code{\link{dat2bernqua}}.}
  \item{lmr.dist}{This is the value for the \code{lmr.n} of the distribution in \code{para} unless explicitly set through \code{lmr.dist}.}
  \item{lmr.n}{The L-moment ratio number for p-factor optimization.}
  \item{nsim}{The number of simulations to run. Experiments suggest the default is adequate for reasonably small sample sizes---the simulation count can be reduced as \code{n} becomes large.}
  \item{plot.em}{A logical to trigger the diagnostic plot of the simulated errors and a smooth line through these errors.}
  \item{pfactors}{An optional vector of p-factors to loop through for the simulations. The vector computing internall is this is set to \code{NULL} seems to be more than adequate.}
  \item{p.lo}{An computational lower boundary for which the \code{pfactors} by argument or default are truncated to. The default for \code{lo} is to be quite small and does no truncate the default \code{pfactors}.}
  \item{p.hi}{An computational upper boundary for which the \code{pfactors} by argument or default are truncated to. The default for \code{hi} is unity, which is the true upper limit that results in a 0 slope between the \eqn{x_{0:n}} to \eqn{x_{1:n}} or \eqn{x_{n:n}} to \eqn{x_{n+1:n}} order statistics.}
}
\value{
  An \R \code{list} or \code{real} is returned.  If \code{pfactors} is a single value, then the single value for the error statistic is returned, otherwise the list described will be. If the returned \code{pfactor} is \code{NA}, then likely the smooth line did not cross zero and the reason the user should keep \code{plot.em=TRUE} and inspect the plot. Perhaps revisions to the arguments will become evident. The contents of the list are
  \item{pfactor}{The estimated value of \eqn{p} smoothed by \code{lowess} that has an error of zero, see \code{err.stat} as a function of \code{ps}.}
  \item{bounds.type}{\code{Carv}, which is the same bound type as needed by \code{\link{dat2bernqua}} and \cr \code{\link{lmoms.bernstein}}.}
  \item{poly.type}{The given \code{poly.type}.}
  \item{stat.type}{The given \code{stat.type}. The \dQuote{Mean} seems to be preferable.}
  \item{lmom.type}{A string of the L-moment type: \dQuote{Tau3}, \dQuote{Tau4}, \dQuote{Tau5}.}
  \item{fix.lower}{The given fixed lower boundary, which could stay \code{NULL}.}
  \item{fix.upper}{The given fixed upper boundary, which could stay \code{NULL}.}
  \item{source}{An attribute identifying the computational source of the L-moments: \dQuote{pfactor.bernstein}.}
  \item{ps}{The p-factors actually evaluated.}
  \item{err.stat}{The error statistic computed by \code{stat.type} of the simulated \eqn{\hat{\tau_r}} by integration provided by \code{\link{lmoms.bernstein}} minus the \dQuote{true} value \eqn{\tau_r} provided by either \code{para} or given by \code{lmr.dist} where \eqn{r} is \code{lmr.n}.}
  \item{err.smooth}{The \code{lowess}-smoothed values for \code{err.stat} and the \code{pfactor} comes from a linear interpolation of this smooth for the error being zero.}
}
\note{
Repeated application of this function for various \code{n} would result in the analyst having a vector of \eqn{n} and \eqn{p} (\code{pfactor}). The analyst could then fit a regression equation and refine the estimated \eqn{p(n)}. For example, a dual-logarithmic regression is suggested \code{lm(log(p)~log(n))}.

Also, symmetrical data likely see little benefit from optimizing on the symmetry-measuring L-moments Tau3 and Tau5; the analyst might prefer to optimize on peakedness measured by Tau4.
}
\note{
This function is highly experimental and subject to extreme overhaul. Please contact the author if you are an interested party in Bernstein and Kantorovich polynomials.
}
\author{W.H. Asquith}
\references{
Turnbull, B.C., and Ghosh, S.K., 2014, Unimodal density estimation using Bernstein polynomials. Computational Statistics and Data Analysis, v. 72, pp. 13--29.
}
\seealso{\code{\link{lmoms.bernstein}}, \code{\link{dat2bernqua}}, \code{\link{lmoms}} }
\examples{
\dontrun{
pdf("pfactor_exampleB.pdf")
X <- exp(rnorm(200)); para <- parexp(lmoms(X))
# nsim is too small, but makes the following three not take too long
pfactor.bernstein(para, n=20, lmr.n="3", nsim=100, p.lo=.06, p.hi=.3)
pfactor.bernstein(para, n=20, lmr.n="4", nsim=100, p.lo=.06, p.hi=.3)
pfactor.bernstein(para, n=20, lmr.n="5", nsim=100, p.lo=.06, p.hi=.3)
dev.off()
}
\dontrun{
# Try intra-sample p-factor optimization from two perspectives. The 3-parameter
# GEV "over fits" the data and provides the parent.  Then use Tau3 of the fitted
# GEV for peakedness restraint and then use Tau3 of the data. Then repeat but use
# the apparent "exact" value of Tau3 for the true exponential parent.
pdf("pfactor_exampleB.pdf")
lmr <- vec2lmom(c(60,20)); paraA <- parexp(lmr); n <- 40
tr <- lmorph(par2lmom(paraA))$ratios[3]
X <- rlmomco(n, paraA); para <- pargev(lmoms(X))
F <- seq(0.001,0.999, by=0.001)
plot(qnorm(pp(X, a=0.40)), sort(X), type="n", log="y",
      xlab="Standard normal variate", ylab="Quantile",
      xlim=qnorm(range(F)), ylim=range(qlmomco(F,paraA)))
lines(qnorm(F), qlmomco(F, paraA), col=8, lwd=2)
lines(qnorm(F), qlmomco(F, para), lty=2)
points(qnorm(pp(X, a=0.40)), sort(X))

# Make sure to fill in the p-factor when needed!
bc <- list(poly.type = "Bernstein", bound.type="Carv",
           stat.type="Mean", fix.lower=0, fix.upper=NULL, p=NULL)
kc <- list(poly.type = "Kantorovich", bound.type="Carv",
           stat.type="Mean", fix.lower=0, fix.upper=NULL, p=NULL)

# Bernstein
A <- pfactor.bernstein(para,      n=n, nsim=100,              bern.control=bc)
B <- pfactor.bernstein(para, x=X, n=n, nsim=100,              bern.control=bc)
C <- pfactor.bernstein(para,      n=n, nsim=100, lmr.dist=tr, bern.control=bc)
D <- pfactor.bernstein(para, x=X, n=n, nsim=100, lmr.dist=tr, bern.control=bc)
plot(qnorm(pp(X, a=0.40)), sort(X), type="n", log="y",
      xlab="Standard normal variate", ylab="Quantile",
      xlim=qnorm(range(F)), ylim=range(qlmomco(F,paraA)))
lines(qnorm(F), qlmomco(F, paraA), col=8, lwd=2)
lines(qnorm(F), qlmomco(F, para), lty=2)
points(qnorm(pp(X, a=0.40)), sort(X))
      bc$p <- A$pfactor
lines(qnorm(F), dat2bernqua(F,X, bern.control=bc), col=2)
      bc$p <- B$pfactor
lines(qnorm(F), dat2bernqua(F,X, bern.control=bc), col=3)
      bc$p <- C$pfactor
lines(qnorm(F), dat2bernqua(F,X, bern.control=bc), col=2, lty=2)
      bc$p <- D$pfactor
lines(qnorm(F), dat2bernqua(F,X, bern.control=bc), col=3, lty=2)
# Kantorovich
A <- pfactor.bernstein(para,      n=n, nsim=100,              bern.control=kc)
B <- pfactor.bernstein(para, x=X, n=n, nsim=100,              bern.control=kc)
C <- pfactor.bernstein(para,      n=n, nsim=100, lmr.dist=tr, bern.control=kc)
D <- pfactor.bernstein(para, x=X, n=n, nsim=100, lmr.dist=tr, bern.control=kc)
plot(qnorm(pp(X, a=0.40)), sort(X), type="n", log="y",
      xlab="Standard normal variate", ylab="Quantile",
      xlim=qnorm(range(F)), ylim=range(qlmomco(F,paraA)))
lines(qnorm(F), qlmomco(F, paraA), col=8, lwd=2)
lines(qnorm(F), qlmomco(F, para), lty=2)
points(qnorm(pp(X, a=0.40)), sort(X))
      kc$p <- A$pfactor
lines(qnorm(F), dat2bernqua(F,X, bern.control=kc), col=2)
      kc$p <- B$pfactor
lines(qnorm(F), dat2bernqua(F,X, bern.control=kc), col=3)
      kc$p <- C$pfactor
lines(qnorm(F), dat2bernqua(F,X, bern.control=kc), col=2, lty=2)
      kc$p <- D$pfactor
lines(qnorm(F), dat2bernqua(F,X, bern.control=kc), col=3, lty=2)
dev.off()
}
\dontrun{
X <- exp(rnorm(200)); para <- parexp(lmoms(X))
"pfactor.root" <- function(para, p.lo, p.hi, ...) {
    afunc <- function(p, para=NULL, x=NULL, ...) {
      return(pfactor.bernstein(para=para, x=x, pfactors=p, ...)) }
    rt <- uniroot(afunc, c(p.lo, p.hi),
                  tol=0.001, maxiter=30, nsim=500, para=para, ...)
    return(rt)
}
pfactor.root(para, 0.05, 0.15, n=10, lmr.n="4")
pfactor.bernstein(para, n=10, lmr.n="4", nsim=200, p.lo=.05, p.hi=.15)
}
}
\keyword{smoothing (Bernstein-Kantorovich)}
