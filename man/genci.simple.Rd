\name{genci.simple}
\alias{genci}
\alias{genci.simple}
\title{Generate (Estimate) Confidence Intervals for Quantiles of a Parent Distribution}
\description{
This function estimates the lower and upper limits of a specified confidence interval for a vector of nonexceedance probabilities \eqn{F} of a specified parent distribution [quantile function \eqn{Q(F,\theta)} with parameters \eqn{\theta}] using Monte Carlo simulation. The \eqn{F} are specified by the user. The user also provides \eqn{\Theta} of the parent distribution (see \code{\link{lmom2par}}). This function is a wrapper on \code{\link{qua2ci.simple}}; please consult the documentation for that function for further details of the simulations.
}
\usage{
genci.simple(para, n, f=NULL, level=0.90, edist="gno", nsim=1000,
             expand=FALSE, verbose=FALSE, showpar=FALSE, quiet=FALSE)
}
\arguments{
  \item{para}{The parameters from \code{\link{lmom2par}} or similar.}
  \item{n}{The sample size for each Monte Carlo simulation will use.}
  \item{f}{Vector of nonexceedance probabilities (\eqn{0 \le f \le 1}) of the quantiles for which the confidence interval are needed. If \code{NULL}, then the vector as returned by \code{\link{nonexceeds}} is used.}
  \item{level}{The confidence interval (\eqn{0 \le } \code{level} \eqn{ < 1}). The interval is specified as the size of the interval. The default is 0.90 or the 90th percentile. The function will return the 5th (\eqn{(1-0.90)/2}) and 95th (\eqn{1-(1-0.90)/2}) percentile cumulative probability of the error distribution for the parent quantile as specified by the nonexceedance probability argument (\code{f}). This argument is passed unused to \code{\link{qua2ci.simple}}.}
  \item{edist}{The model for the error distribution. Although the Normal (the default) commonly is  assumed in error analyses, it need not be, as support for other distributions supported by \pkg{lmomco} is available. The default is the Generalized Normal so the not only is the Normal possible but asymmetry is also accomodated (\code{\link{lmomgno}}).  For example, if the L-skew  (\eqn{\tau_4}) or L-kurtosis (\eqn{\tau_4}) values depart considerably from those of the Normal (\eqn{\tau_3 = 0} and \eqn{\tau_4 = 0.122602}), then the Generalized Normal or some alternative distribution would likely provide more reliable confidence interval estimation. This argument is passed unused to \code{\link{qua2ci.simple}}.}
  \item{nsim}{The number of simulations (replications) for the sample size \code{n} to perform. Much larger simulation numbers are recommended---see discussion about\cr \code{\link{qua2ci.simple}}. This argument is passed unused to \code{\link{qua2ci.simple}}. Users are encouraged to experiment with \code{\link{qua2ci.simple}} to get a feel for the value of \code{edist} and \code{nsim}.}
  \item{expand}{Should the returned values be expanded to include information relating to the distribution type and L-moments of the distribution at the corresponding nonexceedance probabilities---in other words the information necessary to reconstruct the reported confidence interval. The default is \code{FALSE}. If \code{expand=FALSE} then a single \code{data.frame} of the lower and upper limits along with the true quantile value of the parent is returned. If \code{expand=TRUE}, then a more complicated \code{list} containing multiple \code{data.frame}s is returned.}
  \item{verbose}{The verbosity of the operation of the function. This argument is passed unused to \code{\link{qua2ci.simple}}.}
  \item{showpar}{The parameters of the \code{edist} for each simulation for each \eqn{F} value passed to \code{\link{qua2ci.simple}} are printed. This argument is passed unused to \code{\link{qua2ci.simple}}.}
  \item{quiet}{Suppress incremental counter for a count down of the \eqn{F} values.}
}
\value{
  An \R \code{data.frame} or \code{list} is returned (see discussion of argument \code{expand}). The following elements could be available.

  \item{nonexceed}{A vector of \eqn{F} values, which is returned for convenience so that post operations such as plotting are easily coded.}
  \item{lwr}{The lower value of the confidence interval having nonexceedance probability equal to \code{(1-level)/2}.}
  \item{true}{The true quantile value from \eqn{Q(F,\theta)} for the corresponding \eqn{F} value.}
  \item{upr}{The upper value of the confidence interval having \eqn{F} equal to \code{1-(1-level)/2}.}
  \item{lscale}{The second L-moment (L-scale, \eqn{\lambda_2}) of the distribution of quantiles for the corresponding \eqn{F}. This value is included in the primary returned \code{data.frame} because it measures the fundamental sampling variability.}
  \item{parent}{The paraments of the parent distribution if \code{expand=TRUE}.}
  \item{edist}{The type of error distribution used to model the confidence interval if the argument \code{expand=TRUE} is set.}
  \item{elmoms}{The L-moment of the distribution of quantiles for the corresponding \eqn{F} if the argument  \code{expand=TRUE} is set.}
  \item{epara}{An environment containing the parameter lists of the error distribution fit to the \code{elmoms} for each of the \code{f} if the argument \code{expand=TRUE} is set.}
  \item{ifail}{A failure integer.}
  \item{ifailtext}{Text message associated with \code{ifail}.}
}
\author{ W.H. Asquith}
\seealso{\code{\link{genci}}, \code{\link{gen.freq.curves}}
}
\examples{
\dontrun{
# For all these examples, nsim is way too small.
mean   <- 0; sigma <- 100
parent <- vec2par(c(mean,sigma), type='nor') # make parameter object
f      <- c(0.5, 0.8, 0.9, 0.96, 0.98, 0.99) # nonexceed probabilities
# nsim is small for speed of example not accuracy.
CI     <- genci.simple(parent, n=10, f=f, nsim=20); FF <- CI$nonexceed
plot( FF, CI$true, type='l', lwd=2)
lines(FF, CI$lwr, col=2); lines(FF, CI$upr, col=3)

pdf("twoCIplots.pdf")
# The qnorm() call has been added to produce "normal probability"
# paper on the horizonal axis. The parent is heavy-tailed.
GEV  <- vec2par(c(10000,1500,-0.3), type='gev') # a GEV distribution
CI   <- genci.simple(GEV, n=20, nsim=200, edist='gno')
ymin <- log10(min(CI$lwr[! is.na(CI$lwr)]))
ymax <- log10(max(CI$upr[! is.na(CI$upr)]))
qFF  <- qnorm(CI$nonexceed) 
plot( qFF, log10(CI$true), type='l', ylim=c(ymin,ymax),lwd=2)
lines(qFF, log10(CI$lwr), col=2); lines(qFF, log10(CI$upr), col=3)
# another error distribution model
CI   <- genci.simple(GEV, n=20, nsim=200, edist='aep4')
lines(qFF,log10(CI$lwr),col=2,lty=2); lines(qFF,log10(CI$upr),col=3,lty=2)
dev.off() # }
}
\keyword{distribution (utility)}
