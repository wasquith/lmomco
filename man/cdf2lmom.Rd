\name{cdf2lmom}
\alias{cdf2lmom}
\title{ Compute an L-moment from Cumulative Distribution Function }
\description{
Compute a single L-moment from a cumulative distribution function. This function is sequentially called by \code{\link{cdf2lmoms}} to mimic the output structure for multiple L-moments seen by other L-moment computation functions in \pkg{lmomco}.

For \eqn{r = 1}, the quantile function is actually used for numerical integration to compute the mean. The expression for the mean is
\deqn{
\lambda_1 = \int_0^1 x(F)\; \mathrm{d} F\mbox{,}
}
for quantile function \eqn{x(F)} and nonexceedance probability \eqn{F}. For \eqn{r \ge 2}, the L-moments can be computed from the cumulative distribution function \eqn{F(x)} by
\deqn{
\lambda_r = \frac{1}{r}\sum_{j=0}^{r-2} (-1)^j {r-2 \choose j}{r \choose j+1} \int_{-\infty}^{\infty} \! [F(x)]^{r-j-1}\times [1 - F(x)]^{j+1}\; \mathrm{d}x\mbox{.}
}
This equation is described by Asquith (2011, eq. 6.8), Hosking (1996), and Jones (2004).
}
\usage{
cdf2lmom(r, para, fdepth=0, silent=TRUE, ...)
}
\arguments{
  \item{r}{The order of the L-moment.}
  \item{para}{The parameters from \code{\link{lmom2par}} or similar.}
  \item{fdepth}{The depth of the nonexceedance/exceedance probabilities to determine the lower and upper integration limits for the integration involving \eqn{F(x)} through a call to the \code{\link{par2qua}} function. The default of 0 implies the quantile for \eqn{F=0} and quantile for \eqn{F=1} as the respective lower and upper limits.}
  \item{silent}{A logical to be passed into \code{\link{cdf2lmom}} and then onto the \code{try} functions encompassing the \code{integrate} function calls.}
  \item{...}{Additional arguments to pass to \code{\link{par2qua}} and \code{\link{par2cdf}}.}
}
\value{
  The value for the requested L-moment is returned (\eqn{\lambda_r}).
}
\references{
Asquith, W.H., 2011, Distributional analysis with L-moment statistics using the R environment for statistical computing: Createspace Independent Publishing Platform, ISBN 978--146350841--8.

Hosking, J.R.M., 1996, Some theoretical results concerning L-moments: Research Report RC14492, IBM Research Division, T.J. Watson Research Center, Yorktown Heights, New York.

Jones, M.C., 2004, On some expressions for variance, covariance, skewness and L-moments: Journal of Statistical Planning and Inference, v. 126, pp. 97--106.
}
\author{ W.H. Asquith}
\seealso{\code{\link{cdf2lmoms}}}
\examples{
para <- vec2par(c(.9,.4), type="nor")
cdf2lmom(4, para) # summarize the value
}
\keyword{L-moment (theoretical)}
