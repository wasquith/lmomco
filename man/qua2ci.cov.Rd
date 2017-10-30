\name{qua2ci.cov}
\alias{qua2ci.cov}
\title{Estimate a Confidence Interval for Quantiles of a Parent Distribution using Sample Variance-Covariances of L-moments}
\description{
This function estimates the lower and upper limits of a specified confidence interval for aribitrary quantile values for a sample \eqn{x} and a specified distribution form. The estimation is based on the sample variance-covariance structure of the L-moments (\code{\link{lmoms.cov}}) through a Monte Carlo approach. The quantile values, actually the nonexceedance probabilities (\eqn{F} for \eqn{0 \le F \le 1}), are specified by the user. The user provides type of parent distribution distribution and this form which will be fitted internal to the function.
}
\usage{
qua2ci.cov(x,f, type=NULL, nsim=1000,
                interval=c("confidence", "none"), level=0.90,
                asnorm=FALSE, altlmoms=NULL, flip=NULL, dimless=TRUE,
                usefastlcov=TRUE, nmom=5, getsimlmom=FALSE, verbose=FALSE, ...)
}
\arguments{
  \item{x}{A real value vector.}
  \item{f}{Nonexceedance probabilities (\eqn{0 \le F \le 1}) of the quantiles for which the confidence interval is needed.}
  \item{type}{Three character distribution type (for example, type='gev').}
  \item{nsim}{The number of simulations to perform. Large numbers produce more refined confidence limit estimates at the cost of CPU time. The default is anticipated to be large enough to semi-quantitatively interpret results without too much computational delay. Larger simulation numbers are recommended.}
  \item{interval}{The type of interval to compute. If \code{"none"}, then the simulated quantiles are returned at which point \emph{only} the first value in \eqn{f} or \code{f[1]} will be considered but a warning will be issued to remind the user. This option is nice for making boxplots of the quantile distribution.}
  \item{level}{The confidence interval (\eqn{0 \le } \code{level} \eqn{ < 1}). The interval is specified as the size of the interval for which the default is 0.90 or the 90th percentile. The function will return the 5th [\eqn{(1-0.90)/2}] and 95th [\eqn{(1-(1-0.90)/2)}] percentile cumulative probability of the simulated quantile distribution as specified by the nonexceedance probability argument.}
  \item{asnorm}{Use the mean and standard deviation of the simulated quantiles as parameters of the Normal distribution to estimate the confidence interval. Otherwise, a Bernstein polynomial approximation (\code{\link{dat2bernqua}}) to the empirical distribution of the simulated quantile distribution is used.}
  \item{altlmoms}{Alternative L-moments to rescale the simulated L-moments from the variance-covariance structure of the sample L-moments in \code{x}. These L-moments need to be an \pkg{lmomco} package L-moment object (e.g. \code{\link{lmoms}}). The presence of alternative L-moments will result in \code{dimless=TRUE}.}
  \item{flip}{A flipping or reflection value denoted as \eqn{\eta}. The values in \code{x} are flipped by this value (\eqn{y = \eta - x}) and analysis proceeds with flipped information, and then results are flipped back just prior to returning values with the exception that if \code{getsimlmom=TRUE} then the simultated L-moments are in \dQuote{flipped space.}}
  \item{dimless}{Perform the simulations in dimensionless space meaning that values in \code{x} are converted by \eqn{y = (x-\lambda_1)/\lambda_2} and simulation based on \eqn{y} and scale is returned on output according to the L-moments of \code{x} or the alternative L-moments in \code{altlmoms}. Scale is returned to the simulated L-moments, if returned by \code{getsimlmom=TRUE}, which is not fully parallel with the returned behavior when flipping is involved.}
  \item{usefastlcov}{A logical to use the function \code{Lmomcov()} from the \pkg{Lmoments} package to compute the sample variance-covariance matrices and not the much slower function \code{\link{lmoms.cov}} in the \pkg{lmomco} package.}
  \item{nmom}{The number of L-moments involved. This argument needs to be high enough to permit parameterization of the distribution in \code{type} but computational effort increases as \code{nmom} gets large.  This option is provided in conjunction with \code{getsimlmom=TRUE} to be able to get a \dQuote{wider set} of simulated L-moments returned than precisely required by the distribution. Also, some distributions might as part of their specific fitting algorithms, require inspection of higher L-moments than seemingly required than their numer of parameters suggests.}
  \item{getsimlmom}{A logical controlling whether the simulated L-moment matrix having \code{nsim} rows and \code{nmom} columns is returned instead of confidence limits.}
  \item{verbose}{The verbosity of the operation of the function.}
  \item{...}{Additional arguments to pass such as to \code{\link{lmom2par}}.}
}
\value{
  An \R \code{data.frame} is returned.

  \item{lwr}{The lower value of the confidence interval having nonexceedance probability equal to \eqn{(1-}\code{level}\eqn{)/2}.}
  \item{fit}{The fit of the quantile based on the L-moments of \code{x} and possibly by reflection controlled by \code{flip} or based on the alternative L-moments in \code{altlmoms} and again by the reflection controlled by \code{flip}.}
  \item{upr}{The upper  value of the confidence interval having nonexceedance probability equal to \eqn{1-(1-}\code{level}\eqn{)/2}.}
  \item{qua_med}{The median of the simulated quantiles.}
  \item{qua_mean}{The mean of the simulated quantiles for which the median and mean should be very close if the simulation size is large enough and the quantile distribution is symmetrical.}
  \item{qua_var}{The variance (\eqn{\sigma^2(F)}) of the simulated quantiles.}
  \item{qua_lam2}{The L-scale (\eqn{\lambda_2(F)}) of the simulated quantiles for which \eqn{\sigma^2(F) \approx \pi\times\lambda^2_2(F)}.}
}
\author{ W.H. Asquith}
\seealso{\code{\link{lmoms}}, \code{\link{lmoms.cov}}, \code{\link{qua2ci.simple}}
}
\examples{
\dontrun{
samsize <- 128; nsim <- 2000; f <- 0.999
wei <- parwei(vec2lmom(c(100,75,-.3)))
set.seed(1734); X <- rlmomco(samsize, wei); set.seed(1734)
tmp <- qua2ci.cov(X, f, type="wei", nsim=nsim)
print(tmp) # show results of one 2000 replicated Monte Carlo
# nonexceed     lwr    fit    upr  qua_med  qua_mean  qua_var  qua_lam2
#     0.999   310.4  333.2  360.2    333.6     334.3    227.3    8.4988
set.seed(1734)
qf <- qua2ci.cov(X, f, type="wei", nsim=nsim, interval="none") # another
boxplot(qf)
message(" quantile variance: ", round(tmp$qua_var,  digits=2),
        " compared to ", round(var(qf, na.rm=TRUE), digits=2))
set.seed(1734)
genci.simple(wei, n=samsize, f=f)
# nonexceed     lwr    fit    upr  qua_med  qua_mean  qua_var  qua_lam2
#     0.999   289.7  312.0  337.7    313.5     313.6    213.5    8.2330

#----------------------------------------
# Using X from above example, demonstrate that using dimensionless
# simulation that the results are the same.
set.seed(145); qua2ci.cov(X, 0.1, type="wei") # both outputs same
set.seed(145); qua2ci.cov(X, 0.1, type="wei", dimless=TRUE)
# nonexceed     lwr    fit    upr  qua_med  qua_mean  qua_var  qua_lam2
#       0.1  -78.62 -46.01 -11.39   -43.58    -44.38   416.04     11.54

#----------------------------------------
# Using X again, demonstration application of the flip and notice that just
# simple reversal is occurring and that the Weibull is a reversed GEV.
eta <- 0
set.seed(145); qua2ci.cov(X, 0.9, type="wei", nsim=nsim)
# nonexceed     lwr    fit    upr  qua_med  qua_mean  qua_var  qua_lam2
#       0.9   232.2  244.2  255.9    244.3     244.1    51.91    4.0635
set.seed(145); qua2ci.cov(X, 0.9, type="gev", nsim=nsim, flip=eta)
# nonexceed     lwr    fit    upr  qua_med  qua_mean  qua_var  qua_lam2
#       0.9   232.2  244.2  256.2    244.2     244.3    53.02    4.1088
# The values are slightly different, which likely represents a combination
# of numerics of the variance-covariance matrix because the Monte Carlo
# is seeded the same.

#----------------------------------------
# Using X again, removed dimension and have the function add it back.
lmr <- lmoms(X); Y <- (X - lmr$lambdas[1])/lmr$lambdas[2]
set.seed(145); qua2ci.cov(Y, 0.9, type="wei", altlmoms=lmr, nsim=nsim)
# nonexceed     lwr    fit    upr  qua_med  qua_mean  qua_var  qua_lam2
#       0.9   232.2  244.2  255.9    244.3     244.1    51.91   4.0635}
}
\keyword{quantile function (confidence limit)}
