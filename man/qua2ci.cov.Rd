\name{qua2ci.cov}
\alias{qua2ci.cov}
\title{Estimate a Confidence Interval for Quantiles of a Parent Distribution using Sample Variance-Covariances of L-moments}
\description{
This function estimates the lower and upper limits of a specified confidence interval for aribitrary quantile values for a sample \eqn{x} and a specified distribution form. The estimation is based on the sample variance-covariance structure of the L-moments (\code{\link{lmoms.cov}}) through a Monte Carlo approach. The quantile values, actually the nonexceedance probabilities (\eqn{F} for \eqn{0 \le F \le 1}), are specified by the user. The user provides type of parent distribution distribution and this form which will be fitted internal to the function.
}
\usage{
qua2ci.cov(x,f, type=NULL, nsim=1000,
                interval=c("confidence", "none"), level=0.90, tol=1E-6,
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
  \item{tol}{The tolerance argument of same name and default to feed to \code{MASS::mvrnorm()} and try increasing this tolerance if the error \dQuote{'Sigma' is not positive definite} occurs (see \bold{Note} for more discussion).}
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
\note{
These particular data set needs further evaluation as these particular sample can produce non-positive definite matrix being fed to \code{MASS:mvrnorm()}. It is noted that there are no ties in this data set.
\preformatted{
  test_dat <- c(0.048151736, 0.036753258, 0.034895847, 0.082792447, 0.096984927,
                0.213977450, 0.020264292, 0.269585438, 0.304746113, 0.066339093,
                0.015651114, 0.025122412, 0.184095698, 0.047167958, 0.049824752,
                0.043390768, 0.055228680, 0.009325696, 0.042145010, 0.008113992,
                0.118901521, 0.050399301, 0.049646181, 0.032299402, 0.015229284,
                0.013684668, 0.049371734, 0.068426211, 0.207159600, 0.087228473,
                0.306276783, 0.024870356, 0.016946801, 0.051553444, 0.017654117)
  qua2ci.cov(test_dat, 0.5, type="pe3", tol=1E-6, nmom=5) # fails

  lams <- lmoms(    test_dat)$lambdas
  lamc <- lmoms.cov(test_dat)
  n <- 100
  set.seed(1)
  MV1 <- mvtnorm::rmvnorm(n, mean=lams, sigma=lamc, method="eigen")
  MV1 <- mvtnorm::rmvnorm(n, mean=lams, sigma=lamc, method="chol")
  MV1 <- mvtnorm::rmvnorm(n, mean=lams, sigma=lamc, method="svd")
  colnames(MV1) <- paste0(rep("lam",5),1:5)
  set.seed(1)
  MV2 <- MASS::mvrnorm(n, lams, lamc, tol=5E-2)
  set.seed(1)
  MV3 <- MASS::mvrnorm(n, lams, lamc, tol=Inf)

  summary(MV2-MV3)
  summary(MV1)
  summary(MV2)
  plotlmrdia(lmrdia(), xlim=c(0.3,0.7), ylim=c(0,.6))
  points(MV1[,3]/MV1[,2], MV1[,4]/MV1[,2], col="red",  cex=0.5)
  points(MV2[,3]/MV2[,2], MV2[,4]/MV2[,2], col="blue", cex=0.5)
}

Next we, try focusing on the upper left corner of the matrix, after all we do not need beyond the 3rd moment because the Pearson III is being used.
\preformatted{
  qua2ci.cov(test_dat, 0.5, type="pe3", tol=1E-6, nmom=3) # fails
}

Now try increasing the tolerance setting on the matrix postive definite test in the \code{MASS::mvrnorm()} function.
\preformatted{
  qua2ci.cov(test_dat, 0.5, type="pe3", tol=1E-4, nmom=5) # fails
}

Now try again just focusing on the upper left corner that we really need.
\preformatted{
  set.seed(1)
  qua2ci.cov(test_dat, 0.5, type="pe3", tol=1E-4, nmom=3) # IT WORKS
  # nonexceed     lwr      fit      upr  qua_med qua_mean   qua_var qua_lam2
  #       0.5 0.02762 0.044426 0.061189 0.044322 0.044319 0.0001019 0.005672
}

Let us now try a hack of smoothing the data through the Bernstein polynomial. Perhaps subtle issues in the data can be \dQuote{fixed} by this and the seed has been set to have the \code{MASS::mvrnorm()} see the same seed although the variance-covariance matrix is slightly changing. Notice that the tolerance now returns to the default and that we are requesting up through the 5th L-moment.
\preformatted{
  set.seed(1)
  n <- length(test_dat)
  smth_dat <- dat2bernqua((1:n)/(n+1), test_dat)
  qua2ci.cov(smth_dat, 0.5, type="pe3", tol=1E-6, nmom=5) # IT WORKS
  # nonexceed     lwr      fit     upr  qua_med qua_mean   qua_var  qua_lam2
  #       0.5 0.02864 0.048288 0.06778 0.048406 0.048201 0.0001405 0.0066678
}

A quick look at the smoothing. The author is not advocating for this but this trick might be useful in data-mining scale work where for some samples, we need something back. The user might then consider using the differences \code{upr}\eqn{-}\code{fit} and \code{fit}\eqn{-}\code{lwr} to reconstruct the interval from a fit based on the original sample.
\preformatted{
  plot( (1:n)/(n+1), sort(test_dat))
  lines((1:n)/(n+1), smth_dat, col=2)
}
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
