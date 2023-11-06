\name{theopwms}
\alias{theopwms}
\title{The Theoretical Probability-Weighted Moments using Integration of the Quantile Function }
\description{
Compute the theoretical probability-weighted moments (PWMs) for a distribution. A theoretical PWM in integral form is
\deqn{ \beta_r = \int^1_0 x(F)\,F^r\,\mathrm{d}F \mbox{,}}
where \eqn{x(F)} is the quantile function of the random variable \eqn{X} for nonexceedance probability \eqn{F} and \eqn{r} represents the order of the PWM. This function loops across the above equation for each \code{nmom} set in the argument list. The function \eqn{x(F)} is computed through the \code{\link{par2qua}} function. The distribution type is determined using the \code{type} attribute of the \code{para} argument, which is a parameter object of \pkg{lmomco} (see \code{\link{vec2par}}).
}
\usage{
theopwms(para, nmom=5, minF=0, maxF=1, quafunc=NULL,
               nsim=50000, fold=5,
               silent=TRUE, verbose=FALSE, ...)
}
\arguments{
  \item{para}{A distribution parameter object such as that by \code{\link{lmom2par}} or \code{\link{vec2par}}.}
  \item{nmom}{The number of moments to compute. Default is 5.}
  \item{minF}{The end point of nonexceedance probability in which to perform the integration. Try setting to non-zero (but small) if you have a divergent integral.}
  \item{maxF}{The end point of nonexceedance probability in which to perform the integration. Try setting to non-unity (but close) if you have a divergent integral.}
  \item{quafunc}{An optional and arbitrary quantile function that simply needs to except a nonexceedance probability and the parameter object in \code{para}. This is a feature that permits computation of the PWMs of a quantile function that does not have to be implemented in the greater overhead hassles of the \pkg{lmomco} style. This feature might be useful for estimation of quantile function mixtures or those distributions not otherwise implemented in this package.}
   \item{nsim}{Simulation size for Monte Carlo integration is such is internally deemed necessary (see \code{silent} argument).}
   \item{fold}{The number of fractions or number of folds of \code{nsim}, which in other words, means that \code{nsim} is divided by \code{folds} and a loop creating \code{folds} integrations of \code{nsim/folds} is used from which the mean and mean absolute error of the integrand are computed. This is to try to recover similar output as \code{integrate()}.}
   \item{silent}{The argument of \code{silent} for the \code{try()} operation wrapped on \code{integrate()}. If set true and the integral is probability divergent, Monte Carlo integration is triggered using \code{nsim} and \code{folds}. The user would have to set \code{verbose=TRUE} to then acquire the returned table in \code{integrations} of the integration passes including those are or are not Monte Carlo.}
   \item{verbose}{Toggle verbose output. Because the \R function \code{integrate} is used to perform the numerical integration, it might be useful to see selected messages regarding the numerical integration.}
   \item{...}{Additional arguments to pass.}
}
\value{
  An \R \code{list} is returned.

  \item{betas}{The PWMs. Note that convention is the have a \eqn{\beta_0}, but this is placed in the first index \code{i=1} of the \code{betas} vector.}
  \item{nsim}{Echo of the \code{nsim} argument if and only if at least one Monte Carlo integration was required, otherwise this is set to \dQuote{not needed} on the return.}
  \item{folds}{Echo of the \code{folds} argument if and only if at least one Monte Carlo integration was required, otherwise this is set to \dQuote{not needed} on the return.}
  \item{monte_carlo}{A logical vector of whether one or more Monte Carlo integrations was needed for the \code{r}-th index of the vector during the integrations for the \eqn{r}-th PWM (beta).}
  \item{source}{An attribute identifying the computational source of the probability-weighted moments: \dQuote{theopwms}.}
  \item{integrations}{If \code{verbose=TRUE}, then the results of the integrations are a data frame stored here. Otherwise, \code{integrations} is not present in the list.}
}
\references{
Hosking, J.R.M., 1990, L-moments--Analysis and estimation of distributions using linear combinations of order statistics: Journal of the Royal Statistical Society, Series B, v. 52, p. 105--124.
}
\author{W.H. Asquith}

\seealso{\code{\link{theoLmoms}}, \code{\link{pwm}}, \code{\link{pwm2lmom}}}
\examples{
para     <- vec2par(c(0,1),type='nor') # standard normal
the.pwms <- theopwms(para) # compute PWMs
str(the.pwms)

\dontrun{
  # This example has a divergent integral triggered on the beta0. Monte Carlo (MC)
  # integration is thus triggered. The verbose=TRUE saves numerical or MC
  # integration result table to the return.
  para <- vec2par(c(2,2, 1.8673636098392308, -0.1447286792099476), type="kap")
  pwmkap <- lmom2pwm( lmomkap(para) )
  print(pwmkap$betas) # 0.1155903 1.2153105 0.9304619 0.7282926 0.5938137
  pwmthe <- theopwms(para, nmom=5, verbose=TRUE)
  print(pwmthe$betas) # 0.1235817 1.2153104 0.9304619 0.7282926 0.5938137

  para <- vec2par(c(2,2, 0.9898362024687231, -0.5140894097276032), type="kap")
  pwmkap <- lmom2pwm( lmomkap(para) )
  print(pwmkap$betas) # -0.06452787  1.33177963  1.06818379  0.85911124  0.71308145
  pwmthe <- theopwms(para, nmom=5, verbose=TRUE)
  print(pwmthe$betas) # -0.06901669  1.33177952  1.06818379  0.85911123  0.71308144 }
}
\keyword{probability-weighted moment (theoretical)}

