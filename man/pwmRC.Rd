\name{pwmRC}
\alias{pwmRC}
\title{Sample Probability-Weighted Moments for Right-Tail Censoring }
\description{
Compute the sample Probability-Weighted Moments (PWMs) for right-tail censored data set---that is a data set censored from above. The censoring threshold is denoted as \eqn{T}. The data possess \eqn{m} values that are observed (noncensored, \eqn{< T}) out of a total of \eqn{n} samples. The ratio of \eqn{m} to \eqn{n} is defined as \eqn{\zeta = m/n}, which will play an important role in parameter estimation. The \eqn{\zeta} is interpreted as the probability \eqn{\mathrm{Pr}\lbrace \rbrace} that \eqn{x} is less than the quantile at \eqn{\zeta} nonexceedance probability: (\eqn{\mathrm{Pr}\lbrace x < X(\zeta) \rbrace}). Two types of PWMs are computed for right-tail censored situations. The \dQuote{A}-type PWMs and \dQuote{B}-type PWMs. The A-type PWMs are defined by

\deqn{\beta^A_r = m^{-1}\sum^m_{j=1} {j-1 \choose r} x_{[j:n]}\mbox{,}}

which are the PWMs of the uncensored sample of \eqn{m} observed values. The B-type PWMs are computed from the \dQuote{complete} sample, in which the \eqn{n-m} censored values are replaced by the censoring threshold \eqn{T}. The B-type PWMs are defined by

\deqn{\beta^B_r = n^{-1} \biggl( \sum^m_{j=1} {j-1 \choose r} x_{[j:n]} +
                                \sum^n_{j=m+1} {j-1 \choose r} T \biggr) \mbox{.}}

The two previous expressions are used in the function. These PWMs are readily converted to L-moments by the usual methods (\code{\link{pwm2lmom}}). When there are more than a few censored values, the PWMs are readily computed by computing \eqn{\beta^A_r} and using the expression

\deqn{\beta^B_r = Z\beta^A_r + \frac{1-Z}{r+1}T\mbox{,}}

where

\deqn{Z = \frac{m}{n}\frac{{m-1 \choose r}}{{n-1 \choose r}}\mbox{.}}

The two expressions above are consulted when the \code{checkbetas=TRUE} argument is present. Both sequences of B-type are \code{cat}ed to the terminal. This provides a check on the implementation of the algorithm. The functions \code{\link{Apwm2BpwmRC}} and \code{\link{Bpwm2ApwmRC}} can be used to switch back and forth between the two PWM types given fitted parameters for a distribution in the \pkg{lmomco} package that supports right-tail censoring. Finally, the \code{RC} in the function name is to denote \code{R}ight-tail \code{C}ensoring.
}
\usage{
pwmRC(x, threshold=NULL, nmom=5, sort=TRUE, checkbetas=FALSE)
}
\arguments{
  \item{x}{A vector of data values.}
  \item{threshold}{The right-tail censoring (upper) threshold.}
  \item{nmom}{Number of PWMs to return.}
  \item{sort}{Do the data need sorting? Note that convention is the have a \eqn{\beta_0}, but this is placed in the first index \code{i=1} of the \code{betas} vector.}
  \item{checkbetas}{A cross relation between \eqn{\beta^A_r} and \eqn{\beta^B_r} exists---display the results of the secondary computation of the \eqn{\beta^B_r}. The two displayed vectors should be numerically equal.}
}
\value{
  An \R \code{list} is returned.

  \item{Abetas}{The A-type PWMs. These should be same as \code{pwm()} returns if there is no censoring. Note that convention is the have a \eqn{\beta_0}, but this is placed in the first index \code{i=1} of the \code{betas} vector.}
  \item{Bbetas}{The B-type PWMs. These should be \code{NA} if there is no censoring. Note that convention is the have a \eqn{\beta_0}, but this is placed in the first index \code{i=1} of the \code{betas} vector.}
  \item{source}{Source of the PWMs: \dQuote{pwmRC}.}
  \item{threshold}{The upper censoring threshold.}
  \item{zeta}{The right censoring fraction: \code{numabovethreshold/samplesize}.}
  \item{numabovethreshold}{Number of data points equal to or above the threshold.}
  \item{observedsize}{Number of real data points in the sample (below the threshold).}
  \item{samplesize}{Number of actual sample values.}
}
\details{
There is some ambiguity if the threshold also numerically equals valid data in the data set. In the data for the examples below, which are taken from elsewhere, there are real observations at the censoring level. One can see how a hack is made to marginally decrease or increase the data or the threshold for the computations. This is needed because the code uses
\preformatted{
sapply(x, function(v) { if(v >= T) return(T); return(v) } )
}
to reset the data vector \code{x}. By operating on the data in this fashion one can toy with various levels of the threshold for experimental purposes; this seemed a more natural way for general implementation. The code sets \eqn{n} = \code{length(x)} and \eqn{m} = \code{n - length(x[x == T])}, which also seems natural. The \eqn{\beta^A_r} are computed by dispatching to \code{\link{pwm}}.
}
\references{
Greenwood, J.A., Landwehr, J.M., Matalas, N.C., and Wallis, J.R., 1979, Probability weighted moments---Definition and relation to parameters of several distributions expressable in inverse form: Water Resources Research, v. 15, pp. 1,049--1,054.

Hosking, J.R.M., 1990, L-moments---Analysis and estimation of distributions using linear combinations of order statistics: Journal of the Royal Statistical Society, Series B, v. 52, pp. 105--124.

Hosking, J.R.M., 1995, The use of L-moments in the analysis of censored data, in Recent Advances in Life-Testing and Reliability, edited by N. Balakrishnan, chapter 29, CRC Press, Boca Raton, Fla., pp. 546--560.
}
\author{W.H. Asquith}
\seealso{\code{\link{lmoms}}, \code{\link{pwm2lmom}}, \code{\link{pwm}}, \code{\link{pwmLC}}  }
\examples{
# Data listed in Hosking (1995, table 29.2, p. 551)
H <- c(3,4,5,6,6,7,8,8,9,9,9,10,10,11,11,11,13,13,13,13,13,
       17,19,19,25,29,33,42,42,51.9999,52,52,52)
# 51.9999 was really 52, a real (noncensored) data point.
z <-  pwmRC(H,threshold=52,checkbetas=TRUE)
str(z)
# Hosking(1995) reports that A-type L-moments for this sample are
# lamA1=15.7 and lamAL-CV=.389, and lamAL-skew=.393
pwm2lmom(z$Abetas)
# My version of R reports 15.666, 0.3959, and 0.4030


# See p. 553 of Hosking (1995)
# Data listed in Hosking (1995, table 29.3, p. 553)
D <- c(-2.982, -2.849, -2.546, -2.350, -1.983, -1.492, -1.443,
       -1.394, -1.386, -1.269, -1.195, -1.174, -0.854, -0.620,
       -0.576, -0.548, -0.247, -0.195, -0.056, -0.013,  0.006,
        0.033,  0.037,  0.046,  0.084,  0.221,  0.245,  0.296)
D <- c(D,rep(.2960001,40-28)) # 28 values, but Hosking mentions
                              # 40 values in total
z <-  pwmRC(D,.2960001)
# Hosking reports B-type L-moments for this sample are
# lamB1 = -.516 and lamB2 = 0.523
pwm2lmom(z$Bbetas)
# My version of R reports -.5162 and 0.5218
}
\keyword{probability-weighted moment (sample)}
\keyword{data censoring}
