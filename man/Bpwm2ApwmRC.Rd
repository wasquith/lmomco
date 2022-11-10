\name{Bpwm2ApwmRC}
\alias{Bpwm2ApwmRC}
\title{Conversion between B- and A-Type Probability-Weighted Moments for Right-Tail Censoring of an Appropriate Distribution}
\description{
This function converts \dQuote{B}-type probability-weighted moments (PWMs, \eqn{\beta^B_r}) to the \dQuote{A}-type \eqn{\beta^A_r}. The \eqn{\beta^A_r} are the ordinary PWMs for the \eqn{m} left noncensored or observed values. The \eqn{\beta^B_r} are more complex and use the \eqn{m} observed values and the \eqn{m-n} right-tailed censored values for which the censoring threshold is known. The \dQuote{A}- and \dQuote{B}-type PWMs are described in the documentation for \code{\link{pwmRC}}.

This function uses the defined relation between to two PWM types when the \eqn{\beta^B_r} are known along with the parameters (\code{para}) of a right-tail censored distribution inclusive of the censoring fraction \eqn{\zeta=m/n}. The value \eqn{\zeta} is the right-tail censor fraction or the probability \eqn{\mathrm{Pr}\lbrace \rbrace} that \eqn{x} is less than the quantile at \eqn{\zeta} nonexceedance probability (\eqn{\mathrm{Pr}\lbrace x < X(\zeta) \rbrace}). The relation is
\deqn{\beta^A_{r-1} = \frac{r\beta^B_{r-1} - (1-\zeta^r)X(\zeta)}{r\zeta^r} \mbox{,}}
where \eqn{1 \le r \le n} and \eqn{n} is the number of moments, and \eqn{X(\zeta)} is the value of the quantile function at nonexceedance probability \eqn{\zeta}. Finally, the \code{RC} in the function name is to denote \code{R}ight-tail \code{C}ensoring.
}
\usage{
Bpwm2ApwmRC(Bpwm,para)
}
\arguments{
  \item{Bpwm}{A vector of B-type PWMs: \eqn{\beta^B_r}.}
  \item{para}{The parameters of the distribution from a function such as \code{pargpaRC} in which the \eqn{\beta^B_r} are contained in a \code{list} element titled \code{betas} and the right-tail censoring fraction \eqn{\zeta} is contained in an element titled \code{zeta}.}
}
\value{
  An \R \code{list} is returned.

}
\references{
Hosking, J.R.M., 1995, The use of L-moments in the analysis of censored data,
in Recent Advances in Life-Testing and Reliability, edited by N. Balakrishnan,
chapter 29, CRC Press, Boca Raton, Fla., pp. 546--560.
}
\author{W.H. Asquith}
\seealso{\code{\link{Apwm2BpwmRC}} and \code{\link{pwmRC}} }
\examples{
# Data listed in Hosking (1995, table 29.2, p. 551)
H <- c(3,4,5,6,6,7,8,8,9,9,9,10,10,11,11,11,13,13,13,13,13,
             17,19,19,25,29,33,42,42,51.9999,52,52,52)
      # 51.9999 was really 52, a real (noncensored) data point.
z <-  pwmRC(H,52)
# The B-type PMWs are used for the parameter estimation of the
# Reverse Gumbel distribution. The parameter estimator requires
# conversion of the PWMs to L-moments by pwm2lmom().
para <- parrevgum(pwm2lmom(z$Bbetas),z$zeta) # parameter object
Abetas <- Bpwm2ApwmRC(z$Bbetas,para)
Bbetas <- Apwm2BpwmRC(Abetas$betas,para)
# Assertion that both of the vectors of B-type PWMs should be the same.
str(Bbetas)   # B-type PWMs of the distribution
str(z$Bbetas) # B-type PWMs of the original data
}
\keyword{probability-weighted moment (theoretical)}
\keyword{data censoring}
