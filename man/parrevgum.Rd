\name{parrevgum}
\alias{parrevgum}
\title{Estimate the Parameters of the Reverse Gumbel Distribution}
\description{
This function estimates the parameters of the Reverse Gumbel distribution given the type-B L-moments of the data in an L-moment object such as that returned by
\code{\link{pwmRC}} using \code{\link{pwm2lmom}}. This distribution is important in the analysis of censored data. It is the distribution of a logarithmically transformed 2-parameter Weibull distribution. The relations between distribution parameters and L-moments are
\deqn{\alpha = \lambda^B_2/\lbrace\log(2) + \mathrm{Ei}(-2\log(1-\zeta)) - \mathrm{Ei}(-\log(1-\zeta))\rbrace}
and
\deqn{\xi = \lambda^B_1 + \alpha\lbrace\mathrm{Ei}(-\log(1-\zeta))\rbrace\mbox{,}}
where \eqn{\zeta} is the compliment of the right-tail censoring fraction of the sample or the nonexceedance probability of the right-tail censoring threshold, and \eqn{\mathrm{Ei}(x)} is the exponential integral defined as
\deqn{ \mathrm{Ei}(X) = \int_X^{\infty} x^{-1}e^{-x}\mathrm{d}x \mbox{,}}
where \eqn{\mathrm{Ei}(-\log(1-\zeta)) \rightarrow 0} as \eqn{\zeta \rightarrow 1} and \eqn{\mathrm{Ei}(-\log(1-\zeta))} can not be evaluated as \eqn{\zeta \rightarrow 0}.
}
\usage{
parrevgum(lmom, zeta=1, checklmom=TRUE, ...)
}
\arguments{
  \item{lmom}{An L-moment object created by \code{\link{lmoms}} through \code{\link{pwmRC}} or other L-moment type object. The user intervention of the \code{zeta} differentiates this distribution (and this function) from similar parameter estimation functions in the \pkg{lmomco} package.}
  \item{zeta}{The compliment of the right censoring fraction. Number of samples observed (noncensored) divided by the total number of samples.}
  \item{checklmom}{Should the \code{lmom} be checked for validity using the \code{\link{are.lmom.valid}} function. Normally this should be left as the default and it is very unlikely that the L-moments will not be viable (particularly in the \eqn{\tau_4} and \eqn{\tau_3} inequality). However, for some circumstances or large simulation exercises then one might want to bypass this check.}
  \item{...}{Other arguments to pass.}
}

\value{
  An \R \code{list} is returned.

  \item{type}{The type of distribution: \code{revgum}.}
  \item{para}{The parameters of the distribution.}
  \item{zeta}{The compliment of the right censoring fraction. Number of samples observed (noncensored) divided by the total number of samples.}
  \item{source}{The source of the parameters: \dQuote{parrevgum}.}
}
\references{
Hosking, J.R.M., 1995, The use of L-moments in the analysis of censored data, in Recent Advances in Life-Testing and Reliability, edited by N. Balakrishnan, chapter 29, CRC Press, Boca Raton, Fla., pp. 546--560.
}
\author{W.H. Asquith}
\seealso{\code{\link{lmomrevgum}},
\code{\link{cdfrevgum}}, \code{\link{pdfrevgum}}, \code{\link{quarevgum}}, \code{\link{pwm2lmom}}, \code{\link{pwmRC}}
}
\examples{
# See p. 553 of Hosking (1995)
# Data listed in Hosking (1995, table 29.3, p. 553)
D <- c(-2.982, -2.849, -2.546, -2.350, -1.983, -1.492, -1.443,
       -1.394, -1.386, -1.269, -1.195, -1.174, -0.854, -0.620,
       -0.576, -0.548, -0.247, -0.195, -0.056, -0.013,  0.006,
        0.033,  0.037,  0.046,  0.084,  0.221,  0.245,  0.296)
D <- c(D,rep(.2960001,40-28)) # 28 values, but Hosking mentions
                              # 40 values in total
z <-  pwmRC(D,threshold=.2960001)
str(z)
# Hosking reports B-type L-moments for this sample are
# lamB1 = -.516 and lamB2 = 0.523
btypelmoms <- pwm2lmom(z$Bbetas)
# My version of R reports lamB1 = -0.5162 and lamB2 = 0.5218
str(btypelmoms)
rg.pars <- parrevgum(btypelmoms,z$zeta)
str(rg.pars)
# Hosking reports xi = 0.1636 and alpha = 0.9252 for the sample
# My version of R reports xi = 0.1635 and alpha = 0.9254
}
\keyword{distribution (parameters)}
\keyword{Distribution: Gumbel (reversed)}
