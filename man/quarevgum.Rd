\name{quarevgum}
\alias{quarevgum}
\title{Quantile Function of the Reverse Gumbel Distribution}
\description{
This function computes the quantiles of the Reverse Gumbel distribution given parameters (\eqn{\xi} and \eqn{\alpha}) computed by \code{\link{parrevgum}}. The quantile function  is
\deqn{x(F) = \xi + \alpha\log(-\log(1-F)) \mbox{,}}
where \eqn{x(F)} is the quantile for nonexceedance probability \eqn{F}, \eqn{\xi} is a location parameter, and \eqn{\alpha} is a scale parameter.
}
\usage{
quarevgum(f, para, paracheck=TRUE)
}
\arguments{
  \item{f}{Nonexceedance probability (\eqn{0 \le F \le 1}).}
  \item{para}{The parameters from \code{\link{parrevgum}} or \code{\link{vec2par}}.}
  \item{paracheck}{A logical controlling whether the parameters are checked for validity. Overriding of this check might be extremely important and needed for use of the quantile function in the context of TL-moments with nonzero trimming.}
}
\value{
  Quantile value for nonexceedance probability \eqn{F}.
}
\references{
Hosking, J.R.M., 1995, The use of L-moments in the analysis of censored data,
in Recent Advances in Life-Testing and Reliability, edited by N. Balakrishnan,
chapter 29, CRC Press, Boca Raton, Fla., pp. 546--560.
}
\author{ W.H. Asquith}
\seealso{\code{\link{cdfrevgum}}, \code{\link{pdfrevgum}}, \code{\link{lmomrevgum}}, \code{\link{parrevgum}}}
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
F  <- nonexceeds()
PP <- pp(D) # plotting positions of the data
plot(PP,sort(D),ylim=range(quarevgum(F,rg.pars)))
lines(F,quarevgum(F,rg.pars))
# In the plot notice how the data flat lines at the censoring level,
# but the distribution continues on.  Neat.
}
\keyword{distribution}
\keyword{quantile function}
\keyword{Distribution: Gumbel (reversed)}

