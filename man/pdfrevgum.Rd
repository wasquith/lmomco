\name{pdfrevgum}
\alias{pdfrevgum}
\title{Probability Density Function of the Reverse Gumbel Distribution}
\description{
This function computes the probability density
of the Reverse Gumbel distribution given parameters (\eqn{\xi} and \eqn{\alpha})  computed by \code{\link{parrevgum}}. The probability density function is
\deqn{f(x) = \alpha^{-1} \exp(Y) [\exp(\exp[-\exp(Y)])] \mbox{,}}
where
\deqn{Y = \frac{x - \xi}{\alpha} \mbox{,}}
where \eqn{f(x)} is the probability density for quantile \eqn{x},
\eqn{\xi} is a location parameter, and \eqn{\alpha} is a scale parameter.
}
\usage{
pdfrevgum(x, para)
}
\arguments{
  \item{x}{A real value vector.}
  \item{para}{The parameters from \code{\link{parrevgum}} or \code{\link{vec2par}}.}
}

\value{
  Probability density (\eqn{f}) for \eqn{x}.
}
\references{
Hosking, J.R.M., 1990, L-moments---Analysis and estimation of
distributions using linear combinations of order statistics: Journal
of the Royal Statistical Society, Series B, v. 52, pp. 105--124.

Hosking, J.R.M., 1995, The use of L-moments in the analysis of censored data,
in Recent Advances in Life-Testing and Reliability, edited by N. Balakrishnan, chapter 29, CRC Press, Boca Raton, Fla., pp. 546--560.
}
\author{ W.H. Asquith}
\seealso{\code{\link{cdfrevgum}}, \code{\link{quarevgum}}, \code{\link{lmomrevgum}}, \code{\link{parrevgum}}}
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
# lamB1 = -0.516 and lamB2 = 0.523
btypelmoms <- pwm2lmom(z$Bbetas)
# My version of R reports lamB1 = -0.5162 and lamB2 = 0.5218
str(btypelmoms)
rg.pars <- parrevgum(btypelmoms,z$zeta)
str(rg.pars)
# Hosking reports xi=0.1636 and alpha=0.9252 for the sample
# My version of R reports xi = 0.1635 and alpha = 0.9254
# Now one can continue one with a plotting example.
\dontrun{
F  <- nonexceeds()
PP <- pp(D) # plotting positions of the data
D  <- sort(D)
plot(D,PP)
lines(D,cdfrevgum(D,rg.pars))
# Now finally do the PDF
F <- seq(0.01,0.99,by=.01)
x <- quarevgum(F,rg.pars)
plot(x,pdfrevgum(x,rg.pars),type='l')
}
}
\keyword{distribution}
\keyword{probability density function}
\keyword{Distribution: Gumbel (reversed)}
