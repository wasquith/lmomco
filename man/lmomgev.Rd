\name{lmomgev}
\alias{lmomgev}
\title{L-moments of the Generalized Extreme Value Distribution}
\description{
This function estimates the L-moments of the Generalized Extreme Value distribution given the parameters
(\eqn{\xi}, \eqn{\alpha}, and \eqn{\kappa}) from \code{\link{pargev}}. The L-moments in terms of the parameters are
\deqn{\lambda_1 = \xi + \frac{\alpha}{\kappa}(1-\Gamma(1+\kappa)) \mbox{,}}
\deqn{\lambda_2 = \frac{\alpha}{\kappa}(1-2^{-\kappa})\Gamma(1+\kappa) \mbox{,}}
\deqn{\tau_3 = \frac{2(1-3^{-\kappa})}{1-2^{-\kappa}} - 3 \mbox{, and}}
\deqn{\tau_4 = \frac{5(1-4^{-\kappa})-10(1-3^{-\kappa})+6(1-2^{-\kappa})}{1-2^{-\kappa}} \mbox{.}}
}
\usage{
lmomgev(para)
}
\arguments{
  \item{para}{The parameters of the distribution.}
}
\value{
  An \R \code{list} is returned.

  \item{lambdas}{Vector of the L-moments. First element is
\eqn{\lambda_1}, second element is \eqn{\lambda_2}, and so on.}
  \item{ratios}{Vector of the L-moment ratios. Second element is
\eqn{\tau}, third element is \eqn{\tau_3} and so on.}
  \item{trim}{Level of symmetrical trimming used in the computation, which is \code{0}.}
  \item{leftrim}{Level of left-tail trimming used in the computation, which is \code{NULL}.}
  \item{rightrim}{Level of right-tail trimming used in the computation, which is \code{NULL}.}
  \item{source}{An attribute identifying the computational
source of the L-moments: \dQuote{lmomgev}.}
}
\references{
Hosking, J.R.M., 1990, L-moments---Analysis and estimation of distributions using linear combinations of order statistics: Journal of the Royal Statistical Society, Series B, v. 52, pp. 105--124.

Hosking, J.R.M., 1996, FORTRAN routines for use with the method of L-moments: Version 3, IBM Research Report RC20525, T.J. Watson Research Center, Yorktown Heights, New York.

Hosking, J.R.M., and Wallis, J.R., 1997, Regional frequency analysis---An approach based on L-moments: Cambridge University Press.
}
\author{W.H. Asquith}
\seealso{\code{\link{pargev}}, \code{\link{cdfgev}}, \code{\link{pdfgev}}, \code{\link{quagev}} }
\examples{
lmr <- lmoms(c(123,34,4,654,37,78))
lmomgev(pargev(lmr))

\dontrun{
# The Gumbel is a limiting version of the maxima regardless of parent. The GLO,
# PE3 (twice), and GPA are studied here. A giant number of events to simulate is made.
# Then numbers of events per year before the annual maxima are computed are specified.
# The Gumbel is a limiting version of the maxima regardless of parent. The GLO,
# PE3 (twice), and GPA are studied here. A giant number of events to simulate is made.
# Then numbers of events per year before the annual maxima are computed are specified.
nevents <- 100000
nev_yr <- c(1,2,3,4,5,6,10,15,20,30,50,100,200,500); n <- length(nev_yr)
pdf("Gumbel_in_the_limit.pdf", useDingbats=FALSE)
# Draw the usually L-moment ratio diagram but only show a few of the
# three parameter families.
plotlmrdia(lmrdia(), xlim=c(-.5,0.5), ylim=c(0,0.3), nopoints=TRUE,
           autolegend=TRUE, noaep4=TRUE, nogov=TRUE, xleg=0.1, yleg=0.3)
gum <- lmrdia()$gum # extract the L-skew and L-kurtosis of the Gumbel
points(gum[1], gum[2], pch=10, cex=3, col=2) # draw the Gumbel

para <- parglo(vec2lmom(c(1,.1,0))) # generalized logistic
t3 <- t4 <- rep(NA, n) # define
for(k in 1:n) { # generate GLO time series of annual maxima with k-events per year
   lmr <- lmoms(sapply(1:nevents/nev_yr[k], function(i) max(rlmomco(nev_yr[k], para))))
   t3[k] <- lmr$ratios[3]; t4[k] <- lmr$ratios[4]
}
lines(t3, t4, lwd=0.8); points(t3, t4, lwd=0.8, pch=21, bg=3)

para <- parglo(vec2lmom(c(1,.1,0.3))) # generalized logistic
t3 <- t4 <- rep(NA, n) # define
for(k in 1:n) { # generate GLO time series of annual maxima with k-events per year
   lmr <- lmoms(sapply(1:nevents/nev_yr[k], function(i) max(rlmomco(nev_yr[k], para))))
   t3[k] <- lmr$ratios[3]; t4[k] <- lmr$ratios[4]
}
lines(t3, t4, lwd=0.8); points(t3, t4, lwd=0.8, pch=21, bg=3)

para <- parglo(vec2lmom(c(1,.1,-0.3))) # generalized logistic
t3 <- t4 <- rep(NA, n) # define
for(k in 1:n) { # generate GLO time series of annual maxima with k-events per year
   lmr <- lmoms(sapply(1:nevents/nev_yr[k], function(i) max(rlmomco(nev_yr[k], para))))
   t3[k] <- lmr$ratios[3]; t4[k] <- lmr$ratios[4]
}
lines(t3, t4, lwd=0.8); points(t3, t4, lwd=0.8, pch=21, bg=3)

para <- parpe3(vec2lmom(c(1,.1,.4))) # Pearson type III
t3 <- t4 <- rep(NA, n) # reset
for(k in 1:n) { # generate PE3 time series of annual maxima with k-events per year
   lmr <- lmoms(sapply(1:nevents/k, function(i) max(rlmomco(nev_yr[k], para))))
   t3[k] <- lmr$ratios[3]; t4[k] <- lmr$ratios[4]
}
lines(t3, t4, lwd=0.8); points(t3, t4, lwd=0.8, pch=21, bg=6)

para <- parpe3(vec2lmom(c(1,.1,0))) # Pearson type III
t3 <- t4 <- rep(NA, n) # reset
for(k in 1:n) { # generate another PE3 time series of annual maxima with k-events per year
   lmr <- lmoms(sapply(1:nevents/k, function(i) max(rlmomco(nev_yr[k], para))))
   t3[k] <- lmr$ratios[3]; t4[k] <- lmr$ratios[4]
}
lines(t3, t4, lwd=0.8); points(t3, t4, lwd=0.8, pch=21, bg=6)

para <- parpe3(vec2lmom(c(1,.1,-.4))) # Pearson type III
t3 <- t4 <- rep(NA, n) # reset
for(k in 1:n) { # generate PE3 time series of annual maxima with k-events per year
   lmr <- lmoms(sapply(1:nevents/k, function(i) max(rlmomco(nev_yr[k], para))))
   t3[k] <- lmr$ratios[3]; t4[k] <- lmr$ratios[4]
}
lines(t3, t4, lwd=0.8); points(t3, t4, lwd=0.8, pch=21, bg=6)

para <- pargpa(vec2lmom(c(1,.1,0))) # generalized Pareto
t3 <- t4 <- rep(NA, n) # reset
for(k in 1:n) { # generate GPA time series of annual maxima with k-events per year
   lmr <- lmoms(sapply(1:nevents/k, function(i) max(rlmomco(nev_yr[k], para))))
   t3[k] <- lmr$ratios[3]; t4[k] <- lmr$ratios[4]
}
lines(t3, t4, lwd=0.8); points(t3, t4, lwd=0.8, pch=21, bg=4)

para <- pargpa(vec2lmom(c(1,.1,.4))) # generalized Pareto
t3 <- t4 <- rep(NA, n) # reset
for(k in 1:n) { # generate GPA time series of annual maxima with k-events per year
   lmr <- lmoms(sapply(1:nevents/k, function(i) max(rlmomco(nev_yr[k], para))))
   t3[k] <- lmr$ratios[3]; t4[k] <- lmr$ratios[4]
}
lines(t3, t4, lwd=0.8); points(t3, t4, lwd=0.8, pch=21, bg=4)

para <- pargpa(vec2lmom(c(1,.1,-.4))) # generalized Pareto
t3 <- t4 <- rep(NA, n) # reset
for(k in 1:n) { # generate GPA time series of annual maxima with k-events per year
   lmr <- lmoms(sapply(1:nevents/k, function(i) max(rlmomco(nev_yr[k], para))))
   t3[k] <- lmr$ratios[3]; t4[k] <- lmr$ratios[4]
}
lines(t3, t4, lwd=0.8); points(t3, t4, lwd=0.8, pch=21, bg=4)
dev.off() #}
}
\keyword{L-moment (distribution)}
\keyword{Distribution: Generalized Extreme Value}
