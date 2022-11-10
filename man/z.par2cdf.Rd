\name{z.par2cdf}
\alias{z.par2cdf}
\title{Blipping Cumulative Distribution Functions}
\description{
This function acts as a front end or dispatcher to the distribution-specific cumulative distribution functions but also provides for blipping according to
\deqn{F(x) = 0}
for \eqn{x \le z} and
\deqn{F(x) = p + (1-p)G(x)}
for \eqn{x > z} where \eqn{z} is a threshold value. The \eqn{z} is not tracked as part of the parameter object. This might arguably be a design flaw, but the function will do its best to test whether the \eqn{z} given is compatable (but not necessarily equal to \eqn{\hat{x} = x(0)}) with the quantile function \eqn{x(F)} (\code{\link{z.par2qua}}). Lastly, please refer to the finiteness check in the Examples to see how one might accommodate \eqn{-\infty} for \eqn{F = 0} on a standard normal variate plot.

A recommended practice when working with this function is the insertion of the \eqn{x} value at \eqn{F=p}. Analogous practice is suggested for \code{\link{z.par2qua}} (see that documentation).
}
\usage{
z.par2cdf(x, p, para, z=0, ...)
}
\arguments{
  \item{x}{A real value vector.}
  \item{p}{Nonexceedance probability of the \code{z} value. This probability could simply be the portion of record having zero values if \code{z=0}.}
  \item{para}{The parameters from \code{\link{lmom2par}} or \code{\link{vec2par}}.}
  \item{z}{Threshold value.}
  \item{...}{The additional arguments are passed to the cumulative distribution function such as \code{paracheck=FALSE} for the Generalized Lambda distribution (\code{\link{cdfgld}}).}
}
\value{
  Nonexceedance probability (\eqn{0 \le F \le 1}) for \code{x}.
}
\author{W.H. Asquith}
\references{
Asquith, W.H., 2011, Distributional analysis with L-moment statistics using the R environment for statistical computing: Createspace Independent Publishing Platform, ISBN 978--146350841--8.
}
\seealso{\code{\link{z.par2qua}}, \code{\link{par2cdf}}}
\examples{
set.seed(21)
the.gpa   <- vec2par(c(100,1000,0.1),type='gpa')
fake.data <- rlmomco(30,the.gpa) # simulate some data
fake.data <- sort(c(fake.data,rep(0,10))) # add some zero observations
# going to tick to the inside and title right axis as well, so change some
# plotting parameters
par(mgp=c(3,0.5,0), mar=c(5,4,4,3))
# next compute the parameters for the positive data
gpa.all <- pargpa(lmoms(fake.data))
gpa.nzo <- pargpa(lmoms(fake.data[fake.data > 0]))
n   <- length(fake.data) # sample size
p   <- length(fake.data[fake.data == 0])/n # est. prob of zero value
F   <- nonexceeds(sig6=TRUE); F <- sort(c(F,p)); qF <- qnorm(F)
# The following x vector obviously contains zero, so no need to insert it.
x   <- seq(-100, max(fake.data)) # absurd for x<0, but testing implementation
PP  <- pp(fake.data) # compute plotting positions of sim. sample
plot(fake.data, qnorm(PP), xlim=c(0,4000), yaxt="n", ylab="") # plot the sample
add.lmomco.axis(las=2, tcl=0.5, side=2, twoside=FALSE,
                                        side.type="NPP", otherside.type="SNV")
lines(quagpa(F,gpa.all), qF) # the parent (without zeros)
cdf <- qnorm(z.par2cdf(x,p,gpa.nzo))
cdf[! is.finite(cdf)] <- min(fake.data,qnorm(PP)) # See above documentation
lines(x, cdf,lwd=3) # fitted model with zero conditional
# now repeat the above code over and over again and watch the results
par(mgp=c(3,1,0), mar=c(5,4,4,2)+0.1) # restore defaults
}
\keyword{cumulative distribution function}
\keyword{blipping distributions}

