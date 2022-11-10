\name{z.par2qua}
\alias{z.par2qua}
\title{Blipping Quantile Functions}
\description{
This function acts as a front end or dispatcher to the distribution-specific quantile functions but also provides for blipping for zero (or other) threshold according to
\deqn{x(F) = 0}
for \eqn{0 \le F \le p} and
\deqn{x_G\left(\frac{F-p}{1-p}\right)}
for \eqn{F > p}. This function is generalized for \eqn{z \ne 0}. The \eqn{z} is not tracked as part of the parameter object. This might arguably be a design flaw, but the function will do its best to test whether the \eqn{z} given is compatable (but not necessarily equal to \eqn{\hat{x} = x(0)}) with the quantile function \eqn{x(F)}.

A recommended practice when working with this function when \eqn{F} values are generated for various purposes, such as for graphics, then the value of \eqn{p} should be inserted into the vector, and the vector obviously sorted (see the line using the \code{\link{nonexceeds}} function). This should be considered as well when \code{\link{z.par2cdf}} is used but with the insertion of the \eqn{x} value at \eqn{F=p}.
}
\usage{
z.par2qua(f, p, para, z=0, ...)
}
\arguments{
  \item{f}{Nonexceedance probabilities (\eqn{0 \le F \le 1}).}
  \item{p}{Nonexceedance probability of \code{z} value.}
  \item{para}{The parameters from \code{\link{lmom2par}} or \code{\link{vec2par}}.}
  \item{z}{Threshold value.}
  \item{...}{The additional arguments are passed to the quantile function such as \cr  \code{paracheck = FALSE} for the Generalized Lambda distribution (\code{\link{quagld}}).}
}
\value{
  Quantile value for \eqn{f}.
}
\author{W.H. Asquith}
\references{
Asquith, W.H., 2011, Distributional analysis with L-moment statistics using the R environment for statistical computing: Createspace Independent Publishing Platform, ISBN 978--146350841--8.
}
\seealso{\code{\link{z.par2cdf}}, \code{\link{par2qua}}}
\examples{
# define the real parent (or close)
the.gpa   <- vec2par(c(100,1000,0.1),type='gpa')
fake.data <- rlmomco(30,the.gpa) # simulate some data
fake.data <- sort(c(fake.data, rep(0,10))) # add some zero observations

par(mgp=c(3,0.5,0)) # going to tick to the inside, change some parameters
# next compute the parameters for the positive data
gpa.all <- pargpa(lmoms(fake.data))
gpa.nzo <- pargpa(lmoms(fake.data[fake.data > 0]))
n   <- length(fake.data) # sample size
p   <- length(fake.data[fake.data == 0])/n # est. prob of zero value
F   <- nonexceeds(sig6=TRUE); F <- sort(c(F,p)); qF <- qnorm(F)
PP  <- pp(fake.data) # compute plotting positions of sim. sample
plot(qnorm(PP), fake.data, ylim=c(0,4000), xaxt="n", xlab="") # plot the sample
add.lmomco.axis(las=2, tcl=0.5, twoside=TRUE, side.type="SNV", otherside.type="NA")
lines(qF,quagpa(F,gpa.all)) # the parent (without zeros)
lines(qF,z.par2qua(F,p,gpa.nzo),lwd=3) # fitted model with zero conditional
par(mgp=c(3,1,0)) # restore defaults
}
\keyword{quantile function}
\keyword{blipping distributions}

